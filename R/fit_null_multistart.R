#' Multi-start null model fitting for the robust score test
#'
#' Fits the null model \eqn{B_{k,j} = g(B_k)} by polishing several diverse
#' starting points with an exact reduced damped-Newton solver
#' (\code{fit_null_local()}) and returning the feasible stationary point with
#' the highest profiled log-likelihood.
#'
#' For concave constraints (arithmetic mean, subset mean, single-category
#' reference), the problem has a unique global optimum and a single exact local
#' solve suffices. For the non-convex pseudo-Huber smoothed median (all
#' categories or a subset), local optima correspond to different partitions
#' ("sheets") of the remaining categories into those above and below
#' \eqn{B_{k,j}}. \code{fit_null_multistart()} searches across:
#' \enumerate{
#'   \item the supplied starting estimate \code{B} (with \eqn{I^{-1}}-metric
#'     projection onto the constraint),
#'   \item the zero matrix \code{B = 0} and a flat constrained row
#'     (\code{B[k_constr, ] = 0}, other rows from \code{B}),
#'   \item smoothing continuation in the pseudo-Huber parameter \eqn{d'} from
#'     twice the row spread down to \eqn{d} by factors of \code{d_ratio},
#'   \item up to \code{sheets} distinct sheet candidates from a separable
#'     Poisson surrogate solved by dynamic programming (for constraints with
#'     at most 400 arguments), and
#'   \item a local exchange search around the best solution, reflecting pairs
#'     of categories nearest below and above \eqn{B_{k,j}} across \eqn{B_{k,j}}.
#' }
#'
#' @param B Starting coefficient matrix (\code{p x J}), typically the full-model
#'   estimate.
#' @param Y Observation matrix (\code{n x J}) with augmentations held fixed.
#' @param X Design matrix (\code{n x p}).
#' @param k_constr Row index of \code{B} to constrain (1-based).
#' @param j_constr Column index of \code{B} to constrain (1-based).
#' @param j_ref Column index used for the convenience shift \code{B[, j_ref] = 0}.
#' @param constraint_fn Constraint function \eqn{g}, or a list of length \code{p}
#'   of constraint functions.
#' @param constraint_grad_fn Gradient function \eqn{\nabla g}, or a list of
#'   length \code{p} of gradient functions.
#' @param tol_kkt Stopping tolerance on the scaled reduced Gauss-Newton
#'   decrement. Default is \code{1e-9}.
#' @param c1 Armijo line-search constant. Default is \code{1e-4}.
#' @param maxit Maximum Newton iterations per local solve. Default is \code{500}.
#' @param max_step Maximum sup-norm step size per iteration. Default is \code{2}.
#' @param max_abs_B Maximum allowed absolute value of any entry of \code{B}.
#'   Default is \code{250}.
#' @param max_backtrack Maximum line-search halvings per iteration. Default is
#'   \code{40}.
#' @param metric_projection Logical; project each start onto the constraint in
#'   the \eqn{I^{-1}} metric before running Newton? Default is \code{TRUE}.
#' @param zero_start Logical; include \code{B = 0} as a start? Default is
#'   \code{TRUE}.
#' @param flat_start Logical; include the start with row \code{k_constr} set to
#'   zero and other rows from \code{B}? Default is \code{TRUE}.
#' @param d_ratio Geometric reduction factor (> 1) for smoothing continuation in
#'   the pseudo-Huber parameter \code{d}. Set \code{<= 1} to disable. Default
#'   is \code{2}.
#' @param continuation_kkt KKT tolerance for intermediate continuation solves.
#'   Default is \code{1e-6}.
#' @param sheets Number of distinct DP sheet candidates to polish (set \code{0}
#'   to disable; automatically skipped when the constraint has more than 400
#'   arguments). Default is \code{8}.
#' @param search_width Number of nearest categories below and above
#'   \eqn{B_{k,j}} considered in each round of the exchange search (set \code{0}
#'   to disable). Default is \code{3}.
#' @param search_rounds Maximum number of improving exchange-search rounds.
#'   Default is \code{5}.
#' @param verbose Logical; print progress messages? Default is \code{FALSE}.
#' @param trackB Logical; retained for interface compatibility with
#'   \code{score_test()} (unused; \code{Bs} is returned as \code{NULL}).
#' @param solver Linear solver backend: \code{"auto"} (default; dense Cholesky
#'   for \eqn{p(J-1) \le 250} or \eqn{p(J-1) \le n_{\mathrm{distinct}}},
#'   structured Woodbury + PCG otherwise), \code{"dense"}, or
#'   \code{"structured"}.
#'
#' @return A list compatible with \code{score_test()} containing elements
#'   \code{B}, \code{k_constr}, \code{j_constr}, \code{niter},
#'   \code{converged}, \code{gap}, \code{loglik}, \code{kkt},
#'   \code{local_solves}, \code{newton_steps}, \code{Bs}, and \code{it_df}.
#'
#' @keywords internal
#' @noRd
fit_null_multistart <- function(
  B,
  Y,
  X,
  k_constr,
  j_constr,
  j_ref,
  constraint_fn,
  constraint_grad_fn,
  tol_kkt = 1e-9,
  c1 = 1e-4,
  maxit = 500L,
  max_step = 2,
  max_abs_B = 250,
  max_backtrack = 40L,
  metric_projection = TRUE,
  zero_start = TRUE,
  flat_start = TRUE,
  d_ratio = 2,
  continuation_kkt = 1e-6,
  sheets = 8L,
  search_width = 3L,
  search_rounds = 5L,
  verbose = FALSE,
  trackB = FALSE,
  solver = c("auto", "dense", "structured")
) {
  solver <- match.arg(solver)
  if (is.list(constraint_fn)) {
    constraint_fn <- constraint_fn[[k_constr]]
    constraint_grad_fn <- constraint_grad_fn[[k_constr]]
  }

  Y <- as.matrix(Y)
  X <- as.matrix(X)
  B <- as.matrix(B)
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(B) <- "double"

  p <- ncol(X)
  J <- ncol(Y)

  g_spec <- classify_null_constraint(constraint_fn, constraint_grad_fn, J, j_constr)
  if (g_spec$kind == "other") {
    stop("Unsupported constraint type for multistart null fit; use augmented Lagrangian.")
  }

  # Vacuous constraint (e.g. reference constraint tested at its own column): H = 0
  if (length(g_spec$args) == 0L) {
    B_out <- B - B[, j_ref]
    col_data <- null_collapse_rows(Y, X)
    ev0 <- null_eval_ll(list(X = col_data$X, Y = col_data$Y, m = col_data$m, ll_const = col_data$ll_const), B_out)
    return(list(
      B = B_out,
      k_constr = k_constr,
      j_constr = j_constr,
      niter = 0L,
      converged = TRUE,
      gap = 0,
      loglik = ev0$ll,
      kkt = 0,
      local_solves = 0L,
      newton_steps = 0L,
      Bs = NULL,
      it_df = NULL
    ))
  }

  # Working reference column r != j_constr (same rule as fit_null_symmetric)
  if (j_ref == j_constr) {
    r_idx <- which.max(colSums(Y > 1)[-j_ref])
    r <- r_idx + as.integer(r_idx >= j_ref)
  } else {
    r <- j_ref
  }

  ctx <- null_make_ctx(
    Y_orig = Y,
    X_orig = X,
    k_constr = k_constr,
    j_constr = j_constr,
    r = r,
    g_spec = g_spec,
    solver = solver
  )

  run_local <- function(B_init, ctx_use = ctx, tol_use = tol_kkt) {
    fit_null_local(
      ctx = ctx_use,
      B_start = B_init,
      j_ref = j_ref,
      tol_kkt = tol_use,
      c1 = c1,
      maxit = as.integer(maxit),
      max_step = max_step,
      max_abs_B = max_abs_B,
      max_backtrack = as.integer(max_backtrack),
      metric_projection = metric_projection
    )
  }

  # Primary local solve from supplied B (full-model estimate)
  best <- run_local(B)

  # Concave constraints (mean, subset mean, reference) have a unique global optimum
  if (g_spec$kind != "pseudohuber") {
    return(list(
      B = best$B,
      k_constr = k_constr,
      j_constr = j_constr,
      niter = best$niter,
      converged = best$converged,
      gap = best$gap,
      loglik = best$loglik,
      kkt = best$kkt,
      local_solves = best$local_solves,
      newton_steps = best$newton_steps,
      Bs = NULL,
      it_df = NULL
    ))
  }

  is_better <- function(cand, cur) {
    if (cand$converged && !cur$converged) {
      return(TRUE)
    }
    if (cand$converged == cur$converged &&
        is.finite(cand$loglik) &&
        (!is.finite(cur$loglik) || cand$loglik > cur$loglik + 1e-10 * abs(cur$loglik))) {
      return(TRUE)
    }
    FALSE
  }

  consider <- function(cand, cur) {
    total_solves <- cur$local_solves + cand$local_solves
    out <- if (is_better(cand, cur)) cand else cur
    out$local_solves <- total_solves
    out
  }

  # (1) Zero start (as fit_null_discrete) and flat-row start
  if (isTRUE(zero_start)) {
    best <- consider(run_local(matrix(0, nrow = p, ncol = J)), best)
  }
  if (isTRUE(flat_start)) {
    B_flat <- B
    B_flat[k_constr, ] <- 0
    best <- consider(run_local(B_flat), best)
  }

  # (2) Smoothing continuation in pseudo-Huber d from 2 * spread down to d
  if (is.finite(d_ratio) && d_ratio > 1) {
    d_target <- g_spec$d
    vals <- B[k_constr, g_spec$args]
    spread <- max(vals) - min(vals)
    dd <- 2 * spread
    if (is.finite(dd) && dd > d_target) {
      cur_B <- B
      c_solves <- 0L
      while (dd > d_target * (1 + 1e-7)) {
        ctx_d <- ctx
        ctx_d$g_spec$d <- dd
        f_mid <- run_local(cur_B, ctx_use = ctx_d, tol_use = continuation_kkt)
        c_solves <- c_solves + 1L
        cur_B <- f_mid$B
        dd <- dd / d_ratio
      }
      cand <- run_local(cur_B)
      cand$local_solves <- cand$local_solves + c_solves
      best <- consider(cand, best)
    }
  }

  # (3) Sheet candidates from the separable surrogate (dynamic programming)
  if (as.integer(sheets) > 0L) {
    starts <- null_sheet_starts(ctx, B, count = as.integer(sheets))
    for (st_mat in starts) {
      best <- consider(run_local(st_mat), best)
    }
  }

  # (4) Exchange search around the winning sheet
  sw <- as.integer(search_width)
  sr <- as.integer(search_rounds)
  if (sw > 0L && sr > 0L) {
    args <- g_spec$args
    for (rd in seq_len(sr)) {
      bj <- best$B[k_constr, j_constr]
      b_args <- best$B[k_constr, args]
      below <- args[b_args < bj]
      above <- args[b_args >= bj]
      if (length(below) == 0L || length(above) == 0L) {
        break
      }
      below <- below[order(bj - best$B[k_constr, below])]
      above <- above[order(best$B[k_constr, above] - bj)]
      below <- utils::head(below, sw)
      above <- utils::head(above, sw)

      improved <- NULL
      r_solves <- 0L
      for (u in below) {
        for (v in above) {
          st_mat <- best$B
          st_mat[k_constr, u] <- 2 * bj - st_mat[k_constr, u]
          st_mat[k_constr, v] <- 2 * bj - st_mat[k_constr, v]
          cand <- run_local(st_mat)
          r_solves <- r_solves + 1L
          ref_fit <- if (is.null(improved)) best else improved
          if (is_better(cand, ref_fit)) {
            improved <- cand
          }
        }
      }
      best$local_solves <- best$local_solves + r_solves
      if (!is.null(improved)) {
        improved$local_solves <- best$local_solves
        best <- improved
      } else {
        break
      }
    }
  }

  if (isTRUE(verbose)) {
    message(
      "Multistart null fit (k = ", k_constr, ", j = ", j_constr, "): ",
      best$local_solves, " local solves, ll = ", round(best$loglik, 3),
      ", kkt = ", signif(best$kkt, 3), ", gap = ", signif(best$gap, 3)
    )
  }

  list(
    B = best$B,
    k_constr = k_constr,
    j_constr = j_constr,
    niter = best$niter,
    converged = best$converged,
    gap = best$gap,
    loglik = best$loglik,
    kkt = best$kkt,
    local_solves = best$local_solves,
    newton_steps = best$newton_steps,
    Bs = NULL,
    it_df = NULL
  )
}

# Candidate starts ("sheets") for the non-convex pseudo-Huber constraint from a
# separable Poisson surrogate solved by vectorised dynamic programming.
null_sheet_starts <- function(ctx, B_start, count = 8L) {
  if (ctx$g_spec$kind != "pseudohuber" || count <= 0L) {
    return(list())
  }
  d <- ctx$g_spec$d
  args <- ctx$g_spec$args
  ma <- length(args)
  if (ma < 2L || ma > 400L) {
    return(list())
  }

  ev0 <- null_eval_ll(ctx, B_start)
  if (is.null(ev0$P)) {
    return(list())
  }
  mu <- ctx$row_n * ev0$P

  # Aggregate rows with identical x_{ik} when there are few distinct covariate levels
  xk <- ctx$X[, ctx$k]
  levels <- sort(unique(xk))
  if (length(levels) <= 32L && length(levels) < ctx$m) {
    grp <- match(xk, levels)
    A_all <- unname(rowsum(mu, grp, reorder = TRUE))
    B_all <- unname(rowsum(ctx$Y, grp, reorder = TRUE))
    v_vec <- levels
  } else {
    A_all <- mu
    B_all <- ctx$Y
    v_vec <- xk
  }

  pcost_vec <- function(a_col, b_col, eps_vec) {
    ve <- outer(v_vec, eps_vec)
    em1 <- expm1(ve)
    term1 <- a_col * em1
    term1[a_col == 0] <- 0
    colSums(term1 - b_col * ve)
  }

  psi <- function(u) {
    r <- u / d
    r / sqrt(1 + r * r)
  }

  H <- 0.1
  max_shift <- as.integer(round(1 / H))
  smax <- ma * max_shift + 2L
  width <- 2L * smax + 1L
  s0 <- smax + 1L
  zone <- c(0, 0.5 * d, -0.5 * d, d, -d, 2 * d, -2 * d)

  bl <- B_start[ctx$k, args]
  bj <- B_start[ctx$k, ctx$j]
  m0 <- pseudohuber_median_exact(bl, d = d)
  if (!is.finite(m0) || !is.finite(bj)) {
    return(list())
  }

  n_m <- 24L
  if (abs(m0 - bj) < 1e-14) {
    m_grid <- bj
  } else {
    ext <- 0.5 * (m0 - bj)
    m_grid <- (bj - ext) + (m0 - bj + 2 * ext) * (0:n_m) / n_m
  }

  A_args <- A_all[, args, drop = FALSE]
  B_args <- B_all[, args, drop = FALSE]
  a_j <- A_all[, ctx$j]
  b_j <- B_all[, ctx$j]

  results <- list()
  f_buf <- rep(Inf, width)
  f2_buf <- rep(Inf, width)
  choice <- matrix(0L, nrow = ma, ncol = width)

  for (m_val in m_grid) {
    f <- f_buf
    f2 <- f2_buf
    f[s0] <- 0
    lo <- s0
    hi <- s0

    for (a in seq_len(ma)) {
      uh <- bl[a] - m_val
      u_cands <- c(uh, zone)
      c_cands <- pcost_vec(A_args[, a], B_args[, a], m_val + u_cands - bl[a])
      si_cands <- as.integer(round(psi(u_cands) / H))

      src <- lo:hi
      f_src <- f[src]
      finite_mask <- is.finite(c_cands)
      if (!any(finite_mask)) {
        lo <- 1L
        hi <- 0L
        break
      }
      new_lo <- max(1L, lo + min(si_cands[finite_mask]))
      new_hi <- min(width, hi + max(si_cands[finite_mask]))

      for (ci in which(finite_mask)) {
        dst <- src + si_cands[ci]
        cand_val <- f_src + c_cands[ci]
        hit <- which(cand_val < f2[dst])
        if (length(hit) > 0L) {
          d_hit <- dst[hit]
          f2[d_hit] <- cand_val[hit]
          choice[a, d_hit] <- ci
        }
      }

      f[lo:hi] <- Inf
      tmp <- f
      f <- f2
      f2 <- tmp
      lo <- new_lo
      hi <- new_hi
    }

    if (lo > hi) {
      next
    }

    cand_s <- max(1L, s0 - 1L):min(width, s0 + 1L)
    bs <- cand_s[which.min(f[cand_s])]
    if (!is.finite(f[bs])) {
      next
    }
    cost_j <- pcost_vec(a_j, b_j, m_val - bj)[1L]
    total_cost <- f[bs] + cost_j
    if (!is.finite(total_cost)) {
      next
    }

    # Backtrack optimal offset per argument
    u_out <- numeric(ma)
    s_cur <- bs
    for (a in seq.int(ma, 1L, by = -1L)) {
      ci <- choice[a, s_cur]
      uc <- if (ci == 1L) bl[a] - m_val else zone[ci - 1L]
      u_out[a] <- uc
      s_cur <- s_cur - as.integer(round(psi(uc) / H))
    }

    results[[length(results) + 1L]] <- list(cost = total_cost, m = m_val, u = u_out)
  }

  if (length(results) == 0L) {
    return(list())
  }

  ord <- order(vapply(results, `[[`, numeric(1), "cost"))
  out <- list()
  seen <- character(0)
  for (idx in ord) {
    res <- results[[idx]]
    sig <- paste0(as.integer(res$u > 0), collapse = "")
    if (sig %in% seen) {
      next
    }
    seen <- c(seen, sig)
    st_mat <- B_start
    st_mat[ctx$k, ctx$j] <- res$m
    st_mat[ctx$k, args] <- res$m + res$u
    out[[length(out) + 1L]] <- st_mat
    if (length(out) >= count) {
      break
    }
  }
  out
}
