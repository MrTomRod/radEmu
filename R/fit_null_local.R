#' Exact reduced damped-Newton local solver for the null model (internal)
#'
#' Maximises the profiled Poisson log-likelihood \eqn{\ell(B)} on fixed
#' (augmented) counts \code{Y} subject to \eqn{B_{k,j} = g(B_{k,-j})} by
#' eliminating \eqn{B_{k,j}} via the constraint and running damped Newton on the
#' unconstrained reduced coordinates \eqn{\theta = B_{-r} \setminus \{B_{k,j}\}}.
#'
#' Identical covariate rows in \code{X} are collapsed into exact sufficient
#' statistics (group count totals) before optimisation, preserving the objective
#' up to a constant and leaving the gradient, Fisher information, and optimum
#' unchanged. Small and moderate systems (\eqn{q = p(J - 1) \le 250} or
#' \eqn{q \le n_{\mathrm{distinct}}}) use dense Cholesky factorisation of the
#' exact reduced Hessian \eqn{J^\top I_{-r} J - \nu E_g}; larger systems use
#' structured Woodbury solves \eqn{I_{-r}^{-1} = D^{-1} + D^{-1} U K^{-1} U^\top D^{-1}}
#' with tangent-space preconditioned conjugate gradients for the exact Newton
#' step and Gauss-Newton fallback.
#'
#' @keywords internal
#' @noRd

# Exact sufficient-statistics collapse of duplicate rows of X
null_collapse_rows <- function(Y, X) {
  n <- nrow(Y)
  row_n_orig <- rowSums(Y)
  if (any(!is.finite(row_n_orig) | row_n_orig <= 0)) {
    stop("All rows of Y must have strictly positive finite totals.")
  }
  ll_const_orig <- sum(row_n_orig * (log(row_n_orig) - 1))

  key <- apply(X, 1L, function(r) paste(sprintf("%.17g", r), collapse = "\t"))
  uniq_idx <- which(!duplicated(key))
  if (length(uniq_idx) == n) {
    return(list(
      Y = Y,
      X = X,
      row_n = row_n_orig,
      ll_const = ll_const_orig,
      m = n
    ))
  }

  grp <- match(key, key[uniq_idx])
  X_c <- unname(X[uniq_idx, , drop = FALSE])
  Y_c <- unname(rowsum(Y, grp, reorder = TRUE))
  storage.mode(Y_c) <- "double"
  row_n_c <- rowSums(Y_c)
  list(
    Y = Y_c,
    X = X_c,
    row_n = row_n_c,
    ll_const = ll_const_orig,
    m = nrow(X_c)
  )
}

# Stable evaluation of profiled Poisson log-likelihood and softmax probabilities
null_eval_ll <- function(ctx, B) {
  eta <- ctx$X %*% B
  if (any(!is.finite(eta))) {
    return(list(ll = -Inf, P = NULL))
  }
  mx <- eta[cbind(seq_len(ctx$m), max.col(eta, ties.method = "first"))]
  ex <- exp(eta - mx)
  s <- rowSums(ex)
  if (any(!is.finite(s) | s <= 0)) {
    return(list(ll = -Inf, P = NULL))
  }
  lse <- mx + log(s)
  P <- ex / s
  ll <- sum(ctx$Y * (eta - lse)) + ctx$ll_const
  list(ll = ll, P = P)
}

# Row constraint value and derivatives on row k of B
null_g_value <- function(row_k, g_spec) {
  x <- row_k[g_spec$args]
  if (g_spec$kind == "pseudohuber") {
    pseudohuber_median_exact(x, d = g_spec$d)
  } else if (g_spec$kind == "mean") {
    mean(x)
  } else {
    x[1L]
  }
}

null_g_deriv <- function(row_k, g_spec) {
  x <- row_k[g_spec$args]
  if (g_spec$kind == "pseudohuber") {
    pseudohuber_derivatives_exact(x, d = g_spec$d)
  } else if (g_spec$kind == "mean") {
    ma <- length(x)
    list(
      value = mean(x),
      grad = rep(1 / ma, ma),
      kp = rep(0, ma),
      s_sum = 1,
      k_sum = 0
    )
  } else {
    list(value = x[1L], grad = 1, kp = 0, s_sum = 1, k_sum = 0)
  }
}

# Dense Cholesky solve with R's ridge ladder fallback
null_chol_solve <- function(M, rhs) {
  R <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(R)) {
    return(NULL)
  }
  x <- backsolve(R, backsolve(R, rhs, transpose = TRUE))
  if (all(is.finite(x))) x else NULL
}

null_chol_solve_ladder <- function(M, rhs) {
  x <- null_chol_solve(M, rhs)
  if (!is.null(x)) {
    return(list(x = x, lambda = 0))
  }
  n <- nrow(M)
  scale <- max(abs(diag(M)))
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  lambda <- 1e-4
  while (lambda <= 1e6) {
    M_reg <- M
    diag(M_reg) <- diag(M_reg) + lambda * scale
    x <- null_chol_solve(M_reg, rhs)
    if (!is.null(x)) {
      return(list(x = x, lambda = lambda))
    }
    lambda <- lambda * 10
  }
  NULL
}

# Dense information I_{-r} (q x q), vectorised via BLAS
null_dense_info <- function(ctx, P) {
  P_nr <- P[, ctx$cols, drop = FALSE]
  mu_nr <- ctx$row_n * P_nr
  W_nr <- ctx$sqrt_row_n * P_nr

  # U (q x m): column i is sqrt(N_i) * (p_{i,-r} (x) x_i)
  U <- t(W_nr)[ctx$u_col_rep, , drop = FALSE] * ctx$tX[ctx$u_x_rep, , drop = FALSE]
  I_mat <- -tcrossprod(U)
  # Diagonal blocks I_{ll} = X^T diag(mu_l (1 - p_l)) X in one BLAS call
  I_mat[ctx$diag_block_idx] <- crossprod(ctx$XX, mu_nr * (1 - P_nr))
  I_mat
}

# Structured Woodbury factorisation of I_{-r} + ridge * Id = (D + ridge * Id) - U U^T
null_woodbury_factor <- function(ctx, P, ridge = 0) {
  p <- ctx$p
  m <- ctx$m
  Jm1 <- ctx$Jm1
  P_nr <- P[, ctx$cols, drop = FALSE]
  mu_nr <- ctx$row_n * P_nr
  W_nr <- ctx$sqrt_row_n * P_nr

  D_blocks <- crossprod(ctx$XX, mu_nr) # p^2 x Jm1
  if (ridge > 0) {
    D_blocks[ctx$p_diag_idx, ] <- D_blocks[ctx$p_diag_idx, ] + ridge
  }

  if (p == 1L) {
    d11 <- as.numeric(D_blocks)
    if (any(!is.finite(d11) | d11 <= 0)) {
      return(NULL)
    }
    r11 <- sqrt(d11)
    u1 <- (1 / r11) * t(W_nr * ctx$X[, 1L])
    K <- diag(m) - crossprod(u1)
    R_K <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(R_K)) {
      return(NULL)
    }
    dinv <- function(v_mat) {
      matrix(as.numeric(v_mat) / d11, nrow = 1L, ncol = Jm1)
    }
    dmul <- function(u_mat) {
      matrix(d11 * as.numeric(u_mat), nrow = 1L, ncol = Jm1)
    }
  } else if (p == 2L) {
    d11 <- D_blocks[1L, ]
    d12 <- D_blocks[2L, ]
    d22 <- D_blocks[4L, ]
    if (any(!is.finite(d11) | d11 <= 0)) {
      return(NULL)
    }
    r11 <- sqrt(d11)
    r12 <- d12 / r11
    s22 <- d22 - r12 * r12
    if (any(!is.finite(s22) | s22 <= 0)) {
      return(NULL)
    }
    r22 <- sqrt(s22)
    tW <- t(W_nr)
    x1 <- ctx$X[, 1L]
    x2 <- ctx$X[, 2L]
    u1 <- outer(1 / r11, x1) * tW
    u2 <- ((outer(rep(1, Jm1), x2) - outer(r12 / r11, x1)) / r22) * tW
    K <- diag(m) - crossprod(u1) - crossprod(u2)
    R_K <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(R_K)) {
      return(NULL)
    }
    dinv <- function(v_mat) {
      y1 <- v_mat[1L, ] / r11
      y2 <- (v_mat[2L, ] - r12 * y1) / r22
      x2_v <- y2 / r22
      x1_v <- (y1 - r12 * x2_v) / r11
      rbind(x1_v, x2_v)
    }
    dmul <- function(u_mat) {
      rbind(
        d11 * u_mat[1L, ] + d12 * u_mat[2L, ],
        d12 * u_mat[1L, ] + d22 * u_mat[2L, ]
      )
    }
  } else {
    R_list <- vector("list", Jm1)
    t_U_tilde <- matrix(0, nrow = ctx$q, ncol = m)
    for (idx in seq_len(Jm1)) {
      Rj <- tryCatch(chol(matrix(D_blocks[, idx], p, p)), error = function(e) NULL)
      if (is.null(Rj)) {
        return(NULL)
      }
      R_list[[idx]] <- Rj
      Xi_t <- backsolve(Rj, ctx$tX, transpose = TRUE)
      rows <- ((idx - 1L) * p + 1L):(idx * p)
      t_U_tilde[rows, ] <- rep(W_nr[, idx], each = p) * Xi_t
    }
    K <- diag(m) - crossprod(t_U_tilde)
    R_K <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(R_K)) {
      return(NULL)
    }
    dinv <- function(v_mat) {
      out <- matrix(0, nrow = p, ncol = Jm1)
      for (idx in seq_len(Jm1)) {
        Rj <- R_list[[idx]]
        out[, idx] <- backsolve(Rj, backsolve(Rj, v_mat[, idx], transpose = TRUE))
      }
      out
    }
    dmul <- function(u_mat) {
      out <- matrix(0, nrow = p, ncol = Jm1)
      for (idx in seq_len(Jm1)) {
        out[, idx] <- matrix(D_blocks[, idx], p, p) %*% u_mat[, idx]
      }
      out
    }
  }

  solve_vec <- function(v_vec) {
    y_mat <- dinv(matrix(v_vec, nrow = p, ncol = Jm1))
    t_vec <- rowSums(W_nr * (ctx$X %*% y_mat))
    tp <- backsolve(R_K, backsolve(R_K, t_vec, transpose = TRUE))
    g_mat <- crossprod(ctx$X, tp * W_nr)
    as.vector(y_mat + dinv(g_mat))
  }

  apply_I <- function(u_vec) {
    u_mat <- matrix(u_vec, nrow = p, ncol = Jm1)
    Du <- dmul(u_mat)
    t_vec <- rowSums(W_nr * (ctx$X %*% u_mat))
    Ut <- crossprod(ctx$X, t_vec * W_nr)
    as.vector(Du - Ut)
  }

  list(solve = solve_vec, apply_I = apply_I, ridge = ridge)
}

null_woodbury_ladder <- function(ctx, P) {
  fac <- null_woodbury_factor(ctx, P, ridge = 0)
  if (!is.null(fac)) {
    return(list(fac = fac, lambda = 0))
  }
  P_nr <- P[, ctx$cols, drop = FALSE]
  mu_nr <- ctx$row_n * P_nr
  diag_I <- crossprod(ctx$X * ctx$X, mu_nr * (1 - P_nr))
  scale <- max(abs(diag_I))
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  lambda <- 1e-4
  while (lambda <= 1e6) {
    fac <- null_woodbury_factor(ctx, P, ridge = lambda * scale)
    if (!is.null(fac)) {
      return(list(fac = fac, lambda = lambda))
    }
    lambda <- lambda * 10
  }
  NULL
}

# Exact Newton direction in tangent space {u : a^T u = 0} via preconditioned CG (structured path)
null_newton_pcg <- function(ctx, fac, y2, aty2, a_vec, nu, gth_full, gd) {
  keep <- ctx$g_args_keep
  q_idx <- ctx$g_args_qidx
  gv_sub <- gd$grad[keep]
  kp_sub <- gd$kp[keep]
  s_sum <- gd$s_sum
  k_sum <- gd$k_sum

  wmul <- function(u) {
    out <- fac$apply_I(u)
    if (length(q_idx) > 0L) {
      u_sub <- u[q_idx]
      gu <- sum(gv_sub * u_sub)
      ku <- sum(kp_sub * u_sub)
      e <- kp_sub * u_sub - gv_sub * ku - kp_sub * gu + k_sum * gv_sub * gu
      out[q_idx] <- out[q_idx] - (nu / s_sum) * e
    }
    out
  }

  pre <- function(rv) {
    y1 <- fac$solve(rv)
    y1 - (sum(a_vec * y1) / aty2) * y2
  }

  x <- numeric(ctx$q)
  res <- gth_full
  z <- pre(res)
  pdir <- z
  rz <- sum(res * z)
  rz0 <- rz
  if (!is.finite(rz0) || rz0 <= 0) {
    return(NULL)
  }
  stop_tol <- max(min(1e-10, rz0), 1e-26) * rz0

  for (it in seq_len(100L)) {
    ap <- wmul(pdir)
    curv <- sum(pdir * ap)
    if (!is.finite(curv) || curv <= 0) {
      return(NULL)
    }
    alpha <- rz / curv
    x <- x + alpha * pdir
    res <- res - alpha * ap
    z <- pre(res)
    rz_new <- sum(res * z)
    if (!is.finite(rz_new)) {
      return(NULL)
    }
    if (rz_new <= 0 || rz_new <= stop_tol) {
      break
    }
    beta <- rz_new / rz
    rz <- rz_new
    pdir <- z + beta * pdir
  }

  if (all(is.finite(x))) x else NULL
}

# Build context shared across all local solves for a given (Y, X, k, j, r, g_spec)
null_make_ctx <- function(Y_orig, X_orig, k_constr, j_constr, r, g_spec,
                          solver = c("auto", "dense", "structured")) {
  solver <- match.arg(solver)
  col_data <- null_collapse_rows(Y_orig, X_orig)
  p <- ncol(X_orig)
  J <- ncol(Y_orig)
  Jm1 <- J - 1L
  q <- p * Jm1
  m <- col_data$m

  cols <- seq_len(J)[-r]
  cpos <- integer(J)
  cpos[cols] <- seq_len(Jm1)
  cpos[r] <- NA_integer_

  c0 <- (cpos[j_constr] - 1L) * p + k_constr
  umap <- seq_len(q)[-c0]
  tpos <- integer(q)
  tpos[umap] <- seq_len(q - 1L)

  # Precompute X_ia * X_ib matrix (m x p^2) for single-BLAS D_l computation
  X <- col_data$X
  XX <- matrix(0, nrow = m, ncol = p * p)
  for (b in seq_len(p)) {
    for (a in seq_len(p)) {
      XX[, (b - 1L) * p + a] <- X[, a] * X[, b]
    }
  }

  use_dense <- if (solver == "dense") {
    TRUE
  } else if (solver == "structured") {
    FALSE
  } else {
    q <= 250L || (q <= 600L && q <= m)
  }

  # Block-diagonal indices for dense I_{-r}
  diag_block_idx <- NULL
  u_col_rep <- NULL
  u_x_rep <- NULL
  if (use_dense) {
    a_rep <- rep(seq_len(p), times = p * Jm1)
    b_rep <- rep(rep(seq_len(p), each = p), times = Jm1)
    lp_rep <- rep(seq_len(Jm1), each = p * p)
    offset <- (lp_rep - 1L) * p
    diag_block_idx <- (offset + a_rep) + (offset + b_rep - 1L) * q
    u_col_rep <- rep(seq_len(Jm1), each = p)
    u_x_rep <- rep(seq_len(p), times = Jm1)
  }

  g_args_keep <- which(g_spec$args != r)
  g_args_nr <- g_spec$args[g_args_keep]
  g_args_qidx <- if (length(g_args_nr) > 0L) {
    (cpos[g_args_nr] - 1L) * p + k_constr
  } else {
    integer(0)
  }
  g_args_tidx <- if (length(g_args_qidx) > 0L) tpos[g_args_qidx] else integer(0)

  list(
    Y = col_data$Y,
    X = X,
    tX = t(X),
    XX = XX,
    row_n = col_data$row_n,
    sqrt_row_n = sqrt(col_data$row_n),
    ll_const = col_data$ll_const,
    m = m,
    p = p,
    J = J,
    Jm1 = Jm1,
    q = q,
    k = k_constr,
    j = j_constr,
    r = r,
    cols = cols,
    cpos = cpos,
    c0 = c0,
    umap = umap,
    use_dense = use_dense,
    diag_block_idx = diag_block_idx,
    u_col_rep = u_col_rep,
    u_x_rep = u_x_rep,
    p_diag_idx = seq(1L, p * p, by = p + 1L),
    g_spec = g_spec,
    g_args_keep = g_args_keep,
    g_args_qidx = g_args_qidx,
    g_args_tidx = g_args_tidx
  )
}

# Build constraint gradient a = grad_{B_{-r}} c(B) in q-coordinates
null_constraint_grad_q <- function(ctx, grad_g) {
  a <- numeric(ctx$q)
  a[ctx$c0] <- 1
  if (length(ctx$g_args_qidx) > 0L) {
    a[ctx$g_args_qidx] <- -grad_g[ctx$g_args_keep]
  }
  a
}

# Move B to the feasible set along the I^{-1}-metric steepest direction
null_metric_projection <- function(ctx, B, max_abs_B) {
  for (it in seq_len(50L)) {
    cv <- B[ctx$k, ctx$j] - null_g_value(B[ctx$k, ], ctx$g_spec)
    if (!is.finite(cv) || abs(cv) < 1e-12) {
      break
    }
    ev <- null_eval_ll(ctx, B)
    if (is.null(ev$P)) {
      return(B)
    }
    gd <- null_g_deriv(B[ctx$k, ], ctx$g_spec)
    a <- null_constraint_grad_q(ctx, gd$grad)

    if (ctx$use_dense) {
      I_mat <- null_dense_info(ctx, ev$P)
      sol <- null_chol_solve_ladder(I_mat, a)
      if (is.null(sol)) return(B)
      y <- sol$x
    } else {
      fac_res <- null_woodbury_ladder(ctx, ev$P)
      if (is.null(fac_res)) return(B)
      y <- fac_res$fac$solve(a)
    }

    aty <- sum(a * y)
    if (!is.finite(aty) || aty <= 0 || any(!is.finite(y))) {
      return(B)
    }
    f_step <- -cv / aty
    B_new <- B
    B_new[, ctx$cols] <- B_new[, ctx$cols] + f_step * matrix(y, nrow = ctx$p, ncol = ctx$Jm1)
    if (max(abs(B_new)) > max_abs_B) {
      return(B)
    }
    B <- B_new
  }
  B
}

#' Run one local damped-Newton null fit from a single starting matrix
#'
#' @param ctx Context created by \code{null_make_ctx()}.
#' @param B_start Starting coefficient matrix (\code{p x J}).
#' @param j_ref Caller's reference category index (restored on output).
#' @param tol_kkt KKT convergence tolerance (scaled Gauss-Newton decrement).
#' @param c1 Armijo line search constant.
#' @param maxit Maximum Newton iterations.
#' @param max_step Maximum sup-norm step size in free coordinates.
#' @param max_abs_B Maximum allowed absolute value of entries in \code{B}.
#' @param max_backtrack Maximum line-search halvings.
#' @param metric_projection Logical; project \code{B_start} onto the constraint
#'   in the \eqn{I^{-1}} metric before iterating.
#' @return List with \code{B}, \code{loglik}, \code{niter}, \code{kkt},
#'   \code{gap}, \code{converged}, \code{newton_steps}, \code{local_solves}.
#' @keywords internal
#' @noRd
fit_null_local <- function(ctx,
                           B_start,
                           j_ref,
                           tol_kkt = 1e-9,
                           c1 = 1e-4,
                           maxit = 500L,
                           max_step = 2,
                           max_abs_B = 250,
                           max_backtrack = 40L,
                           metric_projection = TRUE) {
  p <- ctx$p
  J <- ctx$J
  k <- ctx$k
  j <- ctx$j
  r <- ctx$r
  q <- ctx$q
  nth <- q - 1L
  is_ph <- (ctx$g_spec$kind == "pseudohuber")

  # Shift to working reference B[, r] = 0
  B <- B_start - B_start[, r]
  if (metric_projection) {
    B <- null_metric_projection(ctx, B, max_abs_B = max_abs_B)
  }
  B[k, j] <- null_g_value(B[k, ], ctx$g_spec)

  iter <- 0L
  newton_steps <- 0L
  kkt <- Inf

  ev <- null_eval_ll(ctx, B)
  if (!is.finite(ev$ll) || is.null(ev$P)) {
    # Fallback to zero-row start if supplied start overflowed
    B <- matrix(0, nrow = p, ncol = J)
    ev <- null_eval_ll(ctx, B)
  }

  repeat {
    mu <- ctx$row_n * ev$P
    s_mat <- crossprod(ctx$X, ctx$Y - mu) # p x J
    s_q <- as.vector(s_mat[, ctx$cols, drop = FALSE])
    nu <- s_q[ctx$c0]

    if (nth == 0L) {
      kkt <- 0
      break
    }

    gd <- null_g_deriv(B[k, ], ctx$g_spec)
    a_vec <- null_constraint_grad_q(ctx, gd$grad)
    gth_full <- s_q - nu * a_vec
    gth_full[ctx$c0] <- 0

    dir_gn <- NULL
    dir_newton <- NULL
    dec <- Inf

    if (ctx$use_dense) {
      I_mat <- null_dense_info(ctx, ev$P)
      gth_red <- gth_full[ctx$umap]

      v_col <- I_mat[ctx$umap, ctx$c0]
      c00 <- I_mat[ctx$c0, ctx$c0]
      GN <- I_mat[ctx$umap, ctx$umap, drop = FALSE]
      tidx <- ctx$g_args_tidx
      if (length(tidx) > 0L) {
        gnz <- gd$grad[ctx$g_args_keep]
        GN[, tidx] <- GN[, tidx, drop = FALSE] + tcrossprod(v_col, gnz)
        GN[tidx, ] <- GN[tidx, , drop = FALSE] + tcrossprod(gnz, v_col)
        GN[tidx, tidx] <- GN[tidx, tidx, drop = FALSE] + c00 * tcrossprod(gnz)
      }

      gn_sol <- null_chol_solve_ladder(GN, gth_red)
      if (!is.null(gn_sol)) {
        d_gn_red <- gn_sol$x
        dec2 <- sum(gth_red * d_gn_red)
        dec <- sqrt(max(0, dec2))
        dir_gn <- numeric(q)
        dir_gn[ctx$umap] <- d_gn_red
      }

      kkt <- dec / sqrt(1 + abs(ev$ll))
      if (kkt < tol_kkt || iter >= maxit) {
        break
      }

      if (nu == 0 || !is_ph) {
        dir_newton <- dir_gn
      } else if (length(tidx) > 0L) {
        # Exact reduced Hessian: GN - nu * E_g
        gnz <- gd$grad[ctx$g_args_keep]
        kpz <- gd$kp[ctx$g_args_keep]
        E_sub <- (-tcrossprod(gnz, kpz) - tcrossprod(kpz, gnz) +
                    gd$k_sum * tcrossprod(gnz)) / gd$s_sum
        diag(E_sub) <- diag(E_sub) + kpz / gd$s_sum
        H_red <- GN
        H_red[tidx, tidx] <- H_red[tidx, tidx, drop = FALSE] - nu * E_sub
        d_nt_red <- null_chol_solve(H_red, gth_red)
        if (!is.null(d_nt_red)) {
          dir_newton <- numeric(q)
          dir_newton[ctx$umap] <- d_nt_red
        }
      } else {
        dir_newton <- dir_gn
      }
    } else {
      # Structured Woodbury + tangent-space PCG
      fac_res <- null_woodbury_ladder(ctx, ev$P)
      if (!is.null(fac_res)) {
        fac <- fac_res$fac
        y1 <- fac$solve(s_q)
        y2 <- fac$solve(a_vec)
        aty2 <- sum(a_vec * y2)
        if (is.finite(aty2) && aty2 > 0) {
          nup <- sum(a_vec * y1) / aty2
          d_gn <- y1 - nup * y2
          if (all(is.finite(d_gn))) {
            gv <- s_q - nup * a_vec
            dec2 <- sum(gv * d_gn)
            dec <- sqrt(max(0, dec2))
            dir_gn <- d_gn
          }
        }
      }

      kkt <- dec / sqrt(1 + abs(ev$ll))
      if (kkt < tol_kkt || iter >= maxit) {
        break
      }

      if (nu == 0 || !is_ph) {
        if (!is.null(fac_res) && fac_res$lambda == 0) {
          dir_newton <- dir_gn
        }
      } else if (!is.null(fac_res) && fac_res$lambda == 0 && !is.null(dir_gn)) {
        dir_newton <- null_newton_pcg(
          ctx = ctx,
          fac = fac_res$fac,
          y2 = y2,
          aty2 = aty2,
          a_vec = a_vec,
          nu = nu,
          gth_full = gth_full,
          gd = gd
        )
      }
    }

    iter <- iter + 1L
    accepted <- FALSE
    cands <- list()
    if (!is.null(dir_newton)) {
      cands[[length(cands) + 1L]] <- list(dir = dir_newton, is_newton = TRUE)
    }
    if (!is.null(dir_gn)) {
      cands[[length(cands) + 1L]] <- list(dir = dir_gn, is_newton = FALSE)
    }

    slack <- 1e-13 * (1 + abs(ev$ll))
    for (cand in cands) {
      dir <- cand$dir
      dir[ctx$c0] <- 0
      mx <- max(abs(dir))
      if (mx > max_step) {
        dir <- dir * (max_step / mx)
      }
      slope <- sum(gth_full * dir)
      if (!is.finite(slope) || slope <= 0) {
        next
      }
      dir_mat <- matrix(dir, nrow = p, ncol = ctx$Jm1)
      t_step <- 1
      for (bt in seq_len(max_backtrack)) {
        B_new <- B
        B_new[, ctx$cols] <- B_new[, ctx$cols] + t_step * dir_mat
        B_new[k, j] <- null_g_value(B_new[k, ], ctx$g_spec)
        if (all(is.finite(B_new)) && max(abs(B_new)) <= max_abs_B) {
          ev_new <- null_eval_ll(ctx, B_new)
          if (is.finite(ev_new$ll) && ev_new$ll >= ev$ll + c1 * t_step * slope - slack) {
            B <- B_new
            ev <- ev_new
            accepted <- TRUE
            if (cand$is_newton) {
              newton_steps <- newton_steps + 1L
            }
            break
          }
        }
        t_step <- t_step * 0.5
      }
      if (accepted) {
        break
      }
    }

    if (!accepted) {
      break
    }
  }

  # Restore caller's reference column to zero (as fit_null_symmetric does)
  B_out <- B - B[, j_ref]
  g_final <- null_g_value(B_out[k, ], ctx$g_spec)
  gap <- g_final - B_out[k, j]
  viol <- abs(gap)
  ev_final <- null_eval_ll(ctx, B_out)

  list(
    B = B_out,
    loglik = ev_final$ll,
    niter = iter,
    kkt = kkt,
    gap = gap,
    converged = is.finite(kkt) && is.finite(viol) && (kkt < tol_kkt) && (viol < 1e-10),
    newton_steps = newton_steps,
    local_solves = 1L
  )
}
