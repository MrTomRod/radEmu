#' Exact pseudo-Huber median and derivatives (internal helpers)
#'
#' Safeguarded Newton evaluation of the pseudo-Huber smoothed median
#' \eqn{c = \arg\min_c \sum_l d^2 (\sqrt{1 + ((x_l - c)/d)^2} - 1)},
#' its gradient and Hessian parts, and deterministic classification of a
#' constraint function for \code{fit_null_multistart()}.
#'
#' Unlike \code{pseudohuber_median()}, which uses linearly convergent IRLS
#' stopping at \code{1e-8}, safeguarded Newton on the monotone stationarity
#' equation \eqn{F(c) = \sum_l (x_l - c) / \sqrt{1 + ((x_l - c)/d)^2} = 0}
#' reaches machine precision in a handful of iterations. This allows the
#' reduced Newton null solver to certify KKT residuals to \code{1e-9} and
#' constraint residuals below \code{1e-10}.
#'
#' @param x Numeric vector.
#' @param d Positive numeric smoothing parameter.
#' @return \code{pseudohuber_median_exact()} returns a numeric scalar;
#'   \code{pseudohuber_derivatives_exact()} returns a list with elements
#'   \code{value}, \code{grad}, \code{kp}, \code{s_sum}, and \code{k_sum}
#'   representing \eqn{\nabla g = \kappa / S} and
#'   \eqn{\nabla^2 g = (\mathrm{diag}(\kappa') - G \kappa'^\top - \kappa' G^\top + K' G G^\top) / S}.
#' @keywords internal
#' @noRd
pseudohuber_median_exact <- function(x, d = 0.1) {
  if (length(x) == 0L || any(!is.finite(x))) {
    return(NA_real_)
  }
  lo <- min(x)
  hi <- max(x)
  if (lo == hi) {
    return(lo)
  }
  c_val <- stats::median(x)
  tol_scale <- 2 * .Machine$double.eps
  for (iter in seq_len(200L)) {
    r <- (x - c_val) / d
    w <- 1 / sqrt(1 + r * r)
    f <- sum(d * r * w)
    if (f == 0) {
      return(c_val)
    }
    if (f > 0) {
      lo <- c_val
    } else {
      hi <- c_val
    }
    fp <- -sum(w * w * w)
    cn <- c_val - f / fp
    if (!is.finite(cn) || cn <= lo || cn >= hi) {
      cn <- 0.5 * (lo + hi)
    }
    if (abs(cn - c_val) <= tol_scale * max(abs(c_val), 1)) {
      return(cn)
    }
    c_val <- cn
  }
  c_val
}

#' @rdname pseudohuber_median_exact
#' @keywords internal
#' @noRd
pseudohuber_derivatives_exact <- function(x, d = 0.1) {
  c_val <- pseudohuber_median_exact(x, d = d)
  r <- x - c_val
  inv_d <- 1 / d
  w <- 1 / sqrt(1 + (r * inv_d)^2)
  w2 <- w * w
  w3 <- w2 * w
  w5 <- w3 * w2
  s_sum <- sum(w3)
  grad <- w3 / s_sum
  kp <- -3 * r * (inv_d * inv_d) * w5
  k_sum <- sum(kp)
  list(
    value = c_val,
    grad = grad,
    kp = kp,
    s_sum = s_sum,
    k_sum = k_sum
  )
}

#' Classify a row constraint function for the reduced null solver
#'
#' Deterministically identifies whether \code{constraint_fn} /
#' \code{constraint_grad_fn} implement a pseudo-Huber median (all categories or
#' a subset, any \code{d > 0}), an arithmetic mean (all categories or a subset),
#' or a single-category reference constraint. Does not touch the RNG state.
#'
#' @param constraint_fn Single constraint function \eqn{g(x)}.
#' @param constraint_grad_fn Single constraint gradient function \eqn{\nabla g(x)}.
#' @param J Number of categories.
#' @param j_constr Tested category index (1-based).
#' @return A list with elements \code{kind} (\code{"pseudohuber"},
#'   \code{"mean"}, \code{"reference"}, or \code{"other"}), \code{args}
#'   (1-based argument indices with \code{j_constr} removed), \code{cols}
#'   (full support indices), and \code{d} (smoothing parameter for
#'   \code{"pseudohuber"}, else \code{NA_real_}).
#' @keywords internal
#' @noRd
classify_null_constraint <- function(constraint_fn, constraint_grad_fn, J, j_constr) {
  other <- list(kind = "other", args = integer(0), cols = integer(0), d = NA_real_)
  if (!is.function(constraint_fn) || !is.function(constraint_grad_fn) || J < 2L) {
    return(other)
  }

  # Deterministic test vectors in [-2.1, 2.1] so pseudo-Huber weights never underflow
  idx <- seq_len(J)
  v1 <- seq(-2, 2, length.out = J) + 0.1 * sin(idx * 1.7)
  v2 <- 1.7 * rev(v1) + 0.6

  eval_pair <- tryCatch(
    list(
      f1 = as.numeric(constraint_fn(v1)),
      g1 = as.numeric(constraint_grad_fn(v1)),
      f2 = as.numeric(constraint_fn(v2)),
      g2 = as.numeric(constraint_grad_fn(v2))
    ),
    error = function(e) NULL
  )
  if (is.null(eval_pair)) {
    return(other)
  }
  f1 <- eval_pair$f1
  g1 <- eval_pair$g1
  f2 <- eval_pair$f2
  g2 <- eval_pair$g2
  if (length(f1) != 1L || length(f2) != 1L ||
      length(g1) != J || length(g2) != J ||
      !is.finite(f1) || !is.finite(f2) ||
      any(!is.finite(g1)) || any(!is.finite(g2))) {
    return(other)
  }

  ref_set <- attr(constraint_fn, "reference_set")
  if (!is.null(ref_set)) {
    S <- sort(unique(as.integer(ref_set)))
    if (length(S) == 0L || any(is.na(S)) || any(S < 1L | S > J)) {
      return(other)
    }
    if (length(S) < J &&
        (max(abs(g1[-S])) > 1e-12 || max(abs(g2[-S])) > 1e-12)) {
      return(other)
    }
  } else {
    S <- which(g1 != 0 | g2 != 0)
    if (length(S) == 0L) {
      return(other)
    }
  }

  # Single-category reference constraint: g(x) = x[c]
  if (length(S) == 1L) {
    if (abs(f1 - v1[S]) < 1e-10 && abs(f2 - v2[S]) < 1e-10 &&
        abs(g1[S] - 1) < 1e-10 && abs(g2[S] - 1) < 1e-10) {
      return(list(
        kind = "reference",
        args = setdiff(S, j_constr),
        cols = S,
        d = NA_real_
      ))
    }
    return(other)
  }

  # Mean constraint on S
  m1 <- mean(v1[S])
  m2 <- mean(v2[S])
  u_wt <- 1 / length(S)
  if (abs(f1 - m1) < 1e-8 * (1 + abs(m1)) &&
      abs(f2 - m2) < 1e-8 * (1 + abs(m2)) &&
      max(abs(g1[S] - u_wt)) < 1e-8 &&
      max(abs(g2[S] - u_wt)) < 1e-8) {
    return(list(
      kind = "mean",
      args = setdiff(S, j_constr),
      cols = S,
      d = NA_real_
    ))
  }

  # Pseudo-Huber median on S: gather candidate d from closure environment, default 0.1,
  # and direct algebraic recovery from the gradient weights w_i^3 / sum(w^3).
  cands <- numeric(0)
  env <- environment(constraint_fn)
  if (is.environment(env)) {
    for (nm in c("d", "constraint_param")) {
      if (exists(nm, envir = env, inherits = FALSE)) {
        val <- get(nm, envir = env, inherits = FALSE)
        if (is.numeric(val) && length(val) == 1L && is.finite(val) && val > 0) {
          cands <- c(cands, val)
        }
      }
    }
  }
  cands <- c(cands, 0.1)

  gs <- g1[S]
  if (all(gs > 0)) {
    r2 <- (v1[S] - f1)^2
    i_min <- which.min(r2)
    i_max <- which.max(r2)
    if (gs[i_min] > gs[i_max] && r2[i_max] > r2[i_min]) {
      ratio <- (gs[i_min] / gs[i_max])^(2 / 3)
      d2 <- (r2[i_max] - ratio * r2[i_min]) / (ratio - 1)
      if (is.finite(d2) && d2 > 0) {
        cands <- c(cands, sqrt(d2))
      }
    }
  }
  cands <- unique(cands[is.finite(cands) & cands > 0])

  for (d_cand in cands) {
    pf1 <- pseudohuber_median(v1[S], d = d_cand)
    pf2 <- pseudohuber_median(v2[S], d = d_cand)
    pg1 <- dpseudohuber_median_dx(v1[S], d = d_cand)
    pg2 <- dpseudohuber_median_dx(v2[S], d = d_cand)
    if (abs(f1 - pf1) < 1e-6 * (1 + abs(pf1)) &&
        abs(f2 - pf2) < 1e-6 * (1 + abs(pf2)) &&
        max(abs(g1[S] - pg1)) < 1e-5 &&
        max(abs(g2[S] - pg2)) < 1e-5) {
      return(list(
        kind = "pseudohuber",
        args = setdiff(S, j_constr),
        cols = S,
        d = d_cand
      ))
    }
  }

  other
}
