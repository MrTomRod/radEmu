test_that("pseudohuber_median_exact and derivatives match stationarity and finite differences", {
  x <- c(-3.2, -1.1, -0.05, 0.02, 0.08, 0.4, 2.7, 5.1)
  for (d in c(0.01, 0.1, 1.0, 5.0)) {
    gd <- radEmu:::pseudohuber_derivatives_exact(x, d = d)
    r <- (x - gd$value) / d
    stat_resid <- sum(d * r / sqrt(1 + r * r))
    expect_lt(abs(stat_resid), 1e-13)

    # Elimination identity: g(x with x_j := g(x_{-j})) == g(x_{-j})
    for (j in seq_along(x)) {
      xj <- radEmu:::pseudohuber_median_exact(x[-j], d = d)
      x_mod <- x
      x_mod[j] <- xj
      expect_lt(abs(radEmu:::pseudohuber_median_exact(x_mod, d = d) - xj), 1e-13)
    }

    # Gradient and Hessian vs central finite differences
    eps <- 1e-6
    m <- length(x)
    fd_grad <- numeric(m)
    fd_hess <- matrix(0, m, m)
    for (i in seq_len(m)) {
      xp <- x
      xm <- x
      xp[i] <- xp[i] + eps
      xm[i] <- xm[i] - eps
      fd_grad[i] <- (radEmu:::pseudohuber_median_exact(xp, d = d) -
                       radEmu:::pseudohuber_median_exact(xm, d = d)) / (2 * eps)
      gp <- radEmu:::pseudohuber_derivatives_exact(xp, d = d)$grad
      gm <- radEmu:::pseudohuber_derivatives_exact(xm, d = d)$grad
      fd_hess[, i] <- (gp - gm) / (2 * eps)
    }
    hess <- (-tcrossprod(gd$grad, gd$kp) - tcrossprod(gd$kp, gd$grad) +
               gd$k_sum * tcrossprod(gd$grad)) / gd$s_sum
    diag(hess) <- diag(hess) + gd$kp / gd$s_sum

    expect_lt(max(abs(gd$grad - fd_grad)), 1e-7)
    expect_lt(max(abs(hess - fd_hess)), 1e-6)
  }
})

test_that("classify_null_constraint identifies constraints deterministically without touching RNG", {
  set.seed(123)
  seed_before <- .Random.seed
  J <- 12L

  # Default and custom-d pseudo-Huber
  c_ph <- radEmu:::classify_null_constraint(
    function(x) pseudohuber_median(x, d = 0.1),
    function(x) dpseudohuber_median_dx(x, d = 0.1),
    J, 3L
  )
  expect_equal(c_ph$kind, "pseudohuber")
  expect_equal(c_ph$d, 0.1, tolerance = 1e-8)
  expect_equal(c_ph$args, setdiff(seq_len(J), 3L))

  c_ph2 <- radEmu:::classify_null_constraint(
    function(x) pseudohuber_median(x, d = 0.25),
    function(x) dpseudohuber_median_dx(x, d = 0.25),
    J, 3L
  )
  expect_equal(c_ph2$kind, "pseudohuber")
  expect_equal(c_ph2$d, 0.25, tolerance = 1e-6)

  # Subset pseudo-Huber (with and without reference_set attribute)
  rs <- 1:6
  f_sub <- function(x) pseudohuber_median(x[rs], d = 0.1)
  g_sub <- function(x) {
    gr <- rep(0, length(x))
    gr[rs] <- dpseudohuber_median_dx(x[rs], d = 0.1)
    gr
  }
  c_sub <- radEmu:::classify_null_constraint(f_sub, g_sub, J, 3L)
  expect_equal(c_sub$kind, "pseudohuber")
  expect_equal(c_sub$args, setdiff(rs, 3L))

  # Mean and subset mean
  c_mean <- radEmu:::classify_null_constraint(
    mean, function(x) rep(1 / length(x), length(x)), J, 4L
  )
  expect_equal(c_mean$kind, "mean")
  expect_equal(c_mean$args, setdiff(seq_len(J), 4L))

  # Reference constraint
  c_ref <- radEmu:::classify_null_constraint(
    function(x) x[5L],
    function(x) { g <- rep(0, length(x)); g[5L] <- 1; g },
    J, 2L
  )
  expect_equal(c_ref$kind, "reference")
  expect_equal(c_ref$args, 5L)

  # Unknown custom constraint falls back to "other"
  c_other <- radEmu:::classify_null_constraint(
    function(x) log(sum(exp(x))),
    function(x) exp(x) / sum(exp(x)),
    J, 2L
  )
  expect_equal(c_other$kind, "other")

  # RNG state untouched
  expect_identical(.Random.seed, seed_before)
})

test_that("dense and structured local solvers agree on discrete and continuous designs", {
  set.seed(42)
  n <- 16
  J <- 12
  X <- cbind(Intercept = 1, x1 = rnorm(n), x2 = rep(c(0, 1), each = n / 2))
  eta <- X %*% rbind(rnorm(J, 0, 0.4), rnorm(J, 0, 0.8), rnorm(J, 0, 0.6))
  Y <- matrix(rpois(n * J, exp(eta + 3)), nrow = n, ncol = J)
  Y <- Y + 0.5
  B_start <- matrix(0, nrow = ncol(X), ncol = J)
  j_ref <- radEmu:::get_j_ref(Y)

  fn <- function(x) pseudohuber_median(x, d = 0.1)
  gr <- function(x) dpseudohuber_median_dx(x, d = 0.1)

  for (j in c(2L, j_ref)) {
    f_dense <- radEmu:::fit_null_multistart(
      B = B_start, Y = Y, X = X, k_constr = 2L, j_constr = j, j_ref = j_ref,
      constraint_fn = fn, constraint_grad_fn = gr, solver = "dense"
    )
    f_struct <- radEmu:::fit_null_multistart(
      B = B_start, Y = Y, X = X, k_constr = 2L, j_constr = j, j_ref = j_ref,
      constraint_fn = fn, constraint_grad_fn = gr, solver = "structured"
    )
    expect_true(f_dense$converged)
    expect_true(f_struct$converged)
    expect_lt(abs(f_dense$gap), 1e-10)
    expect_lt(abs(f_struct$gap), 1e-10)
    expect_lt(abs(f_dense$loglik - f_struct$loglik), 1e-9 * (1 + abs(f_dense$loglik)))
    expect_lt(max(abs(f_dense$B - f_struct$B)), 1e-5)
  }
})

test_that("multistart satisfies constraint and matches or beats existing algorithms on small problems", {
  poisson_ll <- function(Y, X, B) {
    lm <- X %*% B
    z <- log(rowSums(Y)) - log(rowSums(exp(lm)))
    sum(Y * (lm + z) - exp(lm + z))
  }

  # Small J < 20 (continuous design): must NOT be overwritten by J < 20 AL default
  set.seed(2103)
  n <- 10
  J <- 6
  X <- cbind(Intercept = 1, x1 = rnorm(n))
  Y <- matrix(rpois(n * J, exp(3 + X %*% rbind(rnorm(J, 0, 0.5), c(-1.5, -0.5, 0, 0.2, 0.8, 1.5)))),
              nrow = n, ncol = J)
  colnames(Y) <- paste0("cat", seq_len(J))
  rownames(Y) <- rownames(X) <- paste0("s", seq_len(n))
  tk <- data.frame(k = 2, j = c(1, 3, 6))

  fit_ms <- emuFit(Y = Y, X = X, test_kj = tk, null_fit_alg = "multistart",
                   return_nullB = TRUE, tolerance = 1e-8, verbose = FALSE)
  fit_al <- emuFit(Y = Y, X = X, test_kj = tk, null_fit_alg = "augmented_lagrangian",
                   return_nullB = TRUE, tolerance = 1e-8, B_null_tol = 1e-6,
                   constraint_tol = 1e-7, verbose = FALSE)

  expect_true(all(fit_ms$score_test_hyperparams$converged == "converged"))
  for (i in seq_len(nrow(tk))) {
    B_ms <- fit_ms$null_B[[i]]
    B_al <- fit_al$null_B[[i]]
    j <- tk$j[i]
    # Constraint satisfied to tight tolerance
    resid_exact <- abs(B_ms[2, j] - radEmu:::pseudohuber_median_exact(B_ms[2, ], d = 0.1))
    expect_lt(resid_exact, 1e-10)
    # Null log-likelihood >= augmented Lagrangian
    ll_ms <- poisson_ll(fit_ms$Y_augmented, X, B_ms)
    ll_al <- poisson_ll(fit_al$Y_augmented, X, B_al)
    expect_gte(ll_ms, ll_al - 1e-6 * (1 + abs(ll_al)))
  }

  # Moderate J = 20 (discrete two-group design): compare against discrete, constraint_sandwich, AL
  set.seed(2107)
  n2 <- 16
  J2 <- 20
  X2 <- cbind(Intercept = 1, grp = rep(c(0, 1), each = n2 / 2))
  b1 <- c(rnorm(10, -0.5, 0.8), rnorm(10, 0.8, 0.8))
  Y2 <- matrix(rpois(n2 * J2, exp(2.5 + X2 %*% rbind(rnorm(J2, 0, 0.6), b1))),
               nrow = n2, ncol = J2)
  colnames(Y2) <- paste0("cat", seq_len(J2))
  rownames(Y2) <- rownames(X2) <- paste0("s", seq_len(n2))
  tk2 <- data.frame(k = 2, j = c(2, 10, 18))

  fit_base <- emuFit(Y = Y2, X = X2, run_score_tests = FALSE, tolerance = 1e-8, verbose = FALSE)
  fit_ms2 <- emuFit(Y = Y2, X = X2, fitted_model = fit_base, refit = FALSE,
                    test_kj = tk2, null_fit_alg = "multistart",
                    return_nullB = TRUE, verbose = FALSE)
  fit_disc <- emuFit(Y = Y2, X = X2, fitted_model = fit_base, refit = FALSE,
                     test_kj = tk2, null_fit_alg = "discrete",
                     return_nullB = TRUE, tol_discrete = 1e-5, verbose = FALSE)
  fit_cs <- emuFit(Y = Y2, X = X2, fitted_model = fit_base, refit = FALSE,
                   test_kj = tk2, null_fit_alg = "constraint_sandwich",
                   return_nullB = TRUE, tol_lik = 1e-8, tol_test_stat = 1e-4, verbose = FALSE)

  Ya <- fit_base$Y_augmented
  expect_true(all(fit_ms2$score_test_hyperparams$converged == "converged"))
  for (i in seq_len(nrow(tk2))) {
    B_ms <- fit_ms2$null_B[[i]]
    j <- tk2$j[i]
    expect_lt(abs(B_ms[2, j] - radEmu:::pseudohuber_median_exact(B_ms[2, ], d = 0.1)), 1e-10)
    ll_ms <- poisson_ll(Ya, X2, B_ms)
    ll_disc <- poisson_ll(Ya, X2, fit_disc$null_B[[i]])
    ll_cs <- poisson_ll(Ya, X2, fit_cs$null_B[[i]])
    expect_gte(ll_ms, ll_disc - 1e-6 * (1 + abs(ll_disc)))
    expect_gte(ll_ms, ll_cs - 1e-6 * (1 + abs(ll_cs)))
  }
})

test_that("multistart handles mean, subset pseudo-Huber, reference, and custom constraints", {
  set.seed(99)
  n <- 12
  J <- 10
  X <- cbind(Intercept = 1, x1 = rnorm(n))
  Y <- matrix(rpois(n * J, exp(3 + X %*% rbind(rnorm(J, 0, 0.3), rnorm(J, 0, 0.6)))),
              nrow = n, ncol = J)
  colnames(Y) <- paste0("cat", seq_len(J))
  rownames(Y) <- rownames(X) <- paste0("s", seq_len(n))

  # Mean constraint
  fit_mean <- emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 3),
                     constraint_fn = mean,
                     constraint_grad_fn = function(x) rep(1 / length(x), length(x)),
                     null_fit_alg = "multistart", return_nullB = TRUE, verbose = FALSE)
  expect_equal(fit_mean$score_test_hyperparams$converged, "converged")
  expect_lt(abs(fit_mean$null_B[[1]][2, 3] - mean(fit_mean$null_B[[1]][2, ])), 1e-10)

  # Subset pseudo-Huber constraint
  rs <- 1:6
  f_sub <- function(x) pseudohuber_median(x[rs], d = 0.1)
  attr(f_sub, "reference_set") <- rs
  g_sub <- function(x) {
    gr <- rep(0, length(x))
    gr[rs] <- dpseudohuber_median_dx(x[rs], d = 0.1)
    gr
  }
  fit_sub <- emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 8),
                    constraint_fn = f_sub, constraint_grad_fn = g_sub,
                    null_fit_alg = "multistart", return_nullB = TRUE, verbose = FALSE)
  expect_equal(fit_sub$score_test_hyperparams$converged, "converged")
  expect_lt(abs(fit_sub$null_B[[1]][2, 8] -
                  radEmu:::pseudohuber_median_exact(fit_sub$null_B[[1]][2, rs], d = 0.1)), 1e-10)

  # Single-category reference constraint
  ref_c <- make_reference_constraints(p = 2, j = 2)
  fit_ref <- emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 4),
                    constraint_fn = ref_c$constraint_list,
                    constraint_grad_fn = ref_c$constraint_grad_list,
                    null_fit_alg = "multistart", return_nullB = TRUE, verbose = FALSE)
  expect_equal(fit_ref$score_test_hyperparams$converged, "converged")
  expect_lt(abs(fit_ref$null_B[[1]][2, 4] - fit_ref$null_B[[1]][2, 2]), 1e-10)

  # Custom d pseudo-Huber constraint
  fit_d <- emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 5),
                  constraint_param = 0.25,
                  null_fit_alg = "multistart", return_nullB = TRUE, verbose = FALSE)
  expect_equal(fit_d$score_test_hyperparams$converged, "converged")
  expect_lt(abs(fit_d$null_B[[1]][2, 5] -
                  radEmu:::pseudohuber_median_exact(fit_d$null_B[[1]][2, ], d = 0.25)), 1e-10)

  # Non-standard custom constraint falls back cleanly to augmented Lagrangian
  f_lse <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
  g_lse <- function(x) { e <- exp(x - max(x)); e / sum(e) }
  fit_other <- emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 5),
                      constraint_fn = f_lse, constraint_grad_fn = g_lse,
                      null_fit_alg = "multistart", return_nullB = TRUE, verbose = FALSE)
  expect_equal(fit_other$score_test_hyperparams$converged, "converged")
  expect_lt(abs(fit_other$null_B[[1]][2, 5] - f_lse(fit_other$null_B[[1]][2, ])), 1e-4)
})

test_that("multistart matches or beats all existing null_fit_alg options on wirbel_otu_small", {
  data("wirbel_otu_small", package = "radEmu", envir = environment())
  data("wirbel_sample_small", package = "radEmu", envir = environment())
  Y <- as.matrix(wirbel_otu_small)
  storage.mode(Y) <- "double"
  X <- cbind(Intercept = 1, CRC = as.numeric(wirbel_sample_small$Group == "CRC"))
  rownames(Y) <- rownames(X) <- sprintf("s%03d", seq_len(nrow(Y)))

  poisson_ll <- function(Ya, X, B) {
    lm <- X %*% B
    z <- log(rowSums(Ya)) - log(rowSums(exp(lm)))
    sum(Ya * (lm + z) - exp(lm + z))
  }

  # Taxa 17 and 24 have multiple null local optima under the pseudo-Huber constraint;
  # single-start solvers land well below the best optimum on one or both.
  tk <- data.frame(k = 2, j = c(17L, 24L))
  fit_base <- emuFit(
    Y = Y, X = X, run_score_tests = FALSE,
    tolerance = 1e-10, maxit = 2000, verbose = FALSE
  )
  Ya <- fit_base$Y_augmented

  fit_ms <- emuFit(
    Y = Y, X = X, fitted_model = fit_base, refit = FALSE,
    test_kj = tk, null_fit_alg = "multistart",
    return_nullB = TRUE, verbose = FALSE
  )
  expect_true(all(fit_ms$score_test_hyperparams$converged == "converged"))

  ll_ms <- vapply(seq_len(nrow(tk)), function(i) {
    poisson_ll(Ya, X, fit_ms$null_B[[i]])
  }, numeric(1))

  # Tight constraint satisfaction and finite positive score statistics
  for (i in seq_len(nrow(tk))) {
    B0 <- fit_ms$null_B[[i]]
    expect_lt(abs(B0[2, tk$j[i]] - radEmu:::pseudohuber_median_exact(B0[2, ], d = 0.1)), 1e-10)
    row_i <- fit_ms$coef[fit_ms$coef$covariate == "CRC" & fit_ms$coef$category_num == tk$j[i], ]
    expect_true(is.finite(row_i$score_stat) && row_i$score_stat > 0)
  }

  # Null log-likelihood >= every existing null_fit_alg, and strictly higher on the
  # multi-modal cases where single-start solvers get trapped in a suboptimal sheet.
  for (alg in c("discrete", "constraint_sandwich", "augmented_lagrangian")) {
    fit_other <- emuFit(
      Y = Y, X = X, fitted_model = fit_base, refit = FALSE,
      test_kj = tk, null_fit_alg = alg,
      return_nullB = TRUE, verbose = FALSE
    )
    ll_other <- vapply(seq_len(nrow(tk)), function(i) {
      poisson_ll(Ya, X, fit_other$null_B[[i]])
    }, numeric(1))
    expect_gte(ll_ms[1], ll_other[1] - 1e-6 * (1 + abs(ll_other[1])))
    expect_gte(ll_ms[2], ll_other[2] - 1e-6 * (1 + abs(ll_other[2])))
    # On taxon 17 all three single-start algorithms land >100 LL units below multistart
    expect_gt(ll_ms[1] - ll_other[1], 100)
    if (alg != "discrete") {
      # On taxon 24 constraint_sandwich and augmented_lagrangian land >40 LL units below
      expect_gt(ll_ms[2] - ll_other[2], 40)
    }
  }
})
