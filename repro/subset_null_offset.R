# Subset pseudo-Huber null fit with j_constr outside reference_set:
# omitting the fixed constrained offset X[, k_constr] * B[k_constr, j_constr]
# in micro_fisher() causes wrong-direction steps, oscillating test statistics,
# and non-convergence under tight tolerances.
suppressPackageStartupMessages({
  library(radEmu)
  library(fastEmu)
})

set.seed(11)
n <- 40; J <- 150; rs <- 1:30; j_test <- J
X <- cbind(`(Intercept)` = 1, x = seq(0, 1, length.out = n))
b0 <- rnorm(J, 0, 1.0); b0 <- b0 - mean(b0)
b1 <- rnorm(J, 0, 0.5); b1[j_test] <- 1.2; b1 <- b1 - mean(b1)
Y <- radEmu:::simulate_data(n = n, J = J, X = X, b0 = b0, b1 = b1, distn = "Poisson", mean_z = 2.5)
colnames(Y) <- paste0("taxon", seq_len(J))

mk_cf <- function(idx) {
  f <- function(x) pseudohuber_median(x[idx])
  attr(f, "reference_set") <- idx
  f
}
mk_cg <- function(idx) {
  function(x) {
    g <- rep(0, length(x))
    g[idx] <- dpseudohuber_median_dx(x[idx])
    g
  }
}
cf <- list(mk_cf(rs), mk_cf(rs))
cg <- list(mk_cg(rs), mk_cg(rs))

# 1. Unconstrained fit, then direct fit_null_symmetric_subset() with tight tolerances
fit_unc <- suppressMessages(emuFit(
  Y = Y, X = X, run_score_tests = FALSE, compute_cis = FALSE,
  constraint_fn = cf, constraint_grad_fn = cg, tolerance = 1e-6,
  match_row_names = FALSE
))

cf_typed <- cf
for (k in 1:2) attr(cf_typed[[k]], "constraint_type") <- "symmetric_subset:pseudohuber"
j_ref <- rs[which.max(colSums(fit_unc$Y_augmented > 1)[rs])]

res <- radEmu:::fit_null_symmetric_subset(
  B = fit_unc$B, Y = fit_unc$Y_augmented, X = X,
  X_cup = radEmu:::X_cup_from_X(X, J),
  k_constr = 2, j_constr = j_test, j_ref = j_ref,
  constraint_fn = cf_typed, constraint_grad_fn = cg,
  tol_lik = 1e-8, tol_test_stat = 1e-4, maxit = 120
)
df <- res$it_df
dlik <- diff(df$lik)
dstat <- abs(diff(df$test_stat) / df$test_stat[-nrow(df)])
cat(sprintf(
  "fit_null_symmetric_subset (n=%d, J=%d, |ref|=%d, j=%d, tol_lik=1e-8, tol_test_stat=1e-4, maxit=120):\n",
  n, J, length(rs), j_test
))
cat(sprintf(
  "  converged = %s, iterations = %d, final ll = %.6f, downhill iters = %d (min dlik = %+.3e)\n",
  res$converged, nrow(df), tail(df$lik, 1), sum(dlik < -1e-10), min(dlik)
))
cat("  last 6 iterations:\n")
for (i in (nrow(df) - 5):nrow(df)) {
  cat(sprintf(
    "    iter %3d: ll = %13.6f (dlik = %+10.3e)  stat = %.6f (rel dstat = %.3e)\n",
    i, df$lik[i], df$lik[i] - df$lik[i - 1],
    df$test_stat[i], abs((df$test_stat[i] - df$test_stat[i - 1]) / df$test_stat[i - 1])
  ))
}

# 2. Isolate the single micro_fisher() call on j_constr at this null point (working reference j_ref)
B_work <- res$B
for (k in 1:2) B_work[k, ] <- B_work[k, ] - B_work[k, j_ref]
z_work <- radEmu:::update_z(fit_unc$Y_augmented, X, B_work)
col_ll <- function(b0_val) {
  lm <- as.numeric(X[, -2, drop = FALSE] %*% b0_val + X[, 2] * B_work[2, j_test] + z_work)
  sum(fit_unc$Y_augmented[, j_test] * lm - exp(lm))
}
ll0 <- col_ll(B_work[1, j_test])
step_no_off <- radEmu:::micro_fisher(
  X[, -2, drop = FALSE],
  Yj = fit_unc$Y_augmented[, j_test, drop = FALSE],
  Bj = B_work[-2, j_test, drop = FALSE],
  z = z_work, stepsize = 1, c1 = 1e-2
)
step_with_off <- radEmu:::micro_fisher(
  X[, -2, drop = FALSE],
  Yj = fit_unc$Y_augmented[, j_test, drop = FALSE],
  Bj = B_work[-2, j_test, drop = FALSE],
  z = z_work + as.numeric(X[, 2] * B_work[2, j_test]), stepsize = 1, c1 = 1e-2
)
cat(sprintf(
  "Single micro_fisher() step on j_constr (working B[2, j_constr] = %+.4f):\n",
  B_work[2, j_test]
))
cat(sprintf(
  "  without offset: step = %+.4e, true column dlik = %+.4e (wrong direction)\n",
  step_no_off, col_ll(B_work[1, j_test] + step_no_off) - ll0
))
cat(sprintf(
  "  with offset:    step = %+.4e, true column dlik = %+.4e\n",
  step_with_off, col_ll(B_work[1, j_test] + step_with_off) - ll0
))

# 3. Same problem through emuFit(..., null_fit_alg = "constraint_sandwich")
msgs <- character(0)
fit_emu <- withCallingHandlers(
  emuFit(
    Y = Y, X = X, fitted_model = fit_unc, refit = FALSE, compute_cis = FALSE,
    test_kj = data.frame(k = 2, j = j_test),
    constraint_fn = cf, constraint_grad_fn = cg,
    null_fit_alg = "constraint_sandwich",
    tol_lik = 1e-8, tol_test_stat = 1e-4, maxit_null = 120,
    null_diagnostic_plots = TRUE, match_row_names = FALSE
  ),
  message = function(m) {
    msgs <<- c(msgs, trimws(conditionMessage(m)))
    invokeRestart("muffleMessage")
  }
)
diag_df <- fit_emu$null_diagnostic_plots[[1]]$diagnostics_df
cat(sprintf(
  "emuFit (tight tol): fallback_to_AL = %s, iters = %d, stat = %.6f, p = %.6f\n",
  "constraint_diff" %in% names(diag_df), nrow(diag_df),
  fit_emu$coef$score_stat[j_test], fit_emu$coef$pval[j_test]
))

# 4. Impact on fastEmuFit() (which defaults to fit_null_symmetric_subset() on continuous designs)
fit_fast <- suppressMessages(fastEmuFit(
  Y = Y, X = X, test_kj = data.frame(k = 2, j = 40),
  reference_set = "data_driven",
  tol_lik = 1e-8, tol_test_stat = 1e-5, maxit_null = 200,
  match_row_names = FALSE
))
fit_fast_al <- suppressMessages(fastEmuFit(
  Y = Y, X = X, test_kj = data.frame(k = 2, j = 40),
  reference_set = "data_driven",
  null_fit_alg = "augmented_lagrangian",
  match_row_names = FALSE
))
r_fast <- fit_fast$coef[fit_fast$coef$category == "taxon40", ]
r_al <- fit_fast_al$coef[fit_fast_al$coef$category == "taxon40", ]
cat(sprintf(
  "fastEmuFit (j=40): default CS stat = %.6f (p = %.6f); AL stat = %.6f (p = %.6f)\n",
  r_fast$score_stat, r_fast$pval, r_al$score_stat, r_al$pval
))
