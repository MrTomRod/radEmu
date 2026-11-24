suppressMessages(library(radEmu))
set.seed(5)
n <- 30; J <- 25
grp <- rep(0:1, each = n / 2)
Y <- matrix(rpois(n * J, rep(exp(rnorm(J, 3, 0.7)), each = n) * ifelse(grp == 1, exp(rep(c(rnorm(12, 0, 0.1), rnorm(J - 12, 1, 0.8)), each = n)), 1)), n, J,
            dimnames = list(NULL, paste0("cat", 1:J)))
X <- cbind(`(Intercept)` = 1, g = grp)
rs <- 1:12   # reference set for the covariate row: pseudo-Huber median over taxa 1..12 only
mk <- function() {
  f <- function(x) pseudohuber_median(x[rs]); attr(f, "reference_set") <- rs; f }
g <- function(x) { gr <- rep(0, length(x)); gr[rs] <- dpseudohuber_median_dx(x[rs]); gr }
cf <- list(function(x) pseudohuber_median(x), mk())
cg <- list(function(x) dpseudohuber_median_dx(x), g)
j <- 20   # taxon outside the reference set (taxa 13..25 carry a true group effect, taxa 1..12 almost none)
run <- function(alg) suppressMessages(emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = j),
        constraint_fn = cf, constraint_grad_fn = cg, null_fit_alg = alg, return_nullB = TRUE, verbose = FALSE))
Yaug <- NULL
for (alg in c(NA, "augmented_lagrangian")) {
  fit <- if (is.na(alg)) run(NULL) else run(alg)
  Yaug <- fit$Y_augmented; B <- fit$null_B[[1]]
  lm <- X %*% B; z <- log(rowSums(Yaug)) - log(rowSums(exp(lm))); ll <- sum(Yaug * (lm + z) - exp(lm + z))
  cat(sprintf("null_fit_alg = %-22s stat = %.4f  p = %.4f  null loglik = %.4f\n   null constraint residual with the user's subset constraint: B_kj - g_subset(B_k) = %+.3e ; with the all-taxa median: %+.3e\n",
      if (is.na(alg)) "default (NULL)" else alg, fit$coef$score_stat[fit$coef$covariate == "g" & fit$coef$category == paste0("cat", j)],
      fit$coef$pval[fit$coef$covariate == "g" & fit$coef$category == paste0("cat", j)], ll,
      B[2, j] - pseudohuber_median(B[2, rs], tolerance = 1e-14), B[2, j] - pseudohuber_median(B[2, ], tolerance = 1e-14)))
}
