# wirbel_otu_small (320 x 47), CRC vs control; one score test: k = 2 (CRC), j = 24.
# Tightening tolerances lands ~45 log-likelihood units lower than default settings and flips p across 0.05.
suppressMessages(library(radEmu))
data("wirbel_otu_small", package = "radEmu")
data("wirbel_sample_small", package = "radEmu")
Y <- as.matrix(wirbel_otu_small)
storage.mode(Y) <- "double"
X <- cbind(`(Intercept)` = 1, CRC = as.numeric(wirbel_sample_small$Group == "CRC"))
rownames(Y) <- rownames(X) <- sprintf("s%03d", seq_len(nrow(Y)))
tk <- data.frame(k = 2, j = 24)
tight <- list(tolerance = 1e-10, maxit = 2000, tol_lik = 1e-10, tol_test_stat = 1e-6, maxit_null = 1000,
              tol_discrete = 1e-7, B_null_tol = 1e-8, constraint_tol = 1e-9, inner_tol = 1e-6, inner_maxit = 50)
nll <- function(fit) {
  B <- fit$null_B[[1]]; Ya <- fit$Y_augmented; lm <- X %*% B
  z <- log(rowSums(Ya)) - log(rowSums(exp(lm)))
  sum(Ya * (lm + z) - exp(lm + z))
}
for (nm in c("default", "tight")) {
  msgs <- character(0)
  t0 <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    do.call(emuFit, c(list(Y = Y, X = X, test_kj = tk, return_nullB = TRUE, verbose = FALSE), if (nm == "tight") tight)),
    message = function(m) { msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })
  r <- fit$coef[fit$coef$covariate == "CRC" & fit$coef$category_num == 24, ]
  cat(sprintf("%-8s statistic %.4f  p %.4f  null loglik %.3f   (%.0f s)\n", nm, r$score_stat, r$pval, nll(fit), proc.time()[["elapsed"]] - t0))
  for (m in unique(msgs[!grepl("Row names", msgs)])) cat("   message:", m, "\n")
}
