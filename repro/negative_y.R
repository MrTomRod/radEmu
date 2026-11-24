# emuFit accepts negative Y silently (documented as nonnegative).
suppressPackageStartupMessages(library(radEmu))
set.seed(1)
n <- 10; J <- 6
Y <- matrix(rpois(n * J, 20), n, J, dimnames = list(paste0("s", 1:n), paste0("t", 1:J)))
X <- cbind(intercept = 1, group = rep(0:1, each = 5)); rownames(X) <- rownames(Y)
fit <- function(Y) emuFit(Y = Y, X = X, run_score_tests = FALSE, verbose = FALSE)$coef$estimate
cat("original:     ", round(fit(Y), 3), "\n")
Y[2, 3] <- -30
cat("Y[2,3] = -30: ", round(fit(Y), 3), " (no error or warning)\n")
