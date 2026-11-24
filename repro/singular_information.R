# Unpenalised fit, two groups, one taxon absent from group 2: the information matrix is singular.
suppressPackageStartupMessages(library(radEmu))
set.seed(3101)
n <- 12; J <- 160
X <- cbind(intercept = 1, group = rep(0:1, each = n / 2))
Y <- matrix(rpois(n * J, 2), n, J, dimnames = list(paste0("s", 1:n), paste0("t", 1:J)))
rownames(X) <- rownames(Y)
Y[, 1] <- Y[, 1] + 5                             # a clearly most-prevalent taxon
Y[X[, "group"] == 1, 3] <- 0                     # taxon t3 never seen in group 2

fit <- emuFit(Y = Y, X = X, penalize = FALSE, run_score_tests = FALSE, compute_cis = FALSE, verbose = FALSE)
cat("estimate for t3 (group):", round(fit$B[2, 3], 1), "\n")

# Wald confidence intervals
r <- tryCatch(emuFit(Y = Y, X = X, penalize = FALSE, run_score_tests = FALSE, compute_cis = TRUE,
                     verbose = FALSE),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("Wald:", if (is.character(r)) r else sprintf("OK; se for t3 = %s, for t5 = %.4f", r$coef$se[3], r$coef$se[5]), "\n")

# Score test of taxon t1: watch which ridge lambda get_score_stat ends on
lam <- NA
invisible(trace("get_score_stat", exit = quote(lam <<- lambda), where = asNamespace("radEmu"), print = FALSE))
s <- emuFit(Y = Y, X = X, penalize = FALSE, compute_cis = FALSE,
            test_kj = data.frame(k = 2, j = 1), verbose = FALSE)
cat("score test t1: stat", signif(s$coef$score_stat[1], 6), " ridge lambda used:", lam, "\n")
