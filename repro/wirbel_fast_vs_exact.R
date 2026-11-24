# Comparison on wirbel_otu_small (n = 320, J = 47), test of CRC for every taxon:
# fastEmuFit (default reference_set_size = 30) vs exact emuFit (default settings).
suppressMessages({ library(radEmu); library(fastEmu) })
data("wirbel_otu_small", package = "radEmu")
data("wirbel_sample_small", package = "radEmu")
Y <- as.matrix(wirbel_otu_small)
storage.mode(Y) <- "double"
X <- cbind(`(Intercept)` = 1, CRC = as.numeric(wirbel_sample_small$Group == "CRC"))
rownames(Y) <- rownames(X) <- sprintf("s%03d", seq_len(nrow(Y)))
tk <- data.frame(k = 2, j = seq_len(ncol(Y)))

fit_exact <- suppressMessages(emuFit(Y = Y, X = X, test_kj = tk, verbose = FALSE))
set.seed(77)
fit_fast <- suppressMessages(fastEmuFit(Y = Y, X = X, test_kj = tk, verbose = FALSE))

p_exact <- fit_exact$coef$pval[fit_exact$coef$covariate == "CRC"]
p_fast <- fit_fast$coef$pval[fit_fast$coef$covariate == "CRC"]
d <- abs(log10(p_fast) - log10(p_exact))
cat(sprintf("taxa %d | Spearman of p %.3f | |dlog10 p|: median %.2f, max %.2f, >0.1 in %d taxa | p<0.05: fast %d, exact %d\n",
            length(p_exact), cor(p_fast, p_exact, method = "spearman"), median(d), max(d), sum(d > 0.1),
            sum(p_fast < 0.05), sum(p_exact < 0.05)))
