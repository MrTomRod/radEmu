suppressMessages(library(radEmu))
set.seed(2)
worst <- 0; wr <- 0
for (i in 1:2000) {
  x <- rnorm(sample(5:60, 1), sd = runif(1, 0.3, 3))
  a <- pseudohuber_median(x); b <- pseudohuber_median(x, tolerance = 1e-14)
  if (abs(a - b) > worst) { worst <- abs(a - b); wr <- worst / max(abs(b), 1e-3) }
}
cat(sprintf("2000 random vectors: largest |default - converged| = %.2e (relative to max(|value|, 1e-3): %.1e)\n", worst, wr))
