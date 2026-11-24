suppressMessages({library(radEmu); library(fastEmu)})
set.seed(4)
n <- 20; J <- 10
Y <- matrix(rpois(n * J, 20), n, J, dimnames = list(NULL, paste0("cat", 1:J)))
X <- cbind(`(Intercept)` = 1, x = rep(0:1, each = n / 2))
r <- tryCatch(suppressMessages(fastEmuFit(Y = Y, X = X, constraint_fn = 3, reference_set = 3,
                                          test_kj = data.frame(k = 2, j = 1))), error = function(e) conditionMessage(e))
cat("fastEmuFit(constraint_fn = 3, reference_set = 3): ", if (is.character(r)) paste("ERROR:", r) else "OK", "\n")
