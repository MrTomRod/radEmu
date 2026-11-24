suppressMessages({library(radEmu); library(fastEmu)})
set.seed(1)
n <- 20; J <- 8
Y <- matrix(rpois(n * J, 20), n, J, dimnames = list(NULL, paste0("cat", 1:J)))
X <- cbind(`(Intercept)` = 1, x = rep(0:1, each = n / 2))
show <- function(label, expr) {
  r <- tryCatch({ force(expr); "OK" }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(label, "->", r, "\n")
}
show("emuFit(compute_cis = TRUE,  return_wald_p = TRUE, run_score_tests = FALSE)",
     emuFit(Y = Y, X = X, compute_cis = TRUE,  return_wald_p = TRUE, run_score_tests = FALSE))
show("emuFit(compute_cis = FALSE, return_wald_p = TRUE, run_score_tests = TRUE) ",
     suppressMessages(emuFit(Y = Y, X = X, compute_cis = FALSE, return_wald_p = TRUE, run_score_tests = TRUE,
                             test_kj = data.frame(k = 2, j = 1))))
show("emuFit(compute_cis = FALSE, return_wald_p = TRUE, run_score_tests = FALSE)",
     emuFit(Y = Y, X = X, compute_cis = FALSE, return_wald_p = TRUE, run_score_tests = FALSE))
show("fastEmuFit(compute_cis = FALSE, return_wald_p = TRUE, run_score_tests = FALSE)",
     fastEmuFit(Y = Y, X = X, reference_set_size = 4, compute_cis = FALSE, return_wald_p = TRUE, run_score_tests = FALSE))
show("fastEmuFit(compute_cis = FALSE, return_wald_p = TRUE) [score tests on, default]",
     suppressMessages(fastEmuFit(Y = Y, X = X, reference_set_size = 4, compute_cis = FALSE, return_wald_p = TRUE,
                                 test_kj = data.frame(k = 2, j = 1))))
