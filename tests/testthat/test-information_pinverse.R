# Singular Fisher information: information_treatment = "identifiable" (default) vs "r_ladder" (radEmu <= 2.4.0.0)

make_absent <- function(b_absent = -45, n = 14, J = 7, seed = 7001) {
  set.seed(seed)
  X <- cbind(1, as.numeric(seq_len(n) > n / 2))
  Y <- 0.5 + 9 * matrix(runif(n * J), n, J)
  Y[X[, 2] == 1, 4] <- 0          # category 4 is absent from the second group
  B <- rbind(0.3 * sin(0:(J - 1)), 0.2 * cos((0:(J - 1)) * 1.7))
  B[2, 4] <- b_absent
  list(Y = Y, X = X, B = B, J = J, p = 2, n = n)
}

score_at <- function(d, k, j, cons_fn, cons_grad, treatment, j_ref = d$J) {
  get_score_stat(Y = d$Y, X_cup = X_cup_from_X(d$X, d$J), X = d$X, B = d$B, k_constr = k, j_constr = j,
                 constraint_grad_fn = cons_grad, indexes_to_remove = (j_ref - 1) * d$p + 1:d$p, j_ref = j_ref,
                 J = d$J, n = d$n, p = d$p, information_treatment = treatment)
}

test_that("regular information: no directions removed, identifiable == r_ladder exactly", {
  d <- make_absent(b_absent = 0.1)
  d$Y[d$X[, 2] == 1, 4] <- 1
  I <- f_info(d$Y, B_cup_from_B(d$B - d$B[, d$J]), d$B - d$B[, d$J], d$X, X_cup_from_X(d$X, d$J))
  expect_null(information_pinverse(I, d$p, d$J, d$J))
  a <- score_at(d, 2, 1, pseudohuber_median, dpseudohuber_median_dx, "identifiable")
  b <- score_at(d, 2, 1, pseudohuber_median, dpseudohuber_median_dx, "r_ladder")
  expect_identical(a$score_stat, b$score_stat)
})

test_that("a category absent from a group: removed direction, statistic independent of the reference", {
  d <- make_absent()
  I <- f_info(d$Y, B_cup_from_B(d$B - d$B[, d$J]), d$B - d$B[, d$J], d$X, X_cup_from_X(d$X, d$J))
  pinv <- information_pinverse(I, d$p, d$J, d$J)
  expect_false(is.null(pinv))
  expect_equal(pinv$dropped$taxon, 4)
  s_ref7 <- score_at(d, 2, 1, pseudohuber_median, dpseudohuber_median_dx, "identifiable", j_ref = 7)$score_stat
  s_ref1 <- score_at(d, 2, 1, pseudohuber_median, dpseudohuber_median_dx, "identifiable", j_ref = 1)$score_stat
  s_ref5 <- score_at(d, 2, 1, pseudohuber_median, dpseudohuber_median_dx, "identifiable", j_ref = 5)$score_stat
  expect_true(is.finite(s_ref7))
  # equal up to the size of the part of H in the removed direction (about 4e-9 here)
  expect_equal(s_ref1, s_ref7, tolerance = 1e-6)
  expect_equal(s_ref5, s_ref7, tolerance = 1e-6)
  # R's ladder gives a (slightly) different number
  s_old <- score_at(d, 2, 1, pseudohuber_median, dpseudohuber_median_dx, "r_ladder")$score_stat
  expect_false(isTRUE(all.equal(s_old, s_ref7, tolerance = 1e-8)))
})

test_that("a hypothesis that depends on the unidentified parameter is not testable", {
  d <- make_absent()
  mean_grad <- function(x) rep(1 / length(x), length(x))
  expect_true(is.na(score_at(d, 2, 1, mean, mean_grad, "identifiable")$score_stat))
  expect_true(is.finite(score_at(d, 2, 1, mean, mean_grad, "r_ladder")$score_stat) ||
              inherits(try(score_at(d, 2, 1, mean, mean_grad, "r_ladder"), silent = TRUE), "try-error"))
})

test_that("Wald standard errors: NA for the unidentified category, finite for the others, error with r_ladder", {
  d <- make_absent()
  sub <- c(1, 2, 3, 5, 6)
  fn <- list(function(x) mean(x[sub]), function(x) mean(x[sub]))
  gr <- function(x) { g <- rep(0, length(x)); g[sub] <- 1 / length(sub); g }
  grl <- list(gr, gr)
  tk <- data.frame(k = 2, j = 1:7)
  args <- list(Y = d$Y, X = d$X, X_cup = X_cup_from_X(d$X, d$J), B = d$B, test_kj = tk, constraint_fn = fn,
               constraint_grad_fn = grl, j_ref = 7)
  w <- do.call(micro_wald, c(args, list(information_treatment = "identifiable")))$coefficients
  expect_true(is.na(w$se[4]))
  expect_true(all(is.finite(w$se[-4])) && all(w$se[-4] > 0))
  expect_error(do.call(micro_wald, c(args, list(information_treatment = "r_ladder"))))
})
