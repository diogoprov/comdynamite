calibration_data <- function(n) {
  d <- data.frame(a = stats::runif(n), b = stats::runif(n))
  d$s <- d$a + stats::rnorm(n, 0, 0.05)
  d$t <- d$b + stats::rnorm(n, 0, 0.3)
  d
}

test_that("abc_calibrate() measures coverage and widens intervals that are too narrow", {
  skip_if_not_installed("abcrf")
  set.seed(11)
  # with one statistic per parameter and small leaves the nominal intervals are too narrow
  fit <- abc_fit(calibration_data(2000), parameters = c("a", "b"), statistics = c("s", "t"),
                 ntree = 100)
  cal <- abc_calibrate(fit, calibration_data(300))
  expect_s3_class(cal, "mc_abc_calibration")
  expect_named(cal, c("pit", "coverage", "n"))
  expect_identical(cal$n, 300L)
  expect_named(cal$pit, c("a", "b"))
  expect_true(all(cal$pit$a > 0 & cal$pit$a < 1))
  expect_identical(cal$coverage$level, rep(c(0.5, 0.8, 0.9, 0.95), 2))
  test <- calibration_data(300)
  covered <- function(est) mean(est$truth >= est$low & est$truth <= est$high)
  plain <- predict(fit, test, level = 0.9)
  recal <- predict(fit, test, level = 0.9, calibration = cal)
  plain$truth <- ifelse(plain$parameter == "a", test$a[plain$observation],
                        test$b[plain$observation])
  recal$truth <- plain$truth
  expect_lt(covered(plain), 0.85)
  expect_gt(covered(recal), covered(plain))
  expect_lt(abs(covered(recal) - 0.9), 0.05)
  # the estimate does not change, only the interval
  expect_identical(recal$estimate, plain$estimate)
  expect_gt(mean(recal$high - recal$low), mean(plain$high - plain$low))
})

test_that("calibrated probabilities are the quantiles of the PIT, kept inside the grid", {
  pit <- seq(0.3, 0.7, length.out = 101)
  expect_equal(calibrated_probs(pit, c(0.05, 0.95)),
               unname(stats::quantile(pit, c(0.05, 0.95), type = 8)))
  expect_identical(calibrated_probs(c(0, 0, 1, 1), c(0.05, 0.95)), c(0.001, 0.999))
})

test_that("abc_calibrate() and predict() check their arguments", {
  skip_if_not_installed("abcrf")
  set.seed(12)
  fit <- abc_fit(calibration_data(200), parameters = c("a", "b"), statistics = c("s", "t"),
                 ntree = 20)
  pseudo <- calibration_data(30)
  expect_error(abc_calibrate(fit, as.list(pseudo)), "data frame")
  expect_error(abc_calibrate(fit, pseudo[c("a", "s", "t")]), "Not in `pseudo`: b")
  expect_error(abc_calibrate(fit, pseudo[1:19, ]), "at least 20")
  pseudo$s[1:15] <- NA
  expect_error(abc_calibrate(fit, pseudo), "at least 20")
  expect_error(predict(fit, c(s = 0.5, t = 0.5), calibration = list()), "abc_calibrate")
  cal <- abc_calibrate(fit, calibration_data(30))
  cal$pit$b <- NULL
  expect_error(predict(fit, c(s = 0.5, t = 0.5), calibration = cal), "lacks the parameters b")
})
