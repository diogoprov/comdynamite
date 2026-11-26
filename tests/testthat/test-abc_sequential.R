toy_sim <- function(scenario, shared) {
  c(s = scenario$a + stats::rnorm(1, 0, 0.1), t = scenario$b + stats::rnorm(1, 0, 1))
}
toy_summary <- function(sim, shared) sim
toy_prior <- list(a = c(-10, 10), b = c(0, 5))

test_that("abc_sequential() narrows a wide prior around the posterior", {
  skip_if_not_installed("abcrf")
  # few effective simulations per round: the warning says the intervals may be too narrow
  expect_warning(
    fit <- abc_sequential(toy_prior, toy_sim, observed = c(s = 2, t = 2.5), n = 300,
                          rounds = 3, summarise = toy_summary, ntree = 100, seed = 3),
    "Fewer than 50 effective simulations"
  )
  expect_s3_class(fit, "mc_abc_sequential")
  expect_length(fit$rounds, 3)
  expect_named(fit$rounds[[2]], c("table", "fit", "weights"))
  expect_named(fit$rounds[[2]]$weights, c("a", "b", "w_a", "w_b"))
  expect_equal(sum(fit$rounds[[3]]$weights$w_a), 1)
  # later rounds simulate inside the prior and near the posterior
  later <- fit$rounds[[3]]$table$a
  expect_true(all(later >= -10 & later <= 10))
  expect_lt(stats::sd(later), stats::sd(fit$rounds[[1]]$table$a) / 5)
  first <- predict(fit, round = 1)
  last <- predict(fit)
  expect_identical(last$round, c(3L, 3L))
  expect_named(last, c("round", "parameter", "estimate", "low", "high", "mean", "ess"))
  width <- function(e, p) (e$high - e$low)[e$parameter == p]
  expect_lt(width(last, "a"), width(first, "a"))
  # the exact posterior of a is about normal, with sd 0.1, around 2
  expect_lt(abs(last$estimate[1] - 2), 0.15)
  expect_lt(width(last, "a"), 0.6)
  expect_true(last$low[1] < 2 && last$high[1] > 2)
  expect_output(print(fit), "3 rounds")
  expect_output(print(fit), "Effective simulations per round")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(fit), "ggplot")
})

test_that("abc_sequential() takes a table from the prior as its first round", {
  skip_if_not_installed("abcrf")
  table <- abc_table(toy_prior, toy_sim, n = 150, summarise = toy_summary, seed = 4)
  fit <- suppressWarnings(
    abc_sequential(toy_prior, toy_sim, observed = data.frame(s = 1, t = 1), n = 100,
                   rounds = 2, summarise = toy_summary, table = table, ntree = 50)
  )
  expect_identical(fit$rounds[[1]]$table, table)
  expect_identical(nrow(fit$rounds[[2]]$table), 100L)
  one <- abc_sequential(toy_prior, toy_sim, observed = c(s = 1, t = 1), n = 100, rounds = 1,
                        summarise = toy_summary, table = table, ntree = 50)
  expect_equal(predict(one), predict(fit, round = 1))
})

test_that("setup is called once and shared by all rounds", {
  skip_if_not_installed("abcrf")
  calls <- 0
  setup <- function() {
    calls <<- calls + 1
    list(shift = stats::runif(1))
  }
  sim <- function(scenario, shared) c(s = scenario$a + shared$shift + stats::rnorm(1, 0, 0.1))
  fit <- suppressWarnings(
    abc_sequential(list(a = c(0, 1)), sim, observed = c(s = 1), n = 50, rounds = 2,
                   setup = setup, summarise = toy_summary, ntree = 20)
  )
  expect_identical(calls, 1)
})

test_that("weighted quantiles and the proposal density are right", {
  expect_identical(weighted_quantile(c(3, 1, 2), c(0.2, 0.3, 0.5), c(0.3, 0.5, 0.9)),
                   c(1, 2, 3))
  expect_equal(kernel_sd(c(0, 2), c(0.5, 0.5)), sqrt(2))
  expect_equal(sequential_density(1, c(0, 2), c(0.5, 0.5)), stats::dnorm(1, 0, sqrt(2)))
  set.seed(5)
  previous <- data.frame(a = c(0.2, 0.4), w_a = c(0.5, 0.5))
  draws <- sequential_proposal(previous, list(a = c(0, 0.5)), 500)
  expect_true(all(draws$a >= 0 & draws$a <= 0.5))
})

test_that("abc_sequential() and predict() check their arguments", {
  skip_if_not_installed("abcrf")
  expect_error(abc_sequential(toy_prior, toy_sim, observed = c(2, 2.5), n = 50),
               "named numeric")
  expect_error(abc_sequential(toy_prior, toy_sim, observed = data.frame(s = 1:2, t = 1:2),
                              n = 50), "one observation")
  expect_error(abc_sequential(toy_prior, toy_sim, observed = c(s = 2, t = 2), n = 5),
               "10 or more")
  expect_error(abc_sequential(toy_prior, toy_sim, observed = c(s = 2, t = 2), n = 50,
                              table = data.frame(a = 1, s = 1)), "Not in `table`: b")
  fit <- abc_sequential(toy_prior, toy_sim, observed = c(s = 2, t = 2), n = 50, rounds = 1,
                        summarise = toy_summary, ntree = 20)
  expect_error(predict(fit, round = 2), "1 or less")
  expect_error(abc_sequential(toy_prior, toy_sim, observed = c(s = 2), n = 50, rounds = 1,
                              summarise = toy_summary, ntree = 20), "missing")
})
