toy_table <- function(n = 400) {
  set.seed(3)
  table <- data.frame(a = runif(n), b = runif(n))
  table$s1 <- table$a + rnorm(n, 0, 0.05)
  table$s2 <- rnorm(n)
  table$flat <- 1
  table
}

test_that("abc_table draws within the prior and keeps parameters and statistics", {
  skip_if_not_installed("treeshapesim")
  shared_pool <- function() {
    list(tree = treeshapesim::sim_tree_shape(8), landscape = landscape_implicit(6))
  }
  one_run <- function(scenario, shared) {
    simulate_metacommunity(shared$landscape, species_traits(shared$tree, scenario$lambda),
                           competition_matrix(shared$tree, "none"),
                           dispersal = 10^scenario$log_dispersal, n_steps = 20,
                           seed_steps = 1:3)
  }
  prior <- list(log_dispersal = c(-3, -1), lambda = c(0, 1))
  table <- abc_table(prior, one_run, n = 6, setup = shared_pool, seed = 4)
  expect_s3_class(table, "mc_abc_table")
  expect_equal(nrow(table), 6)
  expect_true(all(table$log_dispersal >= -3 & table$log_dispersal <= -1))
  expect_identical(attr(table, "parameters"), names(prior))
  expect_false("replicate" %in% names(table))
  expect_true("alpha" %in% attr(table, "statistics"))
  expect_equal(abc_table(prior, one_run, n = 6, setup = shared_pool, seed = 4), table)
})

test_that("abc_table checks the prior", {
  f <- function(scenario, shared) NULL
  expect_error(abc_table(list(c(0, 1)), f, 2), "named list")
  expect_error(abc_table(list(a = c(1, 0)), f, 2), "increasing")
  expect_error(abc_table(list(a = 1), f, 2), "pair")
})

test_that("abc_fit recovers an informative parameter and not an uninformed one", {
  skip_if_not_installed("abcrf")
  table <- toy_table()
  fit <- abc_fit(table, parameters = c("a", "b"), ntree = 100)
  expect_s3_class(fit, "mc_abc_fit")
  expect_false("flat" %in% fit$statistics)
  post <- predict(fit, c(s1 = 0.7, s2 = 0, flat = 1))
  expect_equal(post$parameter, c("a", "b"))
  expect_true(all(post$low <= post$estimate & post$estimate <= post$high))
  a <- post[post$parameter == "a", ]
  b <- post[post$parameter == "b", ]
  expect_lt(abs(a$estimate - 0.7), 0.1)
  expect_gt(b$high - b$low, 3 * (a$high - a$low))
  imp <- abc_importance(fit)
  expect_equal(imp$statistic[imp$parameter == "a"][1], "s1")
  expect_equal(sum(imp$importance[imp$parameter == "a"]), 1)
  expect_output(print(fit), "2 parameters")
  two <- predict(fit, data.frame(s1 = c(0.2, 0.8), s2 = 0))
  expect_equal(two$observation, c(1, 1, 2, 2))
  expect_error(predict(fit, c(s1 = 0.5)), "missing")
  expect_error(predict(fit, c(s1 = NA, s2 = 0)), "infinite")
})

test_that("abc_fit checks its input", {
  skip_if_not_installed("abcrf")
  table <- toy_table(50)
  expect_error(abc_fit(table), "parameters")
  expect_error(abc_fit(table, "a", "nope"), "Not in")
  expect_error(abc_fit(table, "a", "flat"), "No statistic")
  expect_error(abc_fit(table[1:8, ], "a", "s1"), "Fewer than 10")
  expect_error(abc_recovery(table, n_test = 5), "parameters")
})

test_that("abc_recovery separates a recoverable from an unrecoverable parameter", {
  skip_if_not_installed("abcrf")
  table <- toy_table()
  rec <- abc_recovery(table, parameters = c("a", "b"), n_test = 60, ntree = 100)
  expect_equal(rec$parameter, c("a", "b"))
  expect_gt(rec$r2[1], 0.8)
  expect_lt(rec$r2[2], 0.2)
  expect_gt(rec$coverage[1], 0.7)
  expect_lt(rec$width[1], rec$width[2])
  est <- attr(rec, "estimates")
  expect_equal(nrow(est), 120)
  rows <- match(est$truth[est$parameter == "a"], table$a)
  expect_false(anyNA(rows))
  expect_equal(est$truth[est$parameter == "b"], table$b[rows])
})
