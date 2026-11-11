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

test_that("abc_choose tells apart distinct models and not identical ones", {
  skip_if_not_installed("abcrf")
  set.seed(5)
  a <- data.frame(s1 = rnorm(200, 0), s2 = rnorm(200), only_a = 1)
  b <- data.frame(s1 = rnorm(200, 3), s2 = rnorm(200))
  choice <- abc_choose(list(neutral = a, niche = b), ntree = 100)
  expect_s3_class(choice, "mc_abc_choice")
  expect_identical(choice$statistics, c("s1", "s2"))
  expect_lt(choice$error, 0.1)
  expect_equal(dim(choice$confusion), c(2, 2))
  post <- predict(choice, data.frame(s1 = c(-0.5, 3.5), s2 = 0), ntree = 100)
  expect_equal(post$model, c("neutral", "niche"))
  expect_true(all(post$probability > 0.5 & post$probability <= 1))
  expect_equal(post$votes_neutral + post$votes_niche, c(100, 100))
  expect_output(print(choice), "out-of-bag error")
  same <- abc_choose(list(x = a[1:2], y = data.frame(s1 = rnorm(200), s2 = rnorm(200))),
                     ntree = 100)
  expect_gt(same$error, 0.3)
  long <- rbind(cbind(model = "neutral", a[1:2]), cbind(model = "niche", b))
  expect_identical(abc_choose(long, ntree = 50)$models, c("neutral", "niche"))
})

test_that("abc_choose checks its input", {
  skip_if_not_installed("abcrf")
  a <- data.frame(s1 = rnorm(20), s2 = 1)
  expect_error(abc_choose(list(a)), "named list")
  expect_error(abc_choose(a), "column `model`")
  expect_error(abc_choose(list(x = a, y = a), statistics = "zz"), "Not in every")
  expect_error(abc_choose(list(x = a, y = a), statistics = "s2"), "No statistic")
  expect_error(abc_choose(list(x = a[1:5, ], y = a)), "Fewer than 10")
})

test_that("abc_ppc flags a statistic the model does not reproduce", {
  skip_if_not_installed("abcrf")
  set.seed(6)
  table <- data.frame(a = runif(600))
  table$s1 <- table$a + rnorm(600, 0, 0.05)
  table$s2 <- table$a^2 + rnorm(600, 0, 0.02)
  fit <- abc_fit(table, "a", "s1", ntree = 100)
  good <- abc_ppc(fit, c(s1 = 0.5, s2 = 0.25), table, n = 200)
  bad <- abc_ppc(fit, c(s1 = 0.5, s2 = 0.9), table, n = 200)
  expect_identical(good$statistic, "s2")
  expect_gt(good$p, 0.2)
  expect_lt(bad$p, 0.05)
  draws <- attr(good, "draws")
  expect_equal(nrow(draws), 200)
  expect_lt(abs(stats::median(draws$a) - 0.5), 0.1)
  all <- abc_ppc(fit, c(s1 = 0.5), table[c("a", "s1")], n = 50)
  expect_identical(all$statistic, "s1")
  sim <- abc_ppc(fit, c(s1 = 0.5, s2 = 0.25), table, n = 20,
                 simulate = function(scenario, shared) scenario$a,
                 summarise = function(sim, shared) c(s2 = sim^2))
  expect_equal(attr(sim, "draws")$s2, attr(sim, "draws")$a^2)
  expect_error(abc_ppc(fit, c(s1 = 0.5, s2 = 0.2), table, statistics = "zz"), "Not in")
  expect_error(abc_ppc(fit, c(s1 = 0.5, s2 = 0.2), table[-1, ]), "not the table")
})
