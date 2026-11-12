test_that("spatial_folds makes contiguous blocks or balanced random folds", {
  grid <- landscape_grid(side = 6)
  blocks <- spatial_folds(grid$coords, k = 4)
  expect_length(blocks, 36)
  expect_setequal(blocks, 1:4)
  expect_identical(spatial_folds(grid$coords, k = 4), blocks)
  # blocks are compact: mean distance within blocks is smaller than at random
  within <- function(f) {
    mean(unlist(lapply(split(seq_along(f), f), function(s) stats::dist(grid$coords[s, ]))))
  }
  random <- spatial_folds(grid$coords, k = 4, method = "random")
  expect_equal(as.vector(table(random)), rep(9, 4))
  expect_lt(within(blocks), within(random))
  expect_error(spatial_folds(grid$coords, k = 1))
})

test_that("fold_summarise gives the statistics of kept and left-out sites", {
  skip_if_not_installed("treeshapesim")
  set.seed(1)
  tree <- treeshapesim::sim_tree_shape(8)
  grid <- landscape_grid(side = 4)
  sim <- simulate_metacommunity(grid, species_traits(tree, 0.5), competition_matrix(tree, "none"),
                                dispersal = 0.05, n_steps = 40)
  folds <- spatial_folds(grid$coords, k = 2)
  stat <- function(x, shared) c(sites = nrow(x$comm), alpha = mean(rowSums(x$comm > 0)))
  by_fold <- fold_summarise(folds, stat)
  out <- by_fold(sim, list())
  expect_named(out, c("f1_train_sites", "f1_train_alpha", "f1_test_sites", "f1_test_alpha",
                      "f2_train_sites", "f2_train_alpha", "f2_test_sites", "f2_test_alpha"))
  expect_equal(unname(out["f1_test_sites"]), sum(folds == 1))
  richness <- rowSums(sim$abundance > 0)
  expect_equal(unname(out["f2_test_alpha"]), mean(richness[folds == 2]))
  expect_identical(by_fold(as_mc_data(sim), list()), out)
  expect_error(fold_summarise(folds[1:3], stat)(sim, list()), "one value per site")
  expect_error(fold_summarise(rep(1, 16)), "two folds")
})

cv_table <- function(n = 500) {
  set.seed(7)
  table <- data.frame(a = runif(n))
  for (f in 1:3) {
    table[[paste0("f", f, "_train_s")]] <- table$a + rnorm(n, 0, 0.05)
    table[[paste0("f", f, "_test_s")]] <- table$a + rnorm(n, 0, 0.1)
    table[[paste0("f", f, "_test_t")]] <- rnorm(n)
  }
  table
}

test_that("abc_cv predicts left-out statistics and flags those the model misses", {
  skip_if_not_installed("abcrf")
  table <- cv_table()
  observed <- c(f1_train_s = 0.5, f1_test_s = 0.52, f1_test_t = 0,
                f2_train_s = 0.5, f2_test_s = 0.48, f2_test_t = 0.1,
                f3_train_s = 0.5, f3_test_s = 1.4, f3_test_t = -0.2)
  cv <- abc_cv(table, observed, parameters = "a", n = 300, ntree = 100)
  expect_s3_class(cv, "mc_abc_cv")
  expect_equal(nrow(cv), 6)
  s <- cv[cv$statistic == "s", ]
  expect_equal(s$inside, c(TRUE, TRUE, FALSE))
  expect_lt(s$p[3], 0.05)
  expect_gt(s$z[3], 3)
  expect_lt(abs(s$median[1] - 0.5), 0.1)
  sm <- summary(cv)
  expect_equal(sm$statistic, c("s", "t"))
  expect_equal(sm$coverage[1], 2 / 3)
  expect_equal(sm$rejected[1], 1 / 3)
  two <- abc_cv(table, rbind(as.data.frame(as.list(observed)),
                             as.data.frame(as.list(observed))), "a", n = 300, ntree = 100)
  expect_equal(nrow(two), 12)
  expect_equal(two$inside[two$observation == 2], cv$inside)
  only_s <- abc_cv(table, observed, parameters = "a", statistics = "s", n = 100, ntree = 50)
  expect_setequal(only_s$statistic, "s")
})

test_that("abc_cv checks its input", {
  skip_if_not_installed("abcrf")
  table <- cv_table(60)
  observed <- c(f1_train_s = 0.5)
  expect_error(abc_cv(table, observed), "parameters")
  expect_error(abc_cv(data.frame(a = 1:20, s = 1:20), c(s = 1), "a"), "no statistics by fold")
  expect_error(abc_cv(table, observed, "a", statistics = "zz"), "None of")
  expect_error(abc_cv(table, observed, "a"), "missing from `observed`")
})
