# Cases of the analyses of trait variation among individuals that the other test files do not
# reach: they fix the behaviour of the functions split into helpers (goodpractice).

refactor_sim <- function(record_steps = NULL) {
  set.seed(11)
  simulate_metacommunity(landscape_implicit(8), species_traits(6),
                         competition_matrix(6, "none"), n_steps = 40, seed_steps = 1:5,
                         record_steps = record_steps)
}

test_that("character displacement needs a positive niche breadth or a width", {
  sim <- refactor_sim()
  sim$traits$sigma <- 0
  expect_error(sample_individuals(sim, displacement = 1), "give `displacement_width`")
  sim$traits$sigma <- NA_real_
  expect_error(sample_individuals(sim, displacement = 1), "give `displacement_width`")
  ind <- sample_individuals(sim, displacement = 1, displacement_width = 2)
  expect_s3_class(ind, "data.frame")
  expect_gt(nrow(ind), 0)
})

test_that("sample_individuals() keeps the design of dates, visits and sessions", {
  sim <- refactor_sim(record_steps = c(20, 40))
  env <- matrix(seq_len(16) / 4, 8, 2, dimnames = list(rownames(sim$abundance), NULL))
  set.seed(3)
  ind <- sample_individuals(sim, n = 3, step = c(20, 40), env = env, visits = 2,
                            population_sd = 0.5, session_sd = 0.2, session = "visit",
                            cover = 0.9, min_n = 2)
  expect_named(ind, c("step", "visit", "site", "species", "individual", "abundance", "env",
                      "optimum", "population_mean", "trait"))
  expect_identical(sort(unique(ind$step)), c(20, 40))
  expect_identical(sort(unique(ind$visit)), 1:2)
  expect_identical(ind$env, env[cbind(match(ind$site, rownames(env)),
                                      match(ind$step, c(20, 40)))])
  counts <- table(paste(ind$step, ind$visit, ind$site, ind$species))
  expect_true(all(counts >= 2 & counts <= 3))
  set.seed(3)
  again <- sample_individuals(sim, n = 3, step = c(20, 40), env = env, visits = 2,
                              population_sd = 0.5, session_sd = 0.2, session = "visit",
                              cover = 0.9, min_n = 2)
  expect_identical(again, ind)
  one <- sample_individuals(sim, fraction = 0.5, step = 20, env = env[, 1, drop = FALSE])
  expect_identical(one$env, unname(env[one$site, 1]))
})

test_that("community_index() names a single function after the call", {
  set.seed(5)
  ind <- data.frame(site = rep(c("a", "b", "c"), each = 12),
                    species = rep(c("x", "y", "z"), 12), trait = rnorm(36))
  named <- community_index(ind, stats::var, n_perm = 0)
  expect_identical(unique(named$index), "stats::var")
  inline <- community_index(ind, function(x) diff(range(x)), n_perm = 0)
  expect_identical(unique(inline$index), "index")
  expect_error(community_index(ind, list(a = stats::var, a = stats::sd)), "distinct names")
  expect_error(community_index(ind, list(stats::var)), "distinct names")
  set.seed(6)
  pops <- community_index(ind, list(v = stats::var, m = mean), null = "regional.pop.prab",
                          n_perm = 9)
  expect_identical(dim(attr(pops, "null")), c(9L, 6L))
  expect_identical(pops$index, rep(c("v", "m"), each = 3))
})

test_that("drift_test() takes W as a named data frame in any order", {
  set.seed(7)
  group <- rep(1:12, each = 5)
  x <- matrix(rnorm(48, sd = 3), 12, 4)[group, ] + matrix(rnorm(240), 60, 4)
  colnames(x) <- paste0("t", 1:4)
  w <- diag(c(4, 3, 2, 1))
  dimnames(w) <- list(colnames(x), colnames(x))
  direct <- drift_test(x, group, W = w)
  shuffled <- drift_test(x, group, W = as.data.frame(w[4:1, 4:1]))
  expect_identical(shuffled$W, direct$W)
  expect_identical(shuffled$slope, direct$slope)
  expect_identical(shuffled$components, direct$components)
  expect_error(drift_test(x, group, W = w[, 4:1]), "symmetric matrix")
  expect_error(drift_test(x, group, W = diag(3)), "symmetric matrix")
  expect_error(drift_test(x, group, W = -w), "positive definite")
})

test_that("plot() of null-model tests stops on values that are not in the table", {
  skip_if_not_installed("ggplot2")
  set.seed(8)
  ind <- data.frame(site = rep(paste0("s", 1:4), each = 30),
                    species = sample(c("a", "b", "c"), 120, replace = TRUE))
  ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(120, sd = 3)
  res <- tstats(ind, n_perm = 19)
  expect_error(plot(res, sites = "s9"), "`sites` has values not in `x`: s9")
  expect_error(plot(res, statistics = "T"), "`statistics` has values not in `x`: T")
  p <- plot(res, type = "ses", hatch = FALSE, sites = c("s1", "s2"), order = "observed")
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 4)
  expect_identical(levels(droplevels(ggplot2::layer_data(p, 3)$PANEL)), c("1", "2", "3"))
})
