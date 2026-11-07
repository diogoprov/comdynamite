pool <- function() list(tree = ape::rcoal(6), landscape = landscape_implicit(5))
one_run <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape, species_traits(shared$tree, 1),
                         competition_matrix(shared$tree, "none"),
                         dispersal = scenario$dispersal, n_steps = 15, seed_steps = 1:3)
}
grid <- data.frame(dispersal = c(0.01, 0.1, 0.3))

test_that("run_scenarios() pairs replicates and repeats with a seed", {
  note <- function(sim, shared) {
    c(tips = sum(shared$tree$edge.length), alpha = mean(rowSums(sim$abundance > 0)))
  }
  out <- run_scenarios(grid, one_run, pool, note, n_rep = 2, seed = 5)
  expect_equal(nrow(out), 6)
  expect_named(out, c("dispersal", "replicate", "tips", "alpha"))
  expect_equal(out$dispersal, rep(grid$dispersal, 2))
  expect_equal(length(unique(out$tips[out$replicate == 1])), 1)
  expect_false(out$tips[1] == out$tips[4])
  expect_identical(out, run_scenarios(grid, one_run, pool, note, n_rep = 2, seed = 5))
})

test_that("run_scenarios() keeps the generator of the session and works without setup", {
  set.seed(9)
  before <- .Random.seed
  out <- run_scenarios(grid[1, , drop = FALSE], function(scenario, shared) {
    simulate_static(species_traits(5), seq(0, 100, length.out = 6))
  })
  expect_identical(.Random.seed, before)
  expect_true("alpha" %in% names(out))
  rm(".Random.seed", envir = globalenv())
  run_scenarios(grid[1, , drop = FALSE], function(scenario, shared) {
    simulate_static(species_traits(5), seq(0, 100, length.out = 6))
  })
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(1)
})

test_that("run_scenarios() gives the same results under a future plan", {
  skip_if_not_installed("future.apply")
  skip_on_os("windows")
  serial <- run_scenarios(grid, one_run, pool, n_rep = 2)
  future::plan(future::multicore, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  expect_equal(run_scenarios(grid, one_run, pool, n_rep = 2), serial)
})

sim_with <- function(history) {
  structure(list(history = history, abundance = history[, , dim(history)[3]]), class = "mc_sim")
}
series <- function(values) {
  array(rep(values, each = 6), c(3, 2, length(values)),
        dimnames = list(paste0("site", 1:3), paste0("sp", 1:2), paste0("step", seq_along(values))))
}

test_that("equilibrium_check() tells settled from changing runs", {
  still <- equilibrium_check(sim_with(series(rep(5, 6))))
  expect_true(all(still$settled))
  expect_equal(still$statistic, c("gamma", "alpha", "abundance"))
  expect_equal(still$early, c(2, 2, 10))
  falling <- series(c(10, 10, 10, 5, 5, 5))
  falling[, 2, 4:6] <- 0
  expect_warning(out <- equilibrium_check(sim_with(falling)), "gamma, alpha, abundance")
  expect_equal(out$change, c(-0.5, -0.5, -0.75))
  expect_true(all(equilibrium_check(sim_with(series(rep(0, 4))))$settled))
  expect_error(equilibrium_check(sim_with(series(1:3))), "four recorded steps")
  expect_error(equilibrium_check(structure(list(history = NULL), class = "mc_sim")), "record_steps")
})

test_that("as.data.frame() gives a long table", {
  history <- series(c(0, 3))
  history["site2", "sp1", 1] <- 7L
  x <- sim_with(history)
  long <- as.data.frame(x)
  expect_named(long, c("site", "species", "step", "abundance"))
  expect_equal(nrow(long), 7)
  expect_equal(long$abundance[long$site == "site2" & long$species == "sp1" & long$step == 1], 7)
  expect_equal(nrow(as.data.frame(x, zeros = TRUE)), 12)
  no_history <- structure(list(history = NULL, abundance = history[, , 2]), class = "mc_sim")
  final <- as.data.frame(no_history)
  expect_true(all(is.na(final$step)))
  expect_equal(sum(final$abundance), sum(history[, , 2]))
})
