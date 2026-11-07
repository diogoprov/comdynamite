sp <- c("a", "b", "c")
r3 <- matrix(3, 6, 3, dimnames = list(NULL, sp))
alpha3 <- diag(0.05, 3)

test_that("without pool and initial, results are the same as before", {
  set.seed(4)
  before <- mc_dynamics(r3, alpha3, dispersal = 0.1, n_steps = 50)
  set.seed(4)
  after <- mc_dynamics(r3, alpha3, dispersal = 0.1, n_steps = 50, pool = NULL, initial = NULL)
  expect_identical(before$abundance, after$abundance)
  expect_null(after$pool)
})

test_that("species never occur outside their pool, in any recorded step", {
  pool <- matrix(TRUE, 6, 3, dimnames = list(NULL, sp))
  pool[1:3, "a"] <- FALSE
  pool[, "c"] <- FALSE
  set.seed(1)
  out <- mc_dynamics(r3, alpha3, dispersal = 0.3, n_steps = 60, record_steps = 1:60,
                     pool = pool * 1)
  expect_true(all(out$history[1:3, "a", ] == 0))
  expect_true(all(out$history[, "c", ] == 0))
  expect_true(all(out$abundance[4:6, "a"] > 0))
  expect_identical(unname(out$pool), unname(pool))
})

test_that("initial is the starting state", {
  initial <- matrix(c(5, 0, 2, 0, 7, 0), 2, 3, dimnames = list(NULL, sp))
  out <- mc_dynamics(matrix(1, 2, 3), matrix(0, 3, 3), dispersal = 0, n_steps = 1,
                     seed_steps = numeric(0), initial = initial, sites = c("s1", "s2"),
                     species = sp)
  expect_true(all(out$abundance[initial == 0] == 0))
  presences <- mc_dynamics(matrix(1, 2, 3), matrix(0, 3, 3), dispersal = 0, n_steps = 1,
                           seed_steps = numeric(0), initial = initial > 0,
                           species = sp)$abundance
  expect_true(all(presences[initial == 0] == 0))
  # with r = 1 and no competition, the mean of the next step is the initial abundance
  set.seed(3)
  means <- Reduce(`+`, lapply(1:400, function(i) {
    mc_dynamics(matrix(1, 2, 3), matrix(0, 3, 3), dispersal = 0, n_steps = 1,
                seed_steps = numeric(0), initial = initial, species = sp)$abundance
  })) / 400
  expect_equal(unname(means), unname(initial), tolerance = 0.1)
})

test_that("pool and initial are matched by name and checked", {
  ls <- landscape_implicit(4)
  traits <- species_traits(ape::rcoal(3, tip.label = sp))
  alpha <- matrix(diag(0.05, 3), 3, dimnames = list(sp, sp))
  pool <- matrix(1, 4, 3, dimnames = list(paste0("site", 4:1), rev(sp)))
  pool["site4", "b"] <- 0
  set.seed(5)
  sim <- simulate_metacommunity(ls, traits, alpha, n_steps = 30, pool = pool,
                                initial = pool[paste0("site", 1:4), sp])
  expect_false(sim$pool["site4", "b"])
  expect_equal(sim$abundance["site4", "b"], 0)
  expect_s3_class(as_mc_data(sim), "mc_data")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = matrix(1, 6, 3)), "species names")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = matrix(1, 6, 2, dimnames = list(NULL, sp[1:2]))),
               "no column for the species c")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = matrix(1, 5, 3, dimnames = list(NULL, sp))),
               "5 rows")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = matrix(1, 6, 3, dimnames = list(letters[1:6], sp))),
               "no row for the sites")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = matrix(2, 6, 3, dimnames = list(NULL, sp))),
               "logical or 0/1")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = matrix(NA, 6, 3, dimnames = list(NULL, sp))),
               "missing values")
  expect_error(mc_dynamics(r3, alpha3, 0, pool = list(1)), "matrix of sites x species")
  expect_error(mc_dynamics(r3, alpha3, 0, initial = matrix(-1, 6, 3, dimnames = list(NULL, sp))),
               "whole numbers")
  expect_error(mc_dynamics(r3, alpha3, 0, initial = matrix(1, 6, 3, dimnames = list(NULL, sp)),
                           pool = matrix(0, 6, 3, dimnames = list(NULL, sp))), "outside `pool`")
  expect_warning(mc_dynamics(r3, alpha3, 0, n_steps = 2,
                             pool = data.frame(a = rep(1, 6), b = 1, c = 1, d = 1)),
                 "not simulated, dropped: d")
})

test_that("a species without any site in the pool dies out", {
  pool <- matrix(TRUE, 6, 3, dimnames = list(NULL, sp))
  pool[, "c"] <- FALSE
  set.seed(6)
  out <- mc_dynamics(r3, alpha3, dispersal = 0.2, n_steps = 30, seed_steps = 1:5,
                     record_steps = 30, pool = pool)
  expect_true(all(out$abundance[, "c"] == 0))
  expect_true(all(out$abundance[, c("a", "b")] > 0))
})
