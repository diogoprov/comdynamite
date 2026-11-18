traits6 <- function() {
  data.frame(optimum = 50, sigma = 1000, r_max = 2, row.names = "sp1")
}

test_that("dispersal_time switches emigration on and off", {
  ls <- landscape_points(cbind(0:1, 0), env = c(50, 50), disp_scale = 1)
  init <- matrix(c(200L, 0L), 2, 1, dimnames = list(NULL, "sp1"))
  # no emigration in the first 10 steps, then the usual rate
  mult <- c(rep(0, 10), rep(1, 10))
  set.seed(1)
  sim <- simulate_metacommunity(ls, traits6(), matrix(0.01), dispersal = 0.2, n_steps = 20,
                                initial = init, seed_steps = numeric(0),
                                record_steps = c(10, 20), dispersal_time = mult)
  expect_equal(sim$history[2, 1, "step10"], 0)
  expect_gt(sim$history[2, 1, "step20"], 0)
  expect_equal(sim$dispersal_time, mult)
})

test_that("dispersal_time as a matrix gives probabilities per species and step", {
  ls <- landscape_points(cbind(0:1, 0), env = c(50, 50), disp_scale = 1)
  tr <- data.frame(optimum = c(50, 50), sigma = 1000, r_max = 2, row.names = c("sp1", "sp2"))
  init <- matrix(c(200L, 0L, 200L, 0L), 2, 2, dimnames = list(NULL, c("sp1", "sp2")))
  p <- cbind(rep(0, 15), rep(0.3, 15))
  set.seed(2)
  sim <- simulate_metacommunity(ls, tr, diag(0.01, 2), dispersal = 0.5, n_steps = 15,
                                initial = init, seed_steps = numeric(0), dispersal_time = p)
  expect_equal(sim$abundance[2, 1], 0)
  expect_gt(sim$abundance[2, 2], 0)
})

test_that("dispersal_by_step builds and checks probabilities", {
  expect_null(dispersal_by_step(NULL, 0.1, 5, 2))
  p <- dispersal_by_step(c(0, 1, 20, 1, 1), c(0.1, 0.2), 5, 2)
  expect_equal(dim(p), c(5L, 2L))
  expect_equal(p[3, ], c(1, 1))
  expect_equal(p[2, ], c(0.1, 0.2))
  expect_equal(dispersal_by_step(matrix(0.3, 6, 1), 0.1, 5, 2), matrix(0.3, 5, 2))
  expect_error(dispersal_by_step(c(1, 1), 0.1, 5, 2), "per step")
  expect_error(dispersal_by_step(c(1, -1, 1, 1, 1), 0.1, 5, 2), "per step")
  expect_error(dispersal_by_step(matrix(2, 5, 2), 0.1, 5, 2), "probabilities")
  expect_error(dispersal_by_step(matrix(0.1, 5, 3), 0.1, 5, 2), "probabilities")
})

test_that("kernels and kernel_step change the routes of dispersal through time", {
  ls <- landscape_points(cbind(0:2, 0), env = rep(50, 3), disp_scale = 1)
  # site 1 sends only to site 2 in regime "dry", only to site 3 in regime "flood"
  dry <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
  flood <- rbind(c(0, 0, 1), c(1, 0, 1), c(1, 0, 0))
  init <- matrix(c(300L, 0L, 0L), 3, 1, dimnames = list(NULL, "sp1"))
  set.seed(3)
  sim <- simulate_metacommunity(ls, traits6(), matrix(0.01), dispersal = 0.1, n_steps = 1,
                                initial = init, seed_steps = numeric(0),
                                kernels = list(dry = dry, flood = flood),
                                kernel_step = "flood", immigration = "multinomial")
  expect_equal(sim$abundance[2, 1], 0)
  expect_gt(sim$abundance[3, 1], 0)
  set.seed(3)
  sim <- simulate_metacommunity(ls, traits6(), matrix(0.01), dispersal = 0.1, n_steps = 1,
                                initial = init, seed_steps = numeric(0),
                                kernels = list(dry, flood), kernel_step = 1,
                                immigration = "multinomial")
  expect_gt(sim$abundance[2, 1], 0)
  expect_equal(sim$abundance[3, 1], 0)
})

test_that("kernels_by_step checks its arguments", {
  d <- diag(0, 3) + 1
  expect_null(kernels_by_step(NULL, NULL, 5, 3, d))
  expect_error(kernels_by_step(list(d), NULL, 5, 3, d), "go together")
  expect_error(kernels_by_step(list(d), rep(1, 5), 5, 3, NULL), "spatially explicit")
  expect_error(kernels_by_step(list(diag(2)), rep(1, 5), 5, 3, d), "3 x 3")
  expect_error(kernels_by_step(list(diag(3)), rep(1, 5), 5, 3, d), "destination")
  expect_error(kernels_by_step(list(d), rep(2, 5), 5, 3, d), "index or name")
  expect_error(kernels_by_step(list(d), rep(1, 4), 5, 3, d), "one value per step")
  expect_error(kernels_by_step(list(a = d), rep("b", 5), 5, 3, d), "names")
  k <- kernels_by_step(d, rep(1, 5), 5, 3, d)
  expect_equal(rowSums(k[[1]]), rep(1, 3))
  expect_equal(diag(k[[1]]), rep(0, 3))
})

test_that("kernels accept landscapes and reject other objects", {
  ls <- landscape_points(cbind(0:2, 0), env = rep(50, 3), disp_scale = 1)
  d <- diag(0, 3) + 1
  k <- kernels_by_step(list(ls, d), c(1, 2), 2, 3, d)
  expect_equal(k[[1]], unname(ls$dispersal))
  expect_error(kernels_by_step(list(), 1, 1, 3, d), "list of dispersal matrices")
  expect_error(kernels_by_step("a", 1, 1, 3, d), "list of dispersal matrices")
})

test_that("dispersal_time needs at least two sites", {
  expect_error(mc_dynamics(matrix(5, 1, 1), matrix(0.05), dispersal = 0, n_steps = 2,
                           dispersal_time = matrix(0.2, 2, 1)),
               "two sites")
})
