test_that("landscape_points gives Euclidean or great-circle distances and a kernel", {
  xy <- cbind(x = c(0, 3, 0), y = c(0, 4, 1))
  ls <- landscape_points(xy, env = c(1, 2, 3), disp_scale = 1)
  expect_s3_class(ls, "mc_landscape")
  expect_equal(ls$dist[1, 2], 5)
  expect_equal(unname(rowSums(ls$dispersal)), rep(1, 3))
  expect_equal(unname(diag(ls$dispersal)), rep(0, 3))
  expect_output(print(ls), "points")
  # one degree of latitude is about 111.2 km
  ll <- landscape_points(cbind(c(-72, -72), c(43, 44)), c(0, 1), 10, lonlat = TRUE)
  expect_equal(ll$dist[1, 2], 6371 * pi / 180, tolerance = 1e-6)
  expect_error(landscape_points(xy, env = 1:2, disp_scale = 1))
  expect_error(landscape_points(xy, env = 1:3, disp_scale = 0))
})

test_that("landscape_points drives simulations and the daphniid lakes", {
  skip_if_not_installed("treeshapesim")
  set.seed(1)
  ls <- landscape_points(daphniids$coords[1:20, ], env = seq(0, 100, length.out = 20),
                         disp_scale = 50, lonlat = TRUE)
  tree <- treeshapesim::sim_tree_shape(6)
  sim <- simulate_metacommunity(ls, species_traits(tree, 0.5), competition_matrix(tree, "none"),
                                dispersal = 0.05, n_steps = 20)
  expect_equal(dim(sim$abundance), c(20, 6))
})

test_that("landscape_matrix() keeps a given dispersal matrix, rescaled by rows", {
  w <- matrix(c(0, 2, 2, 1, 0, 3, 1, 1, 0), 3, byrow = TRUE)
  diag(w) <- 5
  ls <- landscape_matrix(w, env = c(1, 2, 3), coords = cbind(1:3, 0), dist = matrix(1, 3, 3))
  expect_s3_class(ls, "mc_landscape")
  expect_equal(diag(ls$dispersal), rep(0, 3))
  expect_equal(rowSums(ls$dispersal), rep(1, 3))
  expect_equal(ls$dispersal[2, ], c(1, 0, 3) / 4)
  expect_output(print(ls), "dispersal matrix")
  expect_error(landscape_matrix(diag(3), env = 1:3), "destination")
  expect_error(landscape_matrix(-w, env = 1:3))
  # a barrier in the matrix keeps populations apart
  block <- matrix(0, 6, 6)
  block[1:3, 1:3] <- block[4:6, 4:6] <- 1
  ls2 <- landscape_matrix(block, env = rep(50, 6))
  set.seed(3)
  tr <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sim <- simulate_metacommunity(ls2, tr, matrix(0.05, dimnames = list("sp1", "sp1")),
                                dispersal = 0.5, n_steps = 50, seed_steps = numeric(0),
                                initial = matrix(c(10, 0, 0, 0, 0, 0), 6, 1,
                                                 dimnames = list(NULL, "sp1")))
  expect_true(all(sim$abundance[4:6, ] == 0))
  expect_true(all(sim$abundance[2:3, ] > 0))
})
