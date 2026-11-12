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
