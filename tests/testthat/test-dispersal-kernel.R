test_that("dispersal_kernel gives row-stochastic matrices of both shapes", {
  d <- unname(as.matrix(dist(cbind(0:5, 0))))
  for (shape in c("exponential", "gaussian")) {
    k <- dispersal_kernel(d, scale = 1, shape = shape)
    expect_equal(rowSums(k), rep(1, 6))
    expect_identical(diag(k), rep(0, 6))
  }
  e <- dispersal_kernel(d, 1)
  g <- dispersal_kernel(d, 1, "gaussian")
  expect_equal(e[1, 2] / e[1, 3], exp(1))
  expect_equal(g[1, 2] / g[1, 3], exp((4 - 1) / 2))
  # thinner tail: less dispersal to the farthest site
  expect_lt(g[1, 6], e[1, 6])
  expect_error(dispersal_kernel(d, 1, "cauchy"))
})

test_that("dispersal_kernel does not underflow for far sites", {
  d <- unname(as.matrix(dist(cbind(c(0, 1000, 2000), 0))))
  k <- dispersal_kernel(d, scale = 1, shape = "gaussian")
  expect_false(anyNA(k))
  expect_identical(k[1, 2], 1)
  expect_equal(k[2, ], c(0.5, 0, 0.5))
})

test_that("landscape constructors pass the kernel shape", {
  set.seed(1)
  g <- landscape_grid(side = 4, disp_scale = 1, disp_shape = "gaussian")
  expect_identical(g$dispersal, dispersal_kernel(g$dist, 1, "gaussian"))
  p <- landscape_points(cbind(0:3, 0), env = 1:4, disp_scale = 1, disp_shape = "gaussian")
  expect_identical(p$dispersal, dispersal_kernel(p$dist, 1, "gaussian"))
  set.seed(1)
  e <- landscape_grid(side = 4, disp_scale = 1)
  expect_identical(e$dispersal, dispersal_kernel(e$dist, 1))
})
