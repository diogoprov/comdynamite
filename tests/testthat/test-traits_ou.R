tr4 <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")

test_that("alpha = 0 gives Brownian motion with the phylogenetic covariance", {
  set.seed(1)
  x <- traits_ou(tr4, n_traits = 20000, alpha = 0, sigma2 = 2)
  v <- stats::cov(t(x))
  expect_equal(unname(v), unname(2 * ape::vcv.phylo(tr4)[rownames(x), rownames(x)]),
               tolerance = 0.05)
  expect_equal(unname(rowMeans(x)), rep(0, 4), tolerance = 0.05)
})

test_that("single-optimum OU has the expected tip mean and variance", {
  set.seed(2)
  a <- 1.5; s2 <- 1; th <- 3; T <- 2
  x <- traits_ou(tr4, n_traits = 20000, alpha = a, sigma2 = s2, theta = th, root = 0)
  expect_equal(mean(x["a", ]), th * (1 - exp(-a * T)), tolerance = 0.02)
  expect_equal(stats::var(x["a", ]), s2 / (2 * a) * (1 - exp(-2 * a * T)), tolerance = 0.05)
  # sister tips share their path from the root to node 6
  expect_equal(stats::cov(x["a", ], x["b", ]),
               s2 / (2 * a) * (1 - exp(-2 * a * 1)) * exp(-2 * a * 1), tolerance = 0.1)
})

test_that("multi-peak OU pulls each clade towards its own optimum", {
  set.seed(3)
  reg <- paint_regimes(tr4, nodes = 7)  # clade (c, d)
  expect_equal(sum(reg == "1"), 2)
  x <- traits_ou(tr4, n_traits = 5000, alpha = 5, sigma2 = 0.1,
                 theta = c(`0` = 0, `1` = 10), regimes = reg)
  expect_equal(unname(rowMeans(x)[c("a", "b")]), c(0, 0), tolerance = 0.05)
  expect_equal(unname(rowMeans(x)[c("c", "d")]), c(10, 10) * (1 - exp(-5)), tolerance = 0.05)
  stem <- paint_regimes(tr4, nodes = 7, stem = TRUE)
  expect_equal(sum(stem == "1"), 3)
})

test_that("traits_ou checks its inputs", {
  expect_error(traits_ou(tr4, theta = c(0, 1)))
  expect_error(traits_ou(tr4, theta = c(x = 0, y = 1), regimes = paint_regimes(tr4, 7)))
  expect_error(traits_ou(tr4, alpha = -1))
  expect_error(paint_regimes(tr4, nodes = 2))
  x <- traits_ou(ape::rcoal(10), n_traits = 2)
  expect_equal(dim(x), c(10, 2))
})

test_that("one unnamed optimum applies to all regimes", {
  set.seed(4)
  reg <- paint_regimes(tr4, nodes = 7)
  x <- traits_ou(tr4, n_traits = 5000, alpha = 5, sigma2 = 0.1, theta = 2,
                 regimes = reg, root = 2)
  expect_equal(unname(rowMeans(x)), rep(2, 4), tolerance = 0.05)
})
