tr4 <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")

# covariance of the traits of one tip over many simulations
tip_cov <- function(n, fun, tip = "a") stats::cov(t(replicate(n, fun()[tip, ])))

test_that("rate_matrix builds each single-matrix structure", {
  m <- rate_matrix(c(1, 4), correlation = 0.5)
  expect_equal(unname(m), matrix(c(1, 1, 1, 4), 2), ignore_attr = TRUE)
  expect_equal(dimnames(m), list(c("trait1", "trait2"), c("trait1", "trait2")))
  expect_identical(attr(m, "constraint"), "default")
  expect_equal(unname(rate_matrix(c(1, 2, 3), constraint = "diagonal")), diag(c(1, 2, 3)),
               ignore_attr = TRUE)
  expect_equal(unname(rate_matrix(2, constraint = "equaldiagonal", n_traits = 3)), diag(2, 3),
               ignore_attr = TRUE)
  e <- rate_matrix(2, correlation = 0.25, constraint = "equal", n_traits = 3)
  expect_equal(unname(diag(e)), rep(2, 3))
  expect_equal(e[1, 3], 0.5)
  rho <- matrix(c(1, 0.2, -0.3, 0.2, 1, 0, -0.3, 0, 1), 3)
  expect_equal(unname(stats::cov2cor(rate_matrix(c(1, 2, 3), correlation = rho))), rho,
               ignore_attr = TRUE)
})

test_that("rate_matrix builds each multi-regime structure", {
  p <- rate_matrix(c(1, 2), correlation = 0.5, constraint = "proportional",
                   scale = c(a = 1, b = 3))
  expect_named(p, c("a", "b"))
  expect_equal(p$b, 3 * p$a)
  s <- rate_matrix(c(1, 2), correlation = 0.5, constraint = "shared",
                   eigenvalues = list(a = c(3, 1), b = c(1, 2)))
  va <- eigen(s$a, symmetric = TRUE)
  vb <- eigen(s$b, symmetric = TRUE)
  expect_equal(va$values, c(3, 1))
  expect_equal(sort(vb$values), c(1, 2))
  expect_equal(abs(crossprod(va$vectors, vb$vectors)), diag(2)[, 2:1])
  v <- rate_matrix(c(1, 2), correlation = list(a = 0.5, b = -0.5), constraint = "variance")
  expect_equal(diag(v$a), diag(v$b))
  expect_equal(v$b[1, 2], -v$a[1, 2])
  k <- rate_matrix(variance = list(a = c(1, 2), b = c(4, 1)), correlation = 0.3,
                   constraint = "correlation")
  expect_equal(stats::cov2cor(k$a), stats::cov2cor(k$b))
  expect_equal(unname(diag(k$b)), c(4, 1))
})

test_that("rate_matrix checks its inputs", {
  expect_error(rate_matrix(1, correlation = 0.5, constraint = "equal"), "n_traits")
  expect_error(rate_matrix(c(1, 2), correlation = 0.5, constraint = "diagonal"), "no correlations")
  expect_error(rate_matrix(c(1, 2), constraint = "equal"), "one variance")
  expect_error(rate_matrix(c(1, 2), correlation = 1.5), "\\[-1, 1\\]")
  expect_error(rate_matrix(1, correlation = -0.9, n_traits = 3), "positive semidefinite")
  expect_error(rate_matrix(c(-1, 1)), "zero or more")
  expect_error(rate_matrix(c(1, 2), constraint = "proportional"), "scale")
  expect_error(rate_matrix(c(1, 2), constraint = "shared", eigenvalues = list(a = 1)),
               "2 values")
  expect_error(rate_matrix(c(1, 2), correlation = 0.5, constraint = "variance"),
               "list named by regime")
  expect_error(rate_matrix(list(a = 1, b = 2), constraint = "default", n_traits = 2),
               "only with constraint")
  expect_error(rate_matrix(1, correlation = list(a = 0.5), n_traits = 2),
               "only with constraint \"variance\"")
  expect_error(rate_matrix(1, n_traits = 0), "one positive number")
  expect_error(rate_matrix(1, correlation = 0.5, constraint = "equaldiagonal", n_traits = 2),
               "no correlations")
  expect_error(rate_matrix(c(1, 2), constraint = "equaldiagonal"), "one variance")
  expect_error(rate_matrix(c(1, 2), correlation = matrix(c(1, 0.5, 0.2, 1), 2)),
               "symmetric 2 x 2")
  expect_error(rate_matrix(1, correlation = diag(3), n_traits = 2), "symmetric 2 x 2")
})

test_that("a rate matrix gives correlated traits with covariance R times the path length", {
  set.seed(1)
  r <- rate_matrix(c(1, 0.5), correlation = 0.8)
  v <- tip_cov(20000, function() traits_ou(tr4, sigma2 = r))
  expect_equal(unname(v), unname(r) * 2, tolerance = 0.05, ignore_attr = TRUE)
  x <- traits_ou(tr4, sigma2 = r)
  expect_equal(colnames(x), c("trait1", "trait2"))
  expect_equal(dim(attr(x, "nodes")), c(3, 2))
})

test_that("a diagonal rate matrix with equal rates reproduces scalar sigma2", {
  set.seed(5)
  a <- traits_ou(tr4, n_traits = 3, alpha = 0.5, sigma2 = 2)
  set.seed(5)
  b <- traits_ou(tr4, alpha = 0.5, sigma2 = rate_matrix(2, constraint = "equaldiagonal",
                                                         n_traits = 3))
  expect_equal(a, b)
})

test_that("a rate matrix under OU has the stationary-type covariance of Hansen (1997)", {
  set.seed(2)
  r <- rate_matrix(c(1, 2), correlation = -0.5)
  a <- 1.5
  v <- tip_cov(20000, function() traits_ou(tr4, alpha = a, sigma2 = r))
  expect_equal(unname(v), unname(r) * (1 - exp(-2 * a * 2)) / (2 * a), tolerance = 0.05,
               ignore_attr = TRUE)
})

test_that("one rate matrix per regime, with regimes by branch", {
  set.seed(3)
  reg <- paint_regimes(tr4, nodes = 7)  # clade (c, d), from node 7 down
  r <- rate_matrix(c(1, 1), correlation = list(`0` = 0.9, `1` = -0.9), constraint = "variance")
  sim <- function() traits_ou(tr4, sigma2 = r, regimes = reg)
  # tip c spent one unit of time in each regime, tip a two units in regime 0
  expect_equal(unname(tip_cov(20000, sim, "a")), unname(2 * r$`0`), tolerance = 0.05)
  expect_equal(unname(tip_cov(20000, sim, "c")), unname(r$`0` + r$`1`), tolerance = 0.1)
  # scalar rates per regime
  v <- stats::var(traits_ou(tr4, n_traits = 20000, sigma2 = c(`0` = 1, `1` = 3),
                            regimes = reg)["c", ])
  expect_equal(v, 4, tolerance = 0.05)
})

test_that("a mapped history gives regimes that change within branches", {
  h <- tr4
  h$maps <- lapply(h$edge.length, function(len) c(a = len))
  h$maps[[which(h$edge[, 2] == 3)]] <- c(a = 0.25, b = 0.75)  # tip c
  h$mapped.edge <- t(vapply(h$maps, function(m) c(a = sum(m[names(m) == "a"]),
                                                  b = sum(m[names(m) == "b"])), numeric(2)))
  class(h) <- c("simmap", "phylo")
  set.seed(4)
  x <- traits_ou(h, n_traits = 20000, sigma2 = c(a = 1, b = 5), regimes = h)
  expect_equal(stats::var(x["c", ]), 1.25 + 0.75 * 5, tolerance = 0.05)
  expect_equal(stats::var(x["d", ]), 2, tolerance = 0.05)
  # optima per state also change within the branch
  y <- traits_ou(h, n_traits = 5000, alpha = 4, sigma2 = 0.01, theta = c(a = 0, b = 10),
                 regimes = h)
  expect_equal(mean(y["c", ]), 10 * (1 - exp(-4 * 0.75)), tolerance = 0.02)
  expect_equal(mean(y["d", ]), 0, tolerance = 0.05)
})

test_that("traits_ou checks rate matrices and mapped regimes", {
  reg <- paint_regimes(tr4, nodes = 7)
  expect_error(traits_ou(tr4, sigma2 = matrix(c(1, 2, 0, 1), 2)), "symmetric")
  expect_error(traits_ou(tr4, sigma2 = matrix(c(1, 2, 2, 1), 2)), "semidefinite")
  expect_error(traits_ou(tr4, n_traits = 3, sigma2 = diag(2)), "n_traits")
  expect_error(traits_ou(tr4, sigma2 = list(`0` = diag(2), `1` = 1), regimes = reg), "not both")
  expect_error(traits_ou(tr4, sigma2 = list(`0` = diag(2), `1` = diag(3)), regimes = reg),
               "same size")
  expect_error(traits_ou(tr4, sigma2 = list(`0` = diag(2)), regimes = reg), "every regime")
  expect_error(traits_ou(tr4, regimes = tr4), "simmap")
  expect_error(traits_ou(tr4, sigma2 = -1), "zero or more")
  # a singular matrix (perfectly correlated traits) is allowed
  set.seed(6)
  x <- traits_ou(tr4, sigma2 = matrix(1, 2, 2))
  expect_equal(x[, 1], x[, 2])
})
