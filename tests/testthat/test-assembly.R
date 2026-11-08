test_that("GAUSE removes, from the closest pair, the species with smaller summed distance", {
  # a and b are the closest pair; b is closer to the others (c, d), so it goes first
  x <- c(a = 0, b = 0.5, c = 3, d = 4)
  out <- assemble_community(x, size = 3, rule = "limiting_similarity")
  expect_equal(out$removed, "b")
  expect_setequal(out$community, c("a", "c", "d"))
  # next removal: c and d are now the closest pair; c has the smaller sum
  out2 <- assemble_community(x, size = 2, rule = "limiting_similarity")
  expect_equal(out2$removed, c("b", "c"))
})

test_that("limiting similarity spaces species more evenly than random draws", {
  set.seed(1)
  x <- stats::setNames(stats::runif(40), paste0("s", 1:40))
  nn <- function(v) min(diff(sort(v)))
  ls <- replicate(50, nn(x[assemble_community(x, 10, "limiting_similarity")$community]))
  nt <- replicate(50, nn(x[assemble_community(x, 10, "neutral")$community]))
  expect_gt(mean(ls), mean(nt))
})

test_that("min_distance stops removal at a limit to similarity", {
  x <- c(a = 0, b = 0.1, c = 1, d = 2.05, e = 3)
  out <- assemble_community(x, size = 1, rule = "limiting_similarity", min_distance = 0.5)
  expect_equal(out$removed, c("b"))
  expect_true(min(dist(x[out$community])) >= 0.5)
})

test_that("filtering keeps the species nearest the optimum", {
  x <- cbind(t1 = c(a = -2, b = -0.1, c = 0.2, d = 3), t2 = 0)
  expect_setequal(assemble_community(x, 2, "filtering")$community, c("b", "c"))
  expect_setequal(assemble_community(x, 2, "filtering", optimum = "derived")$community,
                  c("c", "d"))
  expect_setequal(assemble_community(x, 1, "filtering", optimum = c(-2, 0))$community, "a")
})

test_that("multivariate traits and output structure", {
  set.seed(2)
  tr <- ape::rcoal(30)
  x <- traits_ou(tr, n_traits = 3)
  expect_equal(dim(x), c(30, 3))
  for (rule in c("limiting_similarity", "filtering", "neutral")) {
    out <- assemble_community(x, 12, rule)
    expect_length(out$community, 12)
    expect_length(out$removed, 18)
    expect_setequal(c(out$community, out$removed), tr$tip.label)
  }
  expect_error(assemble_community(x, 31))
})

test_that("traits without row names get default species names", {
  set.seed(1)
  x <- matrix(rnorm(12), 6, 2)
  expect_warning(out <- assemble_community(x, size = 3, rule = "neutral"), "no names")
  expect_length(out$community, 3)
  expect_true(all(out$community %in% paste0("sp", 1:6)))
})

test_that("distance changes trait space as documented", {
  set.seed(10)
  x <- cbind(a = stats::rnorm(30), b = stats::rnorm(30, sd = 10))
  rownames(x) <- paste0("s", 1:30)
  # standardized: the same as Euclidean on traits divided by their standard deviations
  z <- sweep(x, 2, apply(x, 2, stats::sd), "/")
  set.seed(1)
  a <- assemble_community(x, 10, distance = "standardized")
  set.seed(1)
  b <- assemble_community(z, 10)
  expect_identical(a, b)
  # mahalanobis: the same as Euclidean on whitened traits
  w <- x %*% solve(chol(stats::cov(x)))
  set.seed(2)
  a <- assemble_community(x, 10, rule = "filtering", optimum = c(0, 0),
                          distance = "mahalanobis")
  set.seed(2)
  b <- assemble_community(w, 10, rule = "filtering", optimum = c(0, 0))
  expect_identical(a, b)
  # a given covariance replaces the one of the pool
  set.seed(3)
  a <- assemble_community(x, 10, distance = "mahalanobis",
                          trait_cov = diag(apply(x, 2, stats::var)))
  set.seed(3)
  b <- assemble_community(x, 10, distance = "standardized")
  expect_equal(a, b)
})

test_that("Mahalanobis filtering keeps species inside the niche ellipse", {
  # two traits along one direction: Euclidean distance favours the short axis, Mahalanobis not
  x <- rbind(along = c(2, 2), across = c(1, -1), rbind(matrix(c(1, 1, -1, -1, 3, 3, -3, -3,
                                                                0.2, -0.2, -0.2, 0.2),
                                                              ncol = 2, byrow = TRUE)))
  rownames(x) <- c("along", "across", paste0("s", 1:6))
  s <- matrix(c(1, 0.95, 0.95, 1), 2)
  euc <- assemble_community(x, 7, rule = "filtering", optimum = c(0, 0))
  mah <- assemble_community(x, 7, rule = "filtering", optimum = c(0, 0),
                            distance = "mahalanobis", trait_cov = s)
  expect_true("across" %in% euc$community)
  expect_false("across" %in% mah$community)
})

test_that("distance checks its inputs", {
  x <- matrix(stats::rnorm(6), 3, dimnames = list(c("a", "b", "c"), NULL))
  expect_error(assemble_community(x[1:2, ], 1, distance = "mahalanobis"),
               "more species than traits")
  expect_error(assemble_community(x, 2, trait_cov = diag(2)), "only with")
  expect_error(assemble_community(x, 2, distance = "mahalanobis", trait_cov = diag(3)),
               "one row and column per trait")
  expect_error(assemble_community(x, 2, distance = "mahalanobis", trait_cov = matrix(1, 2, 2)),
               "singular")
  expect_error(assemble_community(cbind(x, 1), 2, distance = "standardized"), "vary")
})
