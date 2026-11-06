simulated <- function(lambda, n_sites = 60, n_species = 25) {
  tree <- ape::rcoal(n_species)
  traits <- species_traits(tree, lambda = lambda, optima_range = c(-1, 101),
                           rmax_range = c(0, 30))
  env <- runif(n_sites, 0, 100)
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  keep <- rowSums(comm) > 0
  list(comm = comm[keep, ], tree = tree, env = env[keep])
}

test_that("weighted standardisation gives weighted mean 0 and variance 1", {
  set.seed(1)
  w <- runif(30)
  z <- weighted_standardise(matrix(rnorm(60), 30), w)
  expect_equal(colSums(w * z) / sum(w), c(0, 0))
  expect_equal(colSums(w * z^2) / sum(w), c(1, 1))
})

test_that("eigenvectors are orthonormal, one fewer than species at most", {
  set.seed(2)
  tree <- ape::rcoal(12)
  v <- phylo_eigenvectors(tree)
  expect_equal(dim(v), c(12, 11))
  expect_identical(rownames(v), tree$tip.label)
  expect_equal(crossprod(v), diag(11), ignore_attr = TRUE)
  expect_equal(ncol(phylo_eigenvectors(tree, squared = TRUE)), 11)
})

test_that("with one predictor, R2 is the squared fourth-corner correlation", {
  set.seed(3)
  d <- simulated(1)
  comm <- d$comm[, d$tree$tip.label]
  vectors <- phylo_eigenvectors(d$tree)
  r <- rowSums(comm)
  composition <- (comm %*% weighted_standardise(vectors, colSums(comm))) / r
  e <- drop(weighted_standardise(as.matrix(d$env), r))
  fourth_corner <- colSums(r * e * composition) / sum(r)
  n <- nrow(comm)
  expect_equal(pcc_r2(comm, vectors, as.matrix(d$env)),
               1 - (n - 1) / (n - 2) * (1 - fourth_corner^2))
})

test_that("the test detects composition linked to phylogeny and environment", {
  set.seed(4)
  d <- simulated(1)
  out <- pcc_test(d$comm, d$tree, d$env, n_perm = 99)
  expect_s3_class(out, "mc_pcc")
  expect_lt(out$p_site, 0.05)
  expect_equal(out$p_max, max(out$p_site, out$p_taxa))
  expect_gt(out$r2, 0)
  expect_length(out$r2_adjusted, 24)
  expect_output(print(out), "p \\(max\\)")
})

test_that("an unrelated predictor is not rejected by the site permutations", {
  set.seed(5)
  d <- simulated(1)
  p <- replicate(30, pcc_test(d$comm, d$tree, runif(nrow(d$comm)), n_perm = 19)$p_site)
  expect_lt(mean(p <= 0.05), 0.25)
})

test_that("species order does not matter and several predictors are accepted", {
  set.seed(6)
  d <- simulated(1)
  vectors <- phylo_eigenvectors(d$tree)
  set.seed(7)
  a <- pcc_test(d$comm, d$tree, d$env, n_perm = 5)
  set.seed(7)
  b <- pcc_test(d$comm[, rev(colnames(d$comm))], d$tree, d$env, n_perm = 5)
  expect_equal(a$statistic, b$statistic)
  two <- pcc_test(d$comm, d$tree, cbind(d$env, runif(nrow(d$comm))), n_perm = 5)
  expect_length(two$r2_adjusted, ncol(vectors))
})

test_that("site permutations can be supplied", {
  set.seed(8)
  d <- simulated(1)
  n <- nrow(d$comm)
  perms <- t(replicate(9, sample(n)))
  out <- pcc_test(d$comm, d$tree, d$env, n_perm = 5, site_permutations = perms)
  expect_true(out$p_site %in% (1:10 / 10))
  expect_error(pcc_test(d$comm, d$tree, d$env, n_perm = 5, site_permutations = perms[, -1]))
})

test_that("invalid data are rejected", {
  set.seed(9)
  d <- simulated(1)
  empty <- d$comm
  empty[1, ] <- 0
  expect_error(pcc_test(empty, d$tree, d$env, n_perm = 5), "at least one individual")
  expect_error(pcc_test(d$comm[, -1], d$tree, d$env, n_perm = 5))
  expect_error(pcc_test(d$comm, d$tree, d$env[-1], n_perm = 5))
})
