# Cheap checks of sensitivity to noise (G5.9a) and to the seed (G5.9b) for the evolutionary
# and DAMOCLES functions; they run on CRAN and CI.

noise_tree <- function(n) {
  if (requireNamespace("treeshapesim", quietly = TRUE)) {
    tree <- treeshapesim::sim_tree_shape(n)
  } else {
    tree <- ape::rcoal(n)
  }
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  tree
}

jitter_rel <- function(x, size = 1e-10) x * (1 + stats::runif(length(x), -size, size))

test_that("DAMOCLES likelihoods and fits are insensitive to trivial noise in branch lengths", {
  # G5.9a: a relative jitter of 1e-10 changes log-likelihoods by about 1e-10 times their size
  set.seed(1)
  tree <- noise_tree(40)
  comm <- simulate_damocles(tree, 1, 1, n = 5)
  noisy <- tree
  noisy$edge.length <- jitter_rel(tree$edge.length)
  expect_equal(damocles_loglik(comm, noisy, 1.3, 0.7), damocles_loglik(comm, tree, 1.3, 0.7),
               tolerance = 1e-8)
  a <- damocles_fit(comm, tree, shared = TRUE)
  b <- damocles_fit(comm, noisy, shared = TRUE)
  expect_equal(attr(b, "loglik_total"), attr(a, "loglik_total"), tolerance = 1e-6)
  expect_equal(b$gamma, a$gamma, tolerance = 1e-3)
})

test_that("damocles_test() p-values differ between seeds only by Monte Carlo error", {
  # G5.9b: with 999 null communities, a p-value has a binomial standard error of at most
  # sqrt(p (1 - p) / 999); two seeds must agree within 4 standard errors of their difference.
  set.seed(2)
  tree <- noise_tree(40)
  comm <- simulate_damocles(tree, 1, 1, n = 3)
  set.seed(10)
  a <- damocles_test(comm, tree, n_sim = 999, gamma = 1, mu = 1)
  set.seed(20)
  b <- damocles_test(comm, tree, n_sim = 999, gamma = 1, mu = 1)
  for (col in c("p_mpd_damocles", "p_mntd_damocles", "p_mpd_random", "p_mntd_random")) {
    p <- (a[[col]] + b[[col]]) / 2
    bound <- 4 * sqrt(2 * p * (1 - p) / 999) + 2 / 1000
    expect_true(all(abs(a[[col]] - b[[col]]) <= bound, na.rm = TRUE))
  }
})

test_that("dr_rates() is insensitive to trivial noise in branch lengths", {
  # G5.9a: DR is a smooth function of branch lengths
  set.seed(3)
  tree <- noise_tree(50)
  noisy <- tree
  noisy$edge.length <- jitter_rel(tree$edge.length)
  expect_equal(dr_rates(noisy), dr_rates(tree), tolerance = 1e-8)
  expect_equal(dr_rates(noisy, per = "edge"), dr_rates(tree, per = "edge"), tolerance = 1e-8)
})

test_that("drift_test() is insensitive to trivial noise in the data", {
  # G5.9a: slope and p-value change by far less than their sampling error
  set.seed(4)
  G <- diag(c(4, 2, 1, 0.5, 0.25))
  means <- matrix(stats::rnorm(20 * 5), 20) %*% chol(0.5 * G)
  ind <- means[rep(1:20, each = 10), ] + matrix(stats::rnorm(200 * 5), 200) %*% chol(G)
  a <- drift_test(ind, group = rep(1:20, each = 10))
  b <- drift_test(jitter_rel(ind), group = rep(1:20, each = 10))
  expect_equal(b$slope, a$slope, tolerance = 1e-6)
  expect_equal(b$p_value, a$p_value, tolerance = 1e-6)
})

test_that("persistence() is insensitive to trivial noise in abundances", {
  # G5.9a: presence depends only on abundance above zero, which a relative jitter keeps
  set.seed(5)
  x <- array(stats::rpois(20 * 3 * 10, 1), c(20, 3, 10))
  noisy <- array(jitter_rel(x), dim(x))
  expect_identical(persistence(noisy), persistence(x))
})
