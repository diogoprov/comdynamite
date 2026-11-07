simulated <- function(lambda, n_sites = 60, n_species = 25) {
  tree <- ape::rcoal(n_species)
  traits <- species_traits(tree, lambda = lambda, optima_range = c(-1, 101),
                           rmax_range = c(0, 30))
  env <- runif(n_sites, 0, 100)
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  keep <- rowSums(comm) > 0
  list(comm = comm[keep, ], tree = tree, env = env[keep])
}

# two species 10 apart: t1 and t2 are sister tips of branches of length 5
pair_tree <- ape::read.tree(text = "(t1:5,t2:5);")

test_that("the metric gives the mean distance worked out by hand", {
  comm <- rbind(a = c(t1 = 1, t2 = 0), b = c(0, 1), c = c(1, 1), d = c(3, 1))
  beta <- as.matrix(phylo_beta(comm, pair_tree))
  expect_equal(beta["a", "b"], 10)
  expect_equal(beta["a", "c"], 5)
  expect_equal(beta["c", "b"], 5)
  # weights 3/4 and 1/4 against the single species t2
  expect_equal(beta["d", "b"], 7.5)
  expect_equal(unname(diag(beta)), rep(0, 4))
})

test_that("the metric equals picante::comdist", {
  skip_if_not_installed("picante")
  set.seed(1)
  d <- simulated(1)
  distances <- ape::cophenetic.phylo(d$tree)
  expect_equal(c(phylo_beta(d$comm, d$tree)), c(picante::comdist(d$comm, distances)))
  abundance <- d$comm * rpois(length(d$comm), 5) + d$comm
  expect_equal(c(phylo_beta(abundance, d$tree)),
               c(picante::comdist(abundance, distances, abundance.weighted = TRUE)))
})

test_that("the matrix is not always Euclidean", {
  # a and b lie halfway between c and d, so they should coincide, but they are 5 apart
  comm <- rbind(a = c(t1 = 1, t2 = 1), b = c(1, 1), c = c(1, 0), d = c(0, 1))
  centred <- gower_centre(as.matrix(phylo_beta(comm, pair_tree)))
  expect_lt(min(eigen(centred, symmetric = TRUE)$values), -1)
})

test_that("pseudo-F and R2 equal a distance-based redundancy analysis in vegan", {
  skip_if_not_installed("vegan")
  set.seed(2)
  d <- simulated(1)
  env <- data.frame(a = d$env, b = runif(nrow(d$comm)))
  beta <- phylo_beta(d$comm, d$tree)
  reference <- vegan::adonis2(beta ~ a + b, data = env, by = NULL, permutations = 0)
  n <- nrow(d$comm)
  out <- phylo_beta_test(d$comm, d$tree, env, n_perm = 2)
  expect_equal(out$statistic, reference$F[1])
  expect_equal(out$r2_adjusted, 1 - (1 - reference$R2[1]) * (n - 1) / (n - 3))
})

test_that("scores reproduce the predicted variation", {
  set.seed(3)
  d <- simulated(1)
  env <- cbind(d$env, runif(nrow(d$comm)))
  out <- phylo_beta_test(d$comm, d$tree, env, n_perm = 2)
  expect_lte(ncol(out$scores), 2)
  expect_equal(nrow(out$scores), nrow(d$comm))
  centred <- gower_centre(as.matrix(phylo_beta(d$comm, d$tree)))
  basis <- qr.Q(qr(sweep(env, 2, colMeans(env))))
  expect_equal(sum(out$scores^2), sum(basis * (centred %*% basis)))
})

test_that("the test detects composition linked to phylogeny and environment", {
  set.seed(4)
  d <- simulated(1)
  out <- phylo_beta_test(d$comm, d$tree, d$env, n_perm = 99)
  expect_s3_class(out, "mc_phylo_beta")
  expect_lt(out$p_site, 0.05)
  expect_lt(out$p_species, 0.05)
  expect_equal(out$p_max, max(out$p_site, out$p_species))
  expect_gt(out$r2_adjusted, 0)
  expect_output(print(out), "p \\(max\\)")
})

test_that("an unrelated predictor is not rejected by the site permutations", {
  set.seed(5)
  d <- simulated(1)
  p <- replicate(30, phylo_beta_test(d$comm, d$tree, runif(nrow(d$comm)), n_perm = 19)$p_site)
  expect_lt(mean(p <= 0.05), 0.25)
})

test_that("species order and the scale of predictors do not matter", {
  set.seed(6)
  d <- simulated(1)
  a <- phylo_beta_test(d$comm, d$tree, d$env, n_perm = 2)
  b <- phylo_beta_test(d$comm[, rev(colnames(d$comm))], d$tree, 3 * d$env + 7, n_perm = 2)
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$r2_adjusted, b$r2_adjusted)
})

test_that("site permutations can be supplied", {
  set.seed(7)
  d <- simulated(1)
  n <- nrow(d$comm)
  perms <- t(replicate(9, sample(n)))
  out <- phylo_beta_test(d$comm, d$tree, d$env, n_perm = 5, site_permutations = perms)
  expect_true(out$p_site %in% (1:10 / 10))
  expect_error(phylo_beta_test(d$comm, d$tree, d$env, n_perm = 5,
                               site_permutations = perms[, -1]))
})

test_that("invalid data are rejected", {
  set.seed(8)
  d <- simulated(1)
  empty <- d$comm
  empty[1, ] <- 0
  expect_error(phylo_beta(empty, d$tree), "at least one individual")
  expect_error(phylo_beta_test(d$comm, d$tree, d$env[-1], n_perm = 5))
  expect_error(phylo_beta_test(d$comm, d$tree, cbind(d$env, 2 * d$env), n_perm = 5),
               "collinear")
})
