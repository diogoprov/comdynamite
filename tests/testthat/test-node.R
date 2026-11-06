simulated <- function(lambda, n_sites = 50, n_species = 14) {
  tree <- ape::rcoal(n_species)
  traits <- species_traits(tree, lambda = lambda, optima_range = c(-1, 101),
                           rmax_range = c(0, 30))
  env <- runif(n_sites, 0, 100)
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  keep <- rowSums(comm) > 0
  list(comm = comm[keep, ], tree = tree, env = env[keep])
}

test_that("contrasts halve at each node, as worked out by hand", {
  tree <- ape::read.tree(text = "((a,b),(c,(d,e)));")
  x <- node_contrasts(tree)
  # nodes: root, (a,b) is left out, (c,(d,e)), (d,e) is left out
  expect_equal(ncol(x), 2)
  expect_equal(x[, 1], c(a = -0.5, b = -0.5, c = 0.5, d = 0.25, e = 0.25))
  expect_equal(x[, 2], c(a = 0, b = 0, c = -1, d = 0.5, e = 0.5))
  expect_error(node_contrasts(ape::read.tree(text = "(a,b,c);")), "dichotomous")
})

test_that("statistics equal the original function run in Octave", {
  # expected values from MetacommunityPhylogenetics in Octave 11.3, see
  # data-raw/validation_node.R
  set.seed(11)
  n <- 40
  k <- 12
  inc <- matrix(rbinom(n * k, 1, 0.45), n)
  inc[rowSums(inc) == 0, 1] <- 1
  env <- cbind(runif(n, 10, 30), rnorm(n, 5, 2))
  nodes <- cbind(c(-.5, -.5, .5, .25, .25, rep(0, 7)),
                 c(rep(0, 6), -.5, -.5, .5, .25, .125, .125),
                 c(-.25, -.25, -.25, -.125, -.125, .5, .25, .25, .125, .125, .0625, .0625))
  expected <- rbind(c(0.00783339, -0.05229792, 0.00799555),
                    c(0.03382408, -0.02138597, 0.03483741),
                    c(0.07671692, 0.02680973, 0.07970519))
  fits <- lapply(1:3, function(i) {
    m <- which(nodes[, i] != 0)
    node_fit(inc[, m, drop = FALSE], nodes[m, i], env)
  })
  mine <- t(vapply(fits, function(f) c(f$r2, f$r2_adjusted, f$statistic), numeric(3)))
  expect_equal(mine, expected, tolerance = 1e-6)
})

test_that("partial R2 with one predictor is the R2", {
  set.seed(1)
  d <- simulated(1)
  x <- node_contrasts(d$tree)
  out <- node_test(d$comm, x, d$env, n_perm = 3)
  expect_equal(unname(out$partial_r2[, 1]), out$nodes$r2)
})

test_that("the test detects clades that differ in their response", {
  set.seed(2)
  d <- simulated(1, n_sites = 80)
  out <- node_test(d$comm, node_contrasts(d$tree), d$env, n_perm = 99)
  expect_s3_class(out, "mc_node_test")
  expect_equal(nrow(out$nodes), ncol(node_contrasts(d$tree)))
  expect_lt(min(out$nodes$p_site, na.rm = TRUE), 0.05)
  expect_equal(out$nodes$p_max, pmax(out$nodes$p_site, out$nodes$p_taxa))
  expect_output(print(out), "tested")
})

test_that("species permutations draw from all species", {
  # the first species carry no signal: the original would shuffle only them
  set.seed(3)
  comm <- matrix(rbinom(60 * 12, 1, 0.5), 60, dimnames = list(NULL, paste0("s", 1:12)))
  env <- seq(0, 1, length.out = 60)
  comm[, 11] <- as.integer(env > 0.5)
  comm[, 12] <- as.integer(env <= 0.5)
  nodes <- matrix(0, 12, 1, dimnames = list(colnames(comm), "n"))
  nodes[11:12, 1] <- c(-1, 1)
  out <- node_test(comm, nodes, env, n_perm = 99)
  expect_lte(out$nodes$p_max, 0.02)
})

test_that("several predictors give one partial value each", {
  set.seed(4)
  d <- simulated(1)
  env <- cbind(temperature = d$env, noise = runif(nrow(d$comm)))
  out <- node_test(d$comm, node_contrasts(d$tree), env, n_perm = 5)
  expect_identical(colnames(out$partial_r2), colnames(env))
  expect_equal(dim(out$p_partial), dim(out$partial_r2))
  unnamed <- node_test(d$comm, node_contrasts(d$tree), unname(env), n_perm = 2)
  expect_identical(colnames(unnamed$partial_r2), c("predictor1", "predictor2"))
})

test_that("nodes that cannot be fitted give NA", {
  set.seed(5)
  comm <- matrix(rbinom(30 * 4, 1, 0.6), 30, dimnames = list(NULL, paste0("s", 1:4)))
  comm[, 3:4] <- 0
  comm[1:2, 3] <- 1
  comm[1:2, 4] <- c(1, 0)
  env <- cbind(c(1, 1, runif(28)), c(rep(1, 2), rep(2, 28)))
  nodes <- cbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  rownames(nodes) <- colnames(comm)
  out <- node_test(comm, nodes, env, n_perm = 2)
  expect_true(is.na(out$nodes$p_max[2]))
  expect_false(is.na(out$nodes$r2[1]))
  # one species only: no contrast within the node
  single <- cbind(c(1, 1, 0, 0))
  rownames(single) <- colnames(comm)
  expect_true(is.na(node_test(comm, single, runif(30), n_perm = 2)$nodes$r2))
  # collinear predictors
  expect_true(is.na(node_fit(comm[, 1:2], c(1, -1), cbind(1:30, 2 * (1:30)))$statistic))
  expect_null(weighted_fit(cbind(1:5, 1:5), 1:5, rep(1, 5)))
})

test_that("invalid data are rejected", {
  set.seed(6)
  d <- simulated(1)
  x <- node_contrasts(d$tree)
  expect_error(node_test(d$comm[, -1], x, d$env, n_perm = 5))
  expect_error(node_test(d$comm, x, d$env[-1], n_perm = 5))
})
