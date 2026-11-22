# Quick checks of sensitivity to noise (G5.9a) and to the seed (G5.9b) of the permutation
# tests. G5.9a: a relative jitter of 1e-10 on the predictors and on the branch lengths, with
# the same seed, leaves the statistics equal up to 1e-6 (relative) and the p-values equal up
# to one step of the permutation distribution, 1 / (n_perm + 1), since a statistic tied with
# the observed one can fall on either side. G5.9b: with different seeds, each p-value of a
# permutation test with n_perm permutations has Monte Carlo standard error at most
# sqrt(p (1 - p) / n_perm); p-values must lie within 4 such errors (plus one step) of their
# mean, a bound exceeded with probability below 1e-4 per p-value under a normal approximation.

noise_data <- function(n_sites = 30, n_species = 12) {
  set.seed(31)
  tree <- ape::rcoal(n_species)
  tree$tip.label <- paste0("sp", seq_len(n_species))
  traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101), rmax_range = c(0, 30))
  env <- stats::runif(n_sites, 0, 100)
  comm <- simulate_static(traits, env + stats::rnorm(n_sites, 0, 30), incidence = TRUE)$abundance
  keep <- rowSums(comm) > 0
  comm <- comm[keep, colSums(comm[keep, ]) > 0]
  list(comm = comm, tree = ape::keep.tip(tree, colnames(comm)), env = env[keep])
}

jitter_data <- function(d) {
  set.seed(32)
  d$env <- d$env * (1 + 1e-10 * stats::rnorm(length(d$env)))
  d$tree$edge.length <- d$tree$edge.length * (1 + 1e-10 * stats::rnorm(length(d$tree$edge.length)))
  d
}

expect_seed_stable <- function(p, n_perm) {
  p <- as.matrix(p)
  centre <- rowMeans(p)
  bound <- 4 * sqrt(pmax(centre * (1 - centre), 1 / n_perm) / n_perm) + 1 / (n_perm + 1)
  expect_true(all(abs(p - centre) <= bound))
}

expect_jitter_invariant <- function(run, d, n_perm) {
  jittered <- jitter_data(d)
  set.seed(5)
  a <- run(d)
  set.seed(5)
  b <- run(jittered)
  expect_equal(b$statistic, a$statistic, tolerance = 1e-6)
  expect_true(all(abs(b$p - a$p) <= 1 / (n_perm + 1) + 1e-12))
}

test_that("pcc_test(), phylo_beta_test() and pcps_test() ignore trivial noise", {
  # G5.9a
  d <- noise_data()
  n_perm <- 99
  pcc <- function(d) {
    out <- pcc_test(d$comm, d$tree, d$env, n_perm = n_perm)
    list(statistic = out$statistic, p = c(out$p_site, out$p_species))
  }
  beta <- function(d) {
    out <- phylo_beta_test(d$comm, d$tree, d$env, n_perm = n_perm)
    list(statistic = out$statistic, p = c(out$p_site, out$p_species))
  }
  pcps <- function(d) {
    out <- pcps_test(d$comm, d$tree, d$env, axes = 1, n_perm = n_perm)$tests
    list(statistic = out$statistic, p = c(out$p_site, out$p_species))
  }
  expect_jitter_invariant(pcc, d, n_perm)
  expect_jitter_invariant(beta, d, n_perm)
  expect_jitter_invariant(pcps, d, n_perm)
})

test_that("p-values of the whole-tree tests vary across seeds within Monte Carlo error", {
  # G5.9b: 4 seeds, 199 permutations
  d <- noise_data()
  n_perm <- 199
  p <- sapply(1:4, function(seed) {
    set.seed(seed)
    pcc <- pcc_test(d$comm, d$tree, d$env, n_perm = n_perm)
    beta <- phylo_beta_test(d$comm, d$tree, d$env, n_perm = n_perm)
    pcps <- pcps_test(d$comm, d$tree, d$env, axes = 1, n_perm = n_perm)$tests
    c(pcc$p_site, pcc$p_species, beta$p_site, beta$p_species, pcps$p_site, pcps$p_species)
  })
  expect_seed_stable(p, n_perm)
})

test_that("node_test() and phylo_vector() ignore trivial noise and the seed", {
  # G5.9a and G5.9b
  d <- noise_data()
  n_perm <- 99
  node <- function(d) {
    out <- node_test(d$comm, node_contrasts(d$tree)[, 1:3], d$env, n_perm = n_perm)$nodes
    list(statistic = out$statistic, p = c(out$p_site, out$p_species))
  }
  vector <- function(d) {
    out <- phylo_vector(d$comm, d$tree, d$env, n_perm = n_perm, n_sim = 10)
    list(statistic = out$eigenvector_correlation, p = c(out$p_site, out$p_species))
  }
  expect_jitter_invariant(node, d, n_perm)
  expect_jitter_invariant(vector, d, n_perm)
  n_perm <- 199
  p <- sapply(1:4, function(seed) {
    set.seed(seed)
    out <- node_test(d$comm, node_contrasts(d$tree)[, 1:3], d$env, n_perm = n_perm)$nodes
    vec <- phylo_vector(d$comm, d$tree, d$env, n_perm = n_perm, n_sim = 10)
    c(out$p_site, out$p_species, vec$p_site, vec$p_species)
  })
  expect_seed_stable(p, n_perm)
})

test_that("trait_partition() ignores trivial noise and the seed", {
  # G5.9a and G5.9b, with a transect of 30 sites and polynomial terms of the position as space
  d <- noise_data()
  trait <- stats::setNames(seq_len(ncol(d$comm)) * 7, colnames(d$comm))
  n <- nrow(d$comm)
  space <- stats::poly(seq_len(n), 3)
  rownames(space) <- rownames(d$comm)
  env <- cbind(e = d$env)
  rownames(env) <- rownames(d$comm)
  n_perm <- 99
  run <- function(env, space, trait) {
    f <- trait_partition(d$comm, trait, env, space, n_perm = n_perm)$fractions
    list(statistic = f$statistic, p = c(f$p_site, f$p_species))
  }
  set.seed(5)
  a <- run(env, space, trait)
  set.seed(32)
  jittered <- function(x) x * (1 + 1e-10 * stats::rnorm(length(x)))
  env_j <- env
  env_j[] <- jittered(env)
  space_j <- space
  space_j[] <- jittered(space)
  trait_j <- stats::setNames(jittered(trait), names(trait))
  set.seed(5)
  b <- run(env_j, space_j, trait_j)
  expect_equal(b$statistic, a$statistic, tolerance = 1e-6)
  expect_true(all(abs(b$p - a$p) <= 1 / (n_perm + 1) + 1e-12, na.rm = TRUE))
  n_perm <- 199
  p <- sapply(1:4, function(seed) {
    set.seed(seed)
    run(env, space, trait)$p
  })
  p <- p[stats::complete.cases(p), , drop = FALSE]
  expect_seed_stable(p, n_perm)
})

test_that("torus_permutations() keep distances on the torus whatever the seed", {
  # G5.9b for the restricted null itself: a random subset of moves, drawn with any seed, keeps
  # every distance among sites on the torus, so it changes the draw and not the null model.
  side <- 6
  xy <- as.matrix(expand.grid(seq_len(side), seq_len(side)))
  torus_distance <- function(i) {
    dx <- abs(outer(xy[i, 1], xy[i, 1], "-"))
    dy <- abs(outer(xy[i, 2], xy[i, 2], "-"))
    sqrt(pmin(dx, side - dx)^2 + pmin(dy, side - dy)^2)
  }
  reference <- torus_distance(seq_len(nrow(xy)))
  for (seed in 1:3) {
    set.seed(seed)
    moves <- torus_permutations(side, n_perm = 20)
    for (r in seq_len(nrow(moves))) {
      expect_equal(torus_distance(moves[r, ]), reference)
    }
  }
})
