# Type I error and power of the tests of the fractions of trait_partition(): `a` and `c` with
# the fraction itself as statistic (the original code, and this package up to 0.6.0) against
# the share of the component left by the other set, and the shared fraction `b` with Moran
# spectral randomisation of env, with moves of env on the grid, and with free permutations of
# env, which ignore its autocorrelation. Sites form a complete grid of 8 by 8; `space` holds
# the six broadest Moran's eigenvector maps; species have Gaussian responses to a gradient
# that combines env (weight alpha) and a smooth field in the span of `space` (weight beta).
# Seeds 5001 to 5500 are the holdout set: the choice of statistic and of null models was made
# with seeds 1001 to 1200.
library(comdynamite)
library(parallel)
n_rep <- 500
n_perm <- 199
seeds <- 5000 + seq_len(n_rep)
cores <- min(6, max(1, detectCores() - 2))

side <- 8
n_sites <- side^2
coords <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
mems <- as.matrix(adespatial::dbmem(coords, silent = TRUE))
space <- mems[, 1:6]
colnames(space) <- paste0("MEM", 1:6)
torus <- torus_permutations(side)

field <- function(cols) as.vector(scale(mems[, cols, drop = FALSE] %*% rnorm(length(cols))))
make_env <- function(kind) {
  switch(kind,
         iid = rnorm(n_sites),
         broad = sqrt(0.8) * field(1:6) + sqrt(0.2) * rnorm(n_sites),
         fine = sqrt(0.8) * field(7:15) + sqrt(0.2) * rnorm(n_sites))
}
# with `hidden`, the composition follows another trait than the one analysed, so the null
# hypothesis of the permutation of species is true
generate <- function(alpha, beta, env_kind, k = 20, hidden = FALSE) {
  repeat {
    env <- make_env(env_kind)
    f <- field(1:3)
    g <- if (alpha == 0 && beta == 0) rep(0, n_sites) else {
      (alpha * env + beta * f) / sqrt(alpha^2 + beta^2)
    }
    gradient <- 50 + 20 * g
    trait <- setNames(runif(k, 0, 100), paste0("sp", seq_len(k)))
    driver <- if (hidden) setNames(runif(k, 0, 100), names(trait)) else trait
    lambda <- 5 * exp(-outer(gradient, driver, "-")^2 / (2 * 15^2)) + 0.3
    comm <- matrix(rpois(length(lambda), lambda), nrow(lambda), dimnames = list(NULL, names(trait)))
    if (all(rowSums(comm) > 0) && sum(colSums(comm) > 0) >= 5) break
  }
  list(comm = comm, trait = trait, env = cbind(e = env))
}

scenarios <- data.frame(
  cell = c("none", "space_iid", "space_broad", "env_iid", "env_broad", "both_broad",
           "weak_env_broad", "space_fine", "hidden_env_broad"),
  alpha = c(0, 0, 0, 1, 1, 1, 0.5, 0, 1),
  beta = c(0, 1, 1, 0, 0, 1, 1, 1, 0),
  env_kind = c("iid", "iid", "broad", "iid", "broad", "broad", "broad", "fine", "broad"),
  hidden = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE)

# the test of `a` and `c` up to 0.6.0: the fraction of the total, with the same permutations
old_fraction_test <- function(comm, trait, env, space, n_perm) {
  w <- rowSums(comm)
  std <- function(x) comdynamite:::weighted_standardise(x, w)
  ess <- function(x, y) comdynamite:::explained_ss(x, y, w)
  response <- function(trait) comdynamite:::trait_response(comm, as.matrix(trait), "mean")
  y <- response(trait)
  fractions <- function(y, env, space) {
    env <- std(env)
    space <- std(space)
    abc <- ess(cbind(env, space), y)
    c(a = abc - ess(space, y), c = abc - ess(env, y)) / sum(comm)
  }
  observed <- fractions(y, env, space)
  n <- nrow(comm)
  site <- replicate(n_perm, {
    i <- sample.int(n)
    fractions(y, env[i, , drop = FALSE], space[i, , drop = FALSE])
  })
  species <- replicate(n_perm, fractions(response(trait[sample.int(length(trait))]), env, space))
  p <- function(null, f) comdynamite:::permutation_p(null[f, ], observed[[f]])
  out <- c(a_site = p(site, "a"), a_species = p(species, "a"), c_site = p(site, "c"),
           c_species = p(species, "c"))
  c(out, a_max = max(out[c("a_site", "a_species")]), c_max = max(out[c("c_site", "c_species")]))
}

pick <- function(out, prefix, fields = c("a", "c", "b")) {
  f <- out$fractions
  rownames(f) <- f$fraction
  p <- unlist(lapply(fields, function(x) {
    setNames(f[x, c("p_site", "p_species", "p_max")], paste0(x, c("_site", "_species", "_max")))
  }))
  p <- c(p, r2_ab = f["ab", "r2"], r2_bc = f["bc", "r2"], r2_b = f["b", "r2"],
         b_tested = as.numeric(out$shared$tested), p_env_space = out$shared$p_env_space)
  setNames(p, paste0(prefix, "_", names(p)))
}

one_replicate <- function(seed, sc) {
  set.seed(seed)
  d <- generate(sc$alpha, sc$beta, sc$env_kind, hidden = sc$hidden)
  new <- trait_partition(d$comm, d$trait, d$env, space, n_perm = n_perm, shared_null = "msr",
                         coords = coords)
  old <- old_fraction_test(d$comm, d$trait, d$env, space, n_perm)
  moves <- torus[sample.int(nrow(torus), n_perm), ]
  grid <- trait_partition(d$comm, d$trait, d$env, space, n_perm = n_perm,
                          site_permutations = moves, shared_null = "sites")
  free <- t(replicate(n_perm, sample.int(n_sites)))
  naive <- trait_partition(d$comm, d$trait, d$env, space, n_perm = n_perm,
                           site_permutations = free, shared_null = "sites")
  c(seed = seed, pick(new, "new"), setNames(old, paste0("old_", names(old))),
    pick(grid, "grid"), pick(naive, "naive", "b")[c("naive_b_site", "naive_b_tested")])
}

RNGkind("L'Ecuyer-CMRG")
results <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  sc <- scenarios[i, ]
  rows <- mclapply(seeds, function(s) {
    tryCatch(one_replicate(s, sc), error = function(e) c(seed = s, error = conditionMessage(e)))
  }, mc.cores = cores)
  failed <- vapply(rows, function(r) "error" %in% names(r), logical(1))
  if (any(failed)) print(table(vapply(rows[failed], function(r) r[["error"]], "")))
  cbind(cell = sc$cell, as.data.frame(do.call(rbind, rows[!failed])))
}))

# rejection rates at the 5% level with their standard errors; untestable cases (b when env is
# unrelated to space) are left out of the rate and counted in `untested`
rate <- function(p) {
  p <- p[!is.na(p)]
  r <- mean(p <= 0.05)
  c(rate = r, se = sqrt(r * (1 - r) / length(p)), n = length(p))
}
tests <- grep("_(site|species|max)$", names(results), value = TRUE)
rates <- do.call(rbind, lapply(split(results, results$cell), function(cell) {
  data.frame(cell = cell$cell[1], test = tests, t(sapply(tests, function(t) rate(cell[[t]]))),
             row.names = NULL)
}))
means <- aggregate(cbind(r2_ab = new_r2_ab, r2_bc = new_r2_bc, r2_b = new_r2_b,
                         untested = 1 - new_b_tested, untested_grid = 1 - grid_b_tested) ~ cell,
                   results, mean)
saveRDS(list(rates = rates, means = means, scenarios = scenarios, n_rep = n_rep, n_perm = n_perm,
             seeds = range(seeds)),
        "data-raw/validation_fractions.rds")
