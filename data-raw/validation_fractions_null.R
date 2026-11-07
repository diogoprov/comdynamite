# Restricted null models for the tests of `a` and `c` of trait_partition() with
# `site_permutations`: which set the moves of the torus should move. Two stages, each on its
# own generator, so the choice is made on one and checked on the other.
#
# Stage 1 (tuning): the generator of validation_fractions.R (grid of 8 by 8, Gaussian niches,
# 20 species, 6 Moran's eigenvector maps), seeds 1001 to 1300. The statistics of `a` and `c`
# (share of the component left by the other set) are computed with the package internals for
# every candidate from the same data and the same draws: free permutations of both sets
# together (the default of the function), of `env` alone and of `space` alone, and moves of the
# torus of both together (the restricted null up to 0.6.1), of `env` alone and of `space`
# alone. The permutation of species is the same for all, and `_max` is the larger p-value.
#
# Stage 2 (holdout): the generator of validation_fractions_check.R (counts from Gaussian niches
# on fields of exponential covariance, two env columns, grid of 10 by 10, 8 maps, 25 species),
# seeds 9001 to 9200, run through trait_partition() itself with free permutations and with
# `torus_permutations()`, which from this version move each set alone, plus the old restricted
# null (both sets together) from the internals. The scenario `env_auto_hidden` drives the
# species by `env` through a trait other than the one analysed, so the permutation of species
# holds its null hypothesis there.
#
# Rejection rates at the 5% level with standard errors at the level of the data set. About
# 10 min on 6 cores, half for each stage. Run from the root of the package, with the package
# installed.
library(comdynamite)
library(parallel)
n_perm <- 199
cores <- min(6, max(1, detectCores() - 2))
RNGkind("L'Ecuyer-CMRG")

rate <- function(p) {
  p <- p[!is.na(p)]
  r <- mean(p <= 0.05)
  c(rate = r, se = sqrt(r * (1 - r) / length(p)), n = length(p))
}
rate_table <- function(results) {
  tests <- setdiff(names(results), c("cell", "seed"))
  do.call(rbind, lapply(split(results, results$cell), function(cell) {
    data.frame(cell = cell$cell[1], test = tests, t(sapply(tests, function(t) rate(cell[[t]]))),
               row.names = NULL)
  }))
}
run_cells <- function(scenarios, seeds, one) {
  do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    sc <- scenarios[i, , drop = FALSE]
    rows <- mclapply(seeds, function(s) {
      tryCatch(one(s, sc), error = function(e) c(seed = s, error = conditionMessage(e)))
    }, mc.cores = cores)
    failed <- vapply(rows, function(r) "error" %in% names(r), logical(1))
    if (any(failed)) print(table(vapply(rows[failed], function(r) r[["error"]], "")))
    cbind(cell = sc$cell, as.data.frame(do.call(rbind, rows[!failed])))
  }))
}
# statistics of `a` and `c` for a response y with site weights w, with the rows of env and of
# space that a draw takes
stats_ac <- function(y, w, env, space, i_env, i_space) {
  comp <- sum(comdynamite:::column_ss(y, w))
  e <- comdynamite:::weighted_standardise(env[i_env, , drop = FALSE], w)
  s <- comdynamite:::weighted_standardise(space[i_space, , drop = FALSE], w)
  abc <- comdynamite:::explained_ss(cbind(e, s), y, w)
  ab <- comdynamite:::explained_ss(e, y, w)
  bc <- comdynamite:::explained_ss(s, y, w)
  c(a = comdynamite:::share_left(abc - bc, comp - bc, comp),
    c = comdynamite:::share_left(abc - ab, comp - ab, comp))
}
p_value <- comdynamite:::permutation_p

## Stage 1: tuning ---------------------------------------------------------------------------
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
tuning_scenarios <- data.frame(
  cell = c("none", "space_iid", "space_broad", "env_iid", "env_broad", "both_broad",
           "weak_env_broad", "space_fine", "hidden_env_broad"),
  alpha = c(0, 0, 0, 1, 1, 1, 0.5, 0, 1),
  beta = c(0, 1, 1, 0, 0, 1, 1, 1, 0),
  env_kind = c("iid", "iid", "broad", "iid", "broad", "broad", "broad", "fine", "broad"),
  hidden = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE)

tuning_replicate <- function(seed, sc) {
  set.seed(seed)
  d <- generate(sc$alpha, sc$beta, sc$env_kind, hidden = sc$hidden)
  comm <- d$comm[, colSums(d$comm) > 0, drop = FALSE]
  trait <- d$trait[colnames(comm)]
  w <- rowSums(comm)
  id <- seq_len(nrow(comm))
  y <- comdynamite:::trait_response(comm, as.matrix(trait), "mean")
  obs <- stats_ac(y, w, d$env, space, id, id)
  moves <- torus[sample.int(nrow(torus), n_perm), ]
  free <- t(replicate(n_perm, sample.int(nrow(comm))))
  draws <- list(
    free_both = function(r) stats_ac(y, w, d$env, space, free[r, ], free[r, ]),
    free_env = function(r) stats_ac(y, w, d$env, space, free[r, ], id),
    free_space = function(r) stats_ac(y, w, d$env, space, id, free[r, ]),
    torus_both = function(r) stats_ac(y, w, d$env, space, moves[r, ], moves[r, ]),
    torus_env = function(r) stats_ac(y, w, d$env, space, moves[r, ], id),
    torus_space = function(r) stats_ac(y, w, d$env, space, id, moves[r, ]))
  null <- lapply(draws, function(f) vapply(seq_len(n_perm), f, numeric(2)))
  species <- vapply(seq_len(n_perm), function(r) {
    shuffled <- trait[sample.int(length(trait))]
    stats_ac(comdynamite:::trait_response(comm, as.matrix(shuffled), "mean"), w, d$env, space,
             id, id)
  }, numeric(2))
  p_species <- c(a = p_value(species["a", ], obs[["a"]]), c = p_value(species["c", ], obs[["c"]]))
  out <- c(seed = seed, a_species = p_species[["a"]], c_species = p_species[["c"]])
  for (nm in names(null)) {
    for (f in c("a", "c")) {
      p <- p_value(null[[nm]][f, ], obs[[f]])
      out[paste0(f, "_", nm)] <- p
      out[paste0(f, "_", nm, "_max")] <- max(p, p_species[[f]])
    }
  }
  out
}
tuning <- run_cells(tuning_scenarios, 1000 + seq_len(300), tuning_replicate)

## Stage 2: holdout --------------------------------------------------------------------------
xy <- as.matrix(expand.grid(1:10, 1:10))
n <- 100
D <- as.matrix(dist(xy))
U <- chol(exp(-D / 3) + diag(1e-8, n))
fld <- function() as.vector(scale(drop(crossprod(U, rnorm(n)))))
mem <- suppressWarnings(as.matrix(adespatial::dbmem(xy, MEM.autocor = "positive",
                                                    silent = TRUE)))[, 1:8]
rownames(xy) <- paste0("s", 1:n)
rownames(mem) <- rownames(xy)
torus10 <- torus_permutations(10)
gen <- function(kind) {
  iid <- kind %in% c("null_iid", "space_only_iid")
  e1 <- if (iid) rnorm(n) else fld()
  e2 <- if (iid) rnorm(n) else fld()
  hid <- fld()
  driver <- switch(kind, null_iid = rnorm(n), space_only_iid = hid, space_only_auto = hid,
                   env_auto = e1, env_auto_hidden = e1, both = e1 + hid)
  tr <- setNames(rnorm(25), paste0("sp", 1:25))
  analysed <- if (kind == "env_auto_hidden") setNames(rnorm(25), names(tr)) else tr
  mu <- exp(1 - outer(as.vector(scale(driver)), tr, function(g, t) (g - t)^2 / 2))
  comm <- matrix(rpois(length(mu), 3 * mu), n, dimnames = list(rownames(xy), names(tr)))
  list(comm = comm, tr = analysed, env = `rownames<-`(cbind(e1, e2), rownames(xy)))
}
holdout_scenarios <- data.frame(
  cell = c("null_iid", "space_only_iid", "space_only_auto", "env_auto", "env_auto_hidden", "both"),
  stringsAsFactors = FALSE)

pick <- function(out, prefix) {
  f <- out$fractions
  rownames(f) <- f$fraction
  p <- unlist(lapply(c("a", "c"), function(x) {
    setNames(f[x, c("p_site", "p_species", "p_max")], paste0(x, c("_site", "_species", "_max")))
  }))
  setNames(p, paste0(prefix, "_", names(p)))
}
holdout_replicate <- function(seed, sc) {
  set.seed(seed)
  g <- gen(sc$cell)
  if (!all(rowSums(g$comm) > 0)) stop("empty site")
  comm <- g$comm[, colSums(g$comm) > 0]
  trait <- g$tr[colnames(comm)]
  free <- trait_partition(comm, trait, g$env, mem, n_perm = n_perm)
  moves <- torus10[sample.int(nrow(torus10), n_perm), ]
  grid <- trait_partition(comm, trait, g$env, mem, n_perm = n_perm, site_permutations = moves)
  # the restricted null up to 0.6.1: both sets moved together
  w <- rowSums(comm)
  id <- seq_len(n)
  y <- comdynamite:::trait_response(comm, as.matrix(trait), "mean")
  obs <- stats_ac(y, w, g$env, mem, id, id)
  both <- vapply(seq_len(n_perm), function(r) stats_ac(y, w, g$env, mem, moves[r, ], moves[r, ]),
                 numeric(2))
  p_both <- c(a = p_value(both["a", ], obs[["a"]]), c = p_value(both["c", ], obs[["c"]]))
  species <- grid$fractions$p_species[match(c("a", "c"), grid$fractions$fraction)]
  c(seed = seed, pick(free, "free"), pick(grid, "grid"), old_a_site = p_both[["a"]],
    old_a_max = max(p_both[["a"]], species[1]), old_c_site = p_both[["c"]],
    old_c_max = max(p_both[["c"]], species[2]))
}
holdout <- run_cells(holdout_scenarios, 9000 + seq_len(200), holdout_replicate)

saveRDS(list(tuning = rate_table(tuning), holdout = rate_table(holdout),
             tuning_scenarios = tuning_scenarios, holdout_scenarios = holdout_scenarios,
             n_perm = n_perm, n_rep = c(tuning = 300, holdout = 200),
             seeds = list(tuning = c(1001, 1300), holdout = c(9001, 9200))),
        "data-raw/validation_fractions_null.rds")
