# Time of the permutation tests node_test(), pcc_test(), pcps_test(), phylo_beta_test() and
# trait_partition() as the numbers of sites, species and permutations grow (rOpenSci standard
# RE5.0). Run from the root of the package. One core: the future plan is sequential and
# node_test() has cores = 1. Grid: 25, 50, 100 and 200 sites, 25, 50 and 100 species, 99 and
# 499 permutations, all 24 combinations for every function (no grid point dropped). The four
# fast functions are timed as the shortest of three runs, node_test() once. Data as in
# tests/testthat/test-extended-perm.R: a tree from treeshapesim::sim_tree_shape(), optima that
# evolve on it, a uniform gradient and incidences from simulate_static(); space for
# trait_partition() is the first three dbMEMs of a transect. Timings use
# set.seed() per grid point. The script saves data-raw/benchmark_regression.rds and prints,
# per function, the exponents of the regression of log time on log sites, log species and log
# permutations. It took about 8 minutes on one core of a laptop.
suppressMessages(devtools::load_all(".", quiet = TRUE))
future::plan(future::sequential)

make_data <- function(n_sites, n_species) {
  repeat {
    tree <- treeshapesim::sim_tree_shape(n_species)
    traits <- species_traits(tree, optima_range = c(-1, 101), sigma_mean = 15, sigma_sd = 5)
    env <- stats::runif(n_sites, 0, 100)
    comm <- simulate_static(traits, env, incidence = TRUE)$abundance
    if (all(rowSums(comm) > 0) && all(colSums(comm) > 0)) break
  }
  colnames(comm) <- tree$tip.label
  trait <- stats::setNames(traits$optimum, tree$tip.label)
  space <- as.matrix(adespatial::dbmem(cbind(seq_len(n_sites)), silent = TRUE))[, 1:3]
  rownames(space) <- rownames(comm)
  list(comm = comm, tree = tree, env = env, trait = trait, space = space,
       nodes = node_contrasts(tree))
}

# shortest elapsed time of `times` runs, less sensitive to garbage collection and other load
time_of <- function(f, times = 1) {
  min(vapply(seq_len(times), function(i) system.time(f())[["elapsed"]], numeric(1)))
}
run_all <- function(d, n, times) {
  c(node_test = time_of(function() node_test(d$comm, d$nodes, d$env, n_perm = n, cores = 1)),
    pcc_test = time_of(function() pcc_test(d$comm, d$tree, d$env, n_perm = n), times),
    pcps_test = time_of(function() pcps_test(d$comm, d$tree, d$env, n_perm = n), times),
    phylo_beta_test = time_of(function() phylo_beta_test(d$comm, d$tree, d$env, n_perm = n),
                              times),
    trait_partition = time_of(function() trait_partition(d$comm, d$trait, d$env, d$space,
                                                         n_perm = n), times))
}

# first calls load code and print citations; they are not timed
set.seed(999)
invisible(suppressMessages(run_all(make_data(20, 20), 9, 1)))

grid <- expand.grid(sites = c(25, 50, 100, 200), species = c(25, 50, 100),
                    n_perm = c(99, 499))
start <- Sys.time()
timings <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  set.seed(1000 + i)
  d <- make_data(g$sites, g$species)
  n <- g$n_perm
  seconds <- run_all(d, n, times = 3)
  message(sprintf("%d sites, %d species, %d permutations: %s", g$sites, g$species, n,
                  paste(names(seconds), round(seconds, 2), collapse = ", ")))
  data.frame(g[rep(1, length(seconds)), ], fun = names(seconds), seconds = unname(seconds),
             row.names = NULL)
}))
attr(timings, "machine") <- paste("one core of a laptop,", R.version.string)
attr(timings, "total_minutes") <- as.numeric(difftime(Sys.time(), start, units = "mins"))
saveRDS(timings, "data-raw/benchmark_regression.rds")

exponents <- do.call(rbind, lapply(split(timings, timings$fun), function(x) {
  fit <- stats::lm(log(seconds) ~ log(sites) + log(species) + log(n_perm), data = x)
  b <- stats::coef(fit)
  data.frame(fun = x$fun[1], sites = b[["log(sites)"]], species = b[["log(species)"]],
             n_perm = b[["log(n_perm)"]], r2 = summary(fit)$r.squared,
             largest = max(x$seconds))
}))
print(exponents, digits = 2, row.names = FALSE)
message("Total: ", round(attr(timings, "total_minutes"), 1), " minutes")
