# The bootstrap subset of data-raw/bridge_experiments.R run again after the fix of the null
# generator of population_shifts() (slopes from the fit without displacement, species constants
# from the populations that enter the fit, precision-weighted variance of populations, sampling
# variance of the fitted set net of the push removed).
# Same cells, processes and designs as the bootstrap subset there, new simulations.
#   Rscript data-raw/bridge_bootstrap.R
# Results in data-raw/bridge_bootstrap.rds (about 20 min on 6 cores).

devtools::load_all(quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(2027)
cores <- 6

cells <- expand.grid(dispersal = c(0.001, 0.01, 0.1), extirpation = c(0, 0.01, 0.05),
                     bandwidth = c(0.5, 3), rep = 1:5)
designs <- data.frame(n_sites = c(20, 100), n_ind = c(5, 20))

subsample <- function(ind, n_sites, n_ind) {
  sites <- unique(ind$site)
  keep <- if (n_sites < length(sites)) sample(sites, n_sites) else sites
  ind[ind$site %in% keep & ind$individual <= n_ind, ]
}

one_cell <- function(i) {
  cell <- cells[i, ]
  tree <- treeshapesim::sim_tree_shape(30)
  sim <- simulate_metacommunity(landscape_grid(10, env_bandwidth = cell$bandwidth),
                                species_traits(tree), competition_matrix(tree, "phylo"),
                                dispersal = cell$dispersal, extirpation = cell$extirpation,
                                n_steps = 600)
  if (sum(sim$abundance > 0) < 10) return(NULL)
  rows <- list()
  for (cc in c(0, 3)) {
    full <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = cc,
                               population_sd = 1, intraspecific_sd = 5)
    for (d in seq_len(nrow(designs))) {
      ind <- subsample(full, designs$n_sites[d], designs$n_ind[d])
      fit <- tryCatch(suppressWarnings(population_shifts(ind, nboot = 99)),
                      error = function(e) NULL)
      if (is.null(fit)) next
      rows[[length(rows) + 1]] <- data.frame(i = i, cell, h = 0.5, c = cc, designs[d, ],
                                             as.data.frame(fit), row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

start <- Sys.time()
out <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
  tryCatch(one_cell(i), error = function(e) NULL)
}, mc.cores = cores, mc.set.seed = TRUE)
results <- list(boot = do.call(rbind, out),
                elapsed = difftime(Sys.time(), start, units = "mins"))
saveRDS(results, "data-raw/bridge_bootstrap.rds")
cat("done in", format(results$elapsed), "\n")
