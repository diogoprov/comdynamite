# Cross-validation of ABC across sites (spatial_folds(), fold_summarise(), abc_cv()), in the
# design of data-raw/abc_recovery.R: tree of 30 species, grid of 25 sites, the same priors and
# model with phylogenetic competition. Five folds, as spatial blocks or at random. Questions set
# before running:
# (i) with data from the model itself (pseudo-observations left out of the table), is the
#     coverage of the 90% predictive interval of left-out statistics near 90% for both kinds of
#     folds? This checks abc_cv() itself.
# (ii) with data from a wrong model, in which species track a hidden environmental field with
#     spatial autocorrelation instead of the measured one, do spatial blocks flag the misfit
#     more often than random folds? Expected yes, since random test sites resemble their
#     training neighbours.
# Run from the root of the package. N_REF, N_TEST, N_WRONG, CORES and NTREE change the size;
# OUTPUT the file. SIDE and TORUS change the grid (5 and TRUE by default). Item (iv) of block C
# of dev/ROADMAP.md repeated the study with SIDE=5 TORUS=FALSE, SIDE=10 TORUS=TRUE and SIDE=10
# TORUS=FALSE, and again the default with the statistics of 2026-11-24
# (OUTPUT=data-raw/abc_cv_<side>_<torus or flat>.rds), five folds of 20 sites in
# the larger grid, to learn whether the torus or the small grid hid the advantage of blocks.
# Criterion set before running: blocks are worth recommending if, in the grid of 100 sites
# without torus, they reject the wrong model at least 3 percentage points more often than
# random folds (paired over the same data sets, with the 95% interval of the difference above
# zero), while the coverage with the right model stays at 0.85 or more.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "100"))
n_wrong <- as.integer(Sys.getenv("N_WRONG", "100"))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_cv_5_torus.rds")
side <- as.integer(Sys.getenv("SIDE", "5"))
torus <- as.logical(Sys.getenv("TORUS", "TRUE"))

set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(30),
              landscape = landscape_grid(side = side, torus = torus))
prior <- list(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
              lambda = c(0, 1))
simulate <- function(scenario, shared) {
  landscape <- shared$landscape
  if (!is.null(shared$hidden)) landscape$env <- shared$hidden
  simulate_metacommunity(landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                         competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max),
                         dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                         record_steps = seq(800, 1000, by = 50))
}
folds <- list(B = spatial_folds(fixed$landscape$coords, k = 5, method = "blocks"),
              R = spatial_folds(fixed$landscape$coords, k = 5, method = "random"))
stats_of <- function(x, shared) summary_stats(x, n_null = 19)
by_fold <- lapply(folds, fold_summarise, summarise = stats_of)
# the statistics of both kinds of folds from one simulation; the observed data always carry the
# measured environment, even when the species tracked a hidden one
summarise <- function(sim, shared) {
  data <- as_mc_data(sim)
  data$env <- data.frame(env = fixed$landscape$env)
  unlist(lapply(names(by_fold), function(kind) {
    s <- by_fold[[kind]](data, shared)
    stats::setNames(s, paste0(kind, names(s)))
  }))
}

future::plan(future::multicore, workers = cores)
started <- Sys.time()
table <- abc_table(prior, simulate, n = n_ref + n_test, setup = function() fixed,
                   summarise = summarise, seed = 2032)
# wrong model: each data set has its own hidden field, from a grid drawn with another seed
cell_seed(2033)
wrong_parameters <- as.data.frame(lapply(prior, function(p) stats::runif(n_wrong, p[1], p[2])))
wrong <- do.call(rbind, progress_map(seq_len(n_wrong), function(i) {
  cell_seed(4000 + i)
  shared <- c(fixed, list(hidden = landscape_grid(side = side, torus = torus)$env))
  summarise(simulate(as.list(wrong_parameters[i, ]), shared), shared)
}, "Wrong model"))
future::plan(future::sequential)
minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

reference <- seq_len(nrow(table)) <= n_ref
started <- Sys.time()
results <- do.call(rbind, lapply(names(folds), function(kind) {
  columns <- grep(paste0("^", kind, "f"), names(table), value = TRUE)
  t <- as.data.frame(table)[reference, c(names(prior), columns)]
  names(t) <- sub(paste0("^", kind), "", names(t))
  observed <- list(right = as.data.frame(table)[!reference, columns],
                   wrong = as.data.frame(wrong)[, columns])
  do.call(rbind, lapply(names(observed), function(truth) {
    o <- observed[[truth]]
    names(o) <- sub(paste0("^", kind), "", names(o))
    cv <- abc_cv(t, o, parameters = names(prior), ntree = ntree, ncores = cores, seed = 2034)
    cbind(folds = c(B = "blocks", R = "random")[[kind]], truth = truth,
          as.data.frame(unclass(cv)))
  }))
}))
minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

by_set <- aggregate(cbind(coverage = inside, rejected = p < 0.05) ~ folds + truth + observation,
                    results, mean)
overall <- aggregate(cbind(coverage, rejected) ~ folds + truth, by_set, mean)
# rejection with blocks minus rejection with random folds, paired over the same data sets
paired <- do.call(rbind, lapply(c("right", "wrong"), function(truth) {
  b <- by_set[by_set$folds == "blocks" & by_set$truth == truth, ]
  r <- by_set[by_set$folds == "random" & by_set$truth == truth, ]
  d <- b$rejected[order(b$observation)] - r$rejected[order(r$observation)]
  se <- stats::sd(d) / sqrt(length(d))
  data.frame(truth = truth, difference = mean(d), low = mean(d) - 1.96 * se,
             high = mean(d) + 1.96 * se)
}))
by_statistic <- aggregate(cbind(coverage = inside, rejected = p < 0.05) ~ folds + truth + statistic,
                          results, mean)
saveRDS(list(results = results, overall = overall, paired = paired, by_statistic = by_statistic,
             folds = folds,
             wrong_parameters = wrong_parameters,
             minutes = c(simulation = minutes_sim, forests = minutes_fit),
             settings = list(n_ref = n_ref, n_test = n_test, n_wrong = n_wrong, ntree = ntree,
                             side = side, torus = torus)),
        output, compress = "xz")
print(overall, digits = 3)
print(paired, digits = 3)
print(by_statistic[order(by_statistic$truth, by_statistic$statistic, by_statistic$folds), ],
      digits = 3, row.names = FALSE)
cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
