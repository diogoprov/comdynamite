# Recovery of four parameters by ABC with random forests (abc_fit(), abc_recovery()), in the
# design of the pilot for ABC (data-raw/pilot_identifiability.R): the same tree of 30 species,
# grid of 25 sites, priors, simulation and statistics, with one survey or a series of five. The
# pilot used rejection ABC and its large run (50 000 simulations, with Blomberg's K) recovered
# R2 of 0.38-0.98 but coverage of the 90% interval of 0.64-0.67 in the series view. Criterion
# set before running, the same as the pilot's: recoverable with R2 >= 0.5 and coverage between
# 85% and 95%, in 200 pseudo-observations left out of the table. The parameter draws differ
# from the pilot's, so the comparison is between designs, not between paired draws.
# Run from the root of the package, with abcrf installed. N_REF, N_TEST, CORES and NTREE change
# the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "50000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_recovery.rds")

set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
prior <- list(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
              lambda = c(0, 1))
one_run <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                         competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max),
                         dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                         record_steps = seq(800, 1000, by = 50))
}
future::plan(future::multicore, workers = cores)
started <- Sys.time()
table <- abc_table(prior, one_run, n = n_ref + n_test, setup = function() fixed,
                   summarise = function(sim, shared) summary_stats(sim, shared$tree, n_null = 49),
                   seed = 2027)
future::plan(future::sequential)
minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

temporal <- c("cv_population", "cv_site", "cv_regional", "turnover")
survey <- setdiff(attr(table, "statistics"), temporal)
views <- list(survey = survey, series = c(survey, temporal))
started <- Sys.time()
results <- do.call(rbind, lapply(names(views), function(view) {
  rec <- abc_recovery(table, statistics = views[[view]], n_test = n_test, ntree = ntree,
                      ncores = cores, seed = 2028)
  cbind(view = view, rec,
        recoverable = rec$r2 >= 0.5 & rec$coverage >= 0.85 & rec$coverage <= 0.95)
}))
minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))
fit <- abc_fit(table, statistics = views$series, ntree = ntree, ncores = cores)
importance <- abc_importance(fit)

saveRDS(list(results = results, importance = importance, table = table,
             minutes = c(simulation = minutes_sim, forests = minutes_fit),
             settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree)), output)
print(results, digits = 3)
print(do.call(rbind, lapply(split(importance, importance$parameter), head, 3)), digits = 3)
cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
