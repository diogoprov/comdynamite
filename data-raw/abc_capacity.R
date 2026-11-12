# Is detection still identifiable when the scale of abundance is unknown? Follows
# data-raw/abc_detection.R (same tree, grid, priors, model, three visits, largest count), with
# one more parameter: log_capacity, between -1 and 1, divides every competition coefficient by
# 10^log_capacity, so equilibrium abundances, and the carrying capacity, are multiplied by it
# (Beverton-Holt dynamics, N* = (r - 1) / alpha_ii). Recorded counts scale with detection times
# capacity, so the two can trade off. Questions set before running:
# (i) are detection and capacity recovered separately (R2 >= 0.5 and coverage 85-95%), or only
#     their product, log_detection + log_capacity?
# (ii) do the four parameters of abc_recovery.R keep their recovery?
# Run from the root of the package. N_REF, N_TEST, CORES and NTREE change the size; OUTPUT the
# file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "20000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_capacity.rds")

set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
prior <- list(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
              lambda = c(0, 1), log_detection = log10(c(0.02, 0.5)),
              log_capacity = c(-1, 1))
one_run <- function(scenario, shared) {
  sim <- simulate_metacommunity(shared$landscape,
                                species_traits(shared$tree, scenario$lambda,
                                               sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                                competition_matrix(shared$tree, "phylo",
                                                   inter_max = scenario$inter_max,
                                                   scale = 0.05 / 10^scenario$log_capacity),
                                dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                                record_steps = seq(800, 1000, by = 50))
  observe_detection(sim, 10^scenario$log_detection, visits = 3)
}
future::plan(future::multicore, workers = cores)
started <- Sys.time()
table <- abc_table(prior, one_run, n = n_ref + n_test, setup = function() fixed,
                   summarise = function(sim, shared) summary_stats(sim, n_null = 49),
                   seed = 2036)
future::plan(future::sequential)
minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

started <- Sys.time()
is_test <- seq_len(nrow(table)) > n_ref
table$log_product <- table$log_detection + table$log_capacity
prior$log_product <- c(sum(sapply(prior[c("log_detection", "log_capacity")], `[`, 1)),
                       sum(sapply(prior[c("log_detection", "log_capacity")], `[`, 2)))
parameters <- names(prior)
train <- as.data.frame(table)[!is_test, ]
pseudo <- as.data.frame(table)[is_test, ]
statistics <- attr(table, "statistics")
fit <- abc_fit(train, parameters, statistics, ntree = ntree, ncores = cores)
ok <- rowSums(!is.finite(as.matrix(pseudo[fit$statistics]))) == 0
d <- pseudo[ok, ]
est <- predict(fit, d, ncores = cores)
est$truth <- unlist(lapply(seq_len(nrow(d)), function(i) d[i, parameters]))
results <- do.call(rbind, lapply(parameters, function(p) {
  e <- est[est$parameter == p, ]
  range <- diff(prior[[p]])
  covered <- e$truth >= e$low & e$truth <= e$high
  r2 <- 1 - sum((e$truth - e$estimate)^2) / sum((e$truth - mean(e$truth))^2)
  data.frame(parameter = p, n = nrow(e), r2 = r2, bias = mean(e$estimate - e$truth) / range,
             coverage = mean(covered), width = mean(e$high - e$low) / range,
             recoverable = r2 >= 0.5 & mean(covered) >= 0.85 & mean(covered) <= 0.95)
}))
# correlation between the errors of detection and capacity: strongly negative if they trade off
wide <- reshape(est[est$parameter %in% c("log_detection", "log_capacity"),
                    c("observation", "parameter", "estimate", "truth")],
                idvar = "observation", timevar = "parameter", direction = "wide")
error_correlation <- stats::cor(wide$estimate.log_detection - wide$truth.log_detection,
                                wide$estimate.log_capacity - wide$truth.log_capacity)
importance <- abc_importance(fit)
minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

saveRDS(list(results = results, importance = importance, table = table,
             error_correlation = error_correlation,
             minutes = c(simulation = minutes_sim, forests = minutes_fit),
             settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree)), output)
print(results, digits = 3, row.names = FALSE)
cat("correlation between the errors of detection and capacity:", round(error_correlation, 3), "\n")
print(do.call(rbind, lapply(split(importance, importance$parameter), head, 3)), digits = 3)
cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
