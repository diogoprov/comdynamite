# ABC with imperfect detection (observe_detection()), in the design of data-raw/abc_recovery.R:
# tree of 30 species, grid of 25 sites, the same priors and model, series of five surveys. Each
# individual is detected in each of three visits with probability `detection`, drawn on the
# log scale between 0.02 and 0.5, and the record keeps the largest count. Questions set before
# running:
# (i) with detection in the reference table, are the four parameters still recovered
#     (R2 >= 0.5 and coverage 85-95%, as before), and is detection itself recovered?
# (ii) if detection is ignored, fitting the table with perfect detection of abc_recovery.R to
#     the same imperfectly detected pseudo-observations, how biased are the estimates and how
#     often do the intervals miss the truth?
# Run from the root of the package after abc_recovery.R. N_REF, N_TEST, CORES and NTREE change
# the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "20000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_detection.rds")

set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
prior <- list(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
              lambda = c(0, 1), log_detection = log10(c(0.02, 0.5)))
one_run <- function(scenario, shared) {
  sim <- simulate_metacommunity(shared$landscape,
                                species_traits(shared$tree, scenario$lambda,
                                               sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                                competition_matrix(shared$tree, "phylo",
                                                   inter_max = scenario$inter_max),
                                dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                                record_steps = seq(800, 1000, by = 50))
  observe_detection(sim, 10^scenario$log_detection, visits = 3)
}
future::plan(future::multicore, workers = cores)
started <- Sys.time()
table <- abc_table(prior, one_run, n = n_ref + n_test, setup = function() fixed,
                   summarise = function(sim, shared) summary_stats(sim, n_null = 49),
                   seed = 2035)
future::plan(future::sequential)
minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

started <- Sys.time()
is_test <- seq_len(nrow(table)) > n_ref
train <- table[!is_test, ]
attributes(train)[c("prior", "parameters", "statistics")] <-
  attributes(table)[c("prior", "parameters", "statistics")]
pseudo <- as.data.frame(table)[is_test, ]
four <- setdiff(names(prior), "log_detection")
score <- function(fit, pseudo, parameters, label) {
  ok <- rowSums(!is.finite(as.matrix(pseudo[fit$statistics]))) == 0
  d <- pseudo[ok, ]
  est <- predict(fit, d, ncores = cores)
  est <- est[est$parameter %in% parameters, ]
  est$truth <- unlist(lapply(seq_len(nrow(d)), function(i) d[i, parameters]))
  do.call(rbind, lapply(parameters, function(p) {
    e <- est[est$parameter == p, ]
    range <- diff(prior[[p]])
    covered <- e$truth >= e$low & e$truth <= e$high
    data.frame(fit = label, parameter = p, n = nrow(e),
               r2 = 1 - sum((e$truth - e$estimate)^2) / sum((e$truth - mean(e$truth))^2),
               bias = mean(e$estimate - e$truth) / range, coverage = mean(covered),
               width = mean(e$high - e$low) / range,
               recoverable = 1 - sum((e$truth - e$estimate)^2) /
                 sum((e$truth - mean(e$truth))^2) >= 0.5 &
                 mean(covered) >= 0.85 & mean(covered) <= 0.95)
  }))
}
with_detection <- abc_fit(train, ntree = ntree, ncores = cores)
perfect <- readRDS("data-raw/abc_recovery.rds")$table
perfect_fit <- abc_fit(perfect, ntree = ntree, ncores = cores)
results <- rbind(score(with_detection, pseudo, names(prior), "detection in the table"),
                 score(perfect_fit, pseudo, four, "detection ignored"))
importance <- abc_importance(with_detection)
minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

saveRDS(list(results = results, importance = importance, table = table,
             minutes = c(simulation = minutes_sim, forests = minutes_fit),
             settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree)), output)
print(results, digits = 3, row.names = FALSE)
print(do.call(rbind, lapply(split(importance, importance$parameter), head, 3)), digits = 3)
cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
