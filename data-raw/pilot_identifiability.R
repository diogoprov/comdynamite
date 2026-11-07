# Pilot for ABC: how much of four parameters can be recovered from summary_stats(), with one
# survey or with a series of five? Criterion set before running: a parameter is recoverable with
# R2 >= 0.5 and coverage of the 90% interval between 85% and 95%, in pseudo-observations left
# out of the reference table. The criterion is a choice of the authors, not from the literature.
# Run from the root of the package, with the package abc installed. N_REF and N_TEST change the
# size, CORES the number of workers and TOLS the tolerances. The pilot as planned is the default.
# After seeing that coverage rose with the tolerance, a larger table was run as a sensitivity
# analysis, not part of the plan: N_REF=50000 TOLS=0.002,0.005,0.01,0.05
# OUTPUT=data-raw/pilot_identifiability_large.rds
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "5000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "6"))
output <- Sys.getenv("OUTPUT", "data-raw/pilot_identifiability.rds")
tolerances <- as.numeric(strsplit(Sys.getenv("TOLS", "0.01,0.02,0.05"), ",")[[1]])

set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
n <- n_ref + n_test
parameters <- data.frame(log_dispersal = runif(n, log10(0.001), log10(0.5)),
                         inter_max = runif(n, 0, 1), sigma_mean = runif(n, 2, 30),
                         lambda = runif(n, 0, 1))
limits <- rbind(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
                lambda = c(0, 1))

one_run <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                         competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max),
                         dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                         record_steps = seq(800, 1000, by = 50))
}
# the statistics that are not temporal come from the last step, so one call gives both views
both_views <- function(sim, shared) summary_stats(sim, shared$tree, n_null = 49)
future::plan(future::multicore, workers = cores)
started <- Sys.time()
table <- run_scenarios(parameters, one_run, function() fixed, both_views, n_rep = 1, seed = 2026)
future::plan(future::sequential)
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

names_parameters <- names(parameters)
temporal <- c("cv_population", "cv_site", "cv_regional", "turnover")
survey <- setdiff(names(table), c(names_parameters, "replicate", temporal))
views <- list(survey = survey, series = c(survey, temporal))

is_test <- seq_len(n) > n_ref
weighted_quantile <- function(x, w, p) {
  o <- order(x)
  x[o][pmin(findInterval(p, cumsum(w[o]) / sum(w)) + 1, length(x))]
}
evaluate <- function(view, tol) {
  statistics <- as.matrix(table[, views[[view]]])
  usable <- colMeans(is.finite(statistics)) > 0.95 & apply(statistics, 2, sd, na.rm = TRUE) > 0
  statistics <- statistics[, usable, drop = FALSE]
  complete <- rowSums(!is.finite(statistics)) == 0
  reference <- !is_test & complete
  test <- which(is_test & complete)
  rows <- lapply(test, function(i) {
    fit <- suppressWarnings(abc::abc(target = statistics[i, ],
                                     param = table[reference, names_parameters],
                                     sumstat = statistics[reference, ], tol = tol,
                                     method = "loclinear", transf = rep("logit", 4),
                                     logit.bounds = limits))
    w <- fit$weights
    q <- apply(fit$adj.values, 2, weighted_quantile, w = w, p = c(0.05, 0.5, 0.95))
    data.frame(parameter = names_parameters, truth = unlist(table[i, names_parameters]),
               low = q[1, ], estimate = q[2, ], high = q[3, ], row = i)
  })
  out <- do.call(rbind, rows)
  out$view <- view
  out$tol <- tol
  attr(out, "statistics") <- colnames(statistics)
  attr(out, "used") <- c(reference = sum(reference), test = length(test))
  out
}
runs <- list()
for (view in names(views)) for (tol in tolerances) {
  runs[[paste(view, tol)]] <- evaluate(view, tol)
}
estimates <- do.call(rbind, runs)

# R2 of prediction, with a bootstrap over pseudo-observations for its interval
r2 <- function(d) 1 - sum((d$truth - d$estimate)^2) / sum((d$truth - mean(d$truth))^2)
summarise <- function(d) {
  range <- diff(limits[d$parameter[1], ])
  boot <- replicate(1000, r2(d[sample.int(nrow(d), replace = TRUE), ]))
  covered <- d$truth >= d$low & d$truth <= d$high
  data.frame(view = d$view[1], tol = d$tol[1], parameter = d$parameter[1], n = nrow(d),
             r2 = r2(d), r2_low = quantile(boot, 0.025), r2_high = quantile(boot, 0.975),
             bias = mean(d$estimate - d$truth) / range, coverage = mean(covered),
             coverage_se = sqrt(mean(covered) * (1 - mean(covered)) / nrow(d)),
             width = mean(d$high - d$low) / range,
             recoverable = r2(d) >= 0.5 & mean(covered) >= 0.85 & mean(covered) <= 0.95)
}
set.seed(1)
results <- do.call(rbind, lapply(split(estimates,
                                       list(estimates$view, estimates$tol, estimates$parameter)),
                                 summarise))
rownames(results) <- NULL
results <- results[order(results$tol, results$parameter, results$view), ]

saveRDS(list(results = results, estimates = estimates, table = if (n_ref <= 5000) table,
             minutes = minutes,
             used = lapply(runs, attr, "used"), statistics = lapply(runs, attr, "statistics")),
        output)
cat("minutes:", round(minutes, 1), "\n")
print(lapply(runs, attr, "used")[c(1, 4)])
print(results, digits = 2, row.names = FALSE)
