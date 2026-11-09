# From BAMM output to the regimes and covariates of damocles_fit(), damocles_samples() and
# damocles_partition_test(), with the example data of BAMMtools (version 2.1.12; GPL, only used
# here): the tree of whales and its BAMM event data.
# 1. Regimes per edge from each posterior sample (edata$eventVectors), and the mean branch rate
#    of speciation in each sample (getMarginalBranchRateMatrix()), both in the order of
#    tree$edge.
# 2. Communities simulated under DAMOCLES with colonisation three times faster in the regime
#    of the best shift configuration than elsewhere; then damocles_samples() over 100 posterior
#    samples, with regimes and with the branch rates as a covariate, and
#    damocles_partition_test() with the best configuration.
# RPANDA (ClaDS) was archived on CRAN; with branch-specific rates from ClaDS in the order of
# tree$edge the path is the same as for the covariate here.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(BAMMtools)
})
cores <- as.integer(Sys.getenv("CORES", "6"))
RNGkind("L'Ecuyer-CMRG")
set.seed(2024)

data(whales, events.whales)
quietly <- function(x) {
  utils::capture.output(value <- suppressMessages(x))
  value
}
edata <- quietly(getEventData(whales, events.whales, burnin = 0.1, nsamples = 100,
                              verbose = FALSE))
stopifnot(identical(edata$edge, whales$edge))
regimes <- lapply(edata$eventVectors, factor)
lambda <- getMarginalBranchRateMatrix(edata, verbose = FALSE)$lambda_branch_matrix
covariates <- lapply(seq_len(ncol(lambda)), function(i) lambda[, i])
best <- quietly(getBestShiftConfiguration(edata, expectedNumberOfShifts = 1))
best_regimes <- factor(best$eventVectors[[1]])
cat("1. regimes per posterior sample:", paste(names(table(lengths(lapply(regimes, levels)))),
                                             collapse = ", "),
    "regimes; best configuration", nlevels(best_regimes), "regimes,",
    paste(table(best_regimes), collapse = " / "), "edges\n")

height <- max(ape::node.depth.edgelength(whales))
background <- names(which.max(table(best_regimes[whales$edge[, 1] == ape::Ntip(whales) + 1])))
gamma <- ifelse(best_regimes == background, 0.5, 1.5) / height
comm <- simulate_damocles(whales, gamma, 1 / height, n = 10)
by_regime <- damocles_samples(comm, whales, regimes = regimes, vary = "gamma", cores = cores)
by_rate <- damocles_samples(comm, whales, edge_covariate = covariates, vary = "gamma",
                            cores = cores)
partition <- damocles_partition_test(comm, whales, best_regimes, n_rand = 99, vary = "gamma",
                                     cores = cores)
cat("2. regimes of 100 posterior samples\n")
print(summary(by_regime)$comparison, row.names = FALSE)
cat("   mean branch rates of speciation as covariate\n")
print(summary(by_rate)$comparison, row.names = FALSE)
print(summary(by_rate)$parameters, digits = 3)
print(partition)
saveRDS(list(by_regime = by_regime$samples, by_rate = by_rate$samples,
             partition = partition[c("statistic", "df", "p_chisq", "p_random", "null")]),
        "data-raw/damocles_bamm.rds")
