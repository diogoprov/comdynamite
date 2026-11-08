# BAMM (Rabosky 2014) against trait_rates() with shifts as a Poisson process, as under the prior
# of BAMM, and each new regime with a rate of its own. BAMM is a compiled program (GPL), only
# run here; BAMMtools reads its output. Run from the root of the package with
#   BAMM=/path/to/bamm Rscript data-raw/validation_bamm.R
# CORES sets the number of workers and GENERATIONS the length of each chain.
# For each tree: Spearman and log-scale Pearson correlation between the mean rate per branch of
# the posterior and the true rate, and the marginal probability of a shift on branches with and
# without a true shift.
suppressMessages(devtools::load_all(".", quiet = TRUE))
bamm <- Sys.getenv("BAMM", "bamm")
cores <- as.integer(Sys.getenv("CORES", "6"))
generations <- as.integer(Sys.getenv("GENERATIONS", "2000000"))
RNGkind("L'Ecuyer-CMRG")
set.seed(31)

control <- function(dir, priors, seed) {
  lines <- c(
    "modeltype = trait", "treefile = tree.tre", "traitfile = traits.txt",
    "runInfoFilename = run_info.txt", "sampleFromPriorOnly = 0", "runMCMC = 1",
    "loadEventData = 0", "initializeModel = 1", "useGlobalSamplingProbability = 1",
    priors,
    "betaIsTimeVariablePrior = 0", "betaShiftInit = 0", "updateRateBetaShift = 0",
    "updateRateBetaTimeMode = 0", "useObservedMinMaxAsTraitPriors = 1",
    paste("numberOfGenerations =", generations), "mcmcOutfile = mcmc_out.txt",
    paste("mcmcWriteFreq =", generations / 1000), "eventDataOutfile = event_data.txt",
    paste("eventDataWriteFreq =", generations / 1000), "printFreq = 1000000000",
    "acceptanceResetFreq = 1000", "overwrite = 1", "updateRateEventNumber = 1",
    "updateRateEventPosition = 1", "updateRateEventRate = 1", "updateRateBeta0 = 1",
    "updateRateNodeState = 25", "initialNumberEvents = 0", "minCladeSizeForShift = 1",
    "numberOfChains = 1", paste("seed =", seed))
  writeLines(lines, file.path(dir, "control.txt"))
}

one_tree <- function(r) {
  dir <- tempfile("bamm")
  dir.create(dir)
  tree <- treeshapesim::sim_tree_shape(100)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  x <- trait_rates(tree, "shift", shift_rate = 0.3, relative = FALSE,
                   shift_size = function(n) exp(sample(c(-1, 1), n, replace = TRUE) * log(6)))
  ape::write.tree(tree, file.path(dir, "tree.tre"))
  utils::write.table(data.frame(names(x$tips), x$tips), file.path(dir, "traits.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  prior_file <- file.path(dir, "priors.txt")
  BAMMtools::setBAMMpriors(tree, traits = file.path(dir, "traits.txt"), outfile = prior_file,
                           suppressWarning = TRUE)
  priors <- grep("=", readLines(prior_file), value = TRUE)
  priors <- grep("^#", priors, value = TRUE, invert = TRUE)
  control(dir, priors, sample.int(1e6, 1))
  old <- setwd(dir)  # BAMM reads and writes in the working directory
  on.exit(setwd(old))
  status <- system2(bamm, c("-c", "control.txt"), stdout = FALSE, stderr = FALSE)
  if (status != 0 || !file.exists("event_data.txt")) {
    return(data.frame(shifts = nrow(x$shifts), spearman = NA, pearson_log = NA,
                      prob_shift = NA, prob_none = NA))
  }
  events <- BAMMtools::getEventData(tree, "event_data.txt", burnin = 0.2, type = "trait",
                                    verbose = FALSE)
  mean_rate <- BAMMtools::getMeanBranchLengthTree(events, rate = "trait")$phy
  estimate <- mean_rate$edge.length[match(tree$edge[, 2], mean_rate$edge[, 2])]
  probability <- BAMMtools::marginalShiftProbsTree(events)
  shift_probability <- probability$edge.length[match(tree$edge[, 2], probability$edge[, 2])]
  shifted <- seq_len(nrow(tree$edge)) %in% x$shifts$edge
  data.frame(shifts = nrow(x$shifts),
             spearman = stats::cor(estimate, x$edge_rate, method = "spearman"),
             pearson_log = stats::cor(log(estimate), log(x$edge_rate)),
             prob_shift = if (any(shifted)) mean(shift_probability[shifted]) else NA,
             prob_none = mean(shift_probability[!shifted]))
}

results <- do.call(rbind, parallel::mclapply(seq_len(20), one_tree, mc.cores = cores,
                                             mc.set.seed = TRUE))
cat("BAMM against trait_rates(shift_rate = 0.3, relative = FALSE), 100 tips:\n")
print(round(colMeans(results, na.rm = TRUE), 2))
cat("failed runs:", sum(is.na(results$spearman)), "of", nrow(results), "\n")
saveRDS(results, "data-raw/validation_bamm.rds")
