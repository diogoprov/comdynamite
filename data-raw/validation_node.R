# Checks node_test(). The statistics of each node were compared with the MATLAB
# function MetacommunityPhylogenetics of Leibold et al. (2010), run in Octave 11.3 on
# the data in tests/testthat/test-node.R: the largest difference was below 1e-15.
# node_contrasts() gave the same 13 contrasts as CreateNodeContrasts on the tree of
# that function. The MATLAB code has no licence, so it is not in this repository.
#
# Below, rejection rates at the 5% level per node, before any correction, with 60
# sites, 24 species, one Yule tree, incidence data, 99 permutations and 300 replicates.
# A site and a species test are each valid only under their own null hypothesis.
# Run from the package root; writes data-raw/validation_node.rds.
library(comdynamite)

RNGkind("L'Ecuyer-CMRG")
set.seed(2010)
tree <- ape::rphylo(24, birth = 1, death = 0)
nodes <- node_contrasts(tree, cherries = FALSE)

one_replicate <- function(scenario) {
  traits <- species_traits(tree, lambda = as.numeric(scenario %in% c("phylogeny", "both")),
                           optima_range = c(-1, 101), sigma_mean = 10, sigma_sd = 10,
                           rmax_range = c(0, 30))
  env <- runif(60, 0, 100)
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  used <- if (scenario %in% c("environment", "both")) env else runif(60, 0, 100)
  keep <- rowSums(comm) > 0
  out <- node_test(comm[keep, ], nodes, used[keep], n_perm = 99)$nodes
  c(max = mean(out$p_max <= 0.05, na.rm = TRUE), site = mean(out$p_site <= 0.05, na.rm = TRUE),
    taxa = mean(out$p_taxa <= 0.05, na.rm = TRUE), joint = mean(out$p_joint <= 0.05, na.rm = TRUE),
    any_node = any(out$p_max <= 0.05, na.rm = TRUE))
}

scenarios <- c("neither", "environment", "phylogeny", "both")
rates <- do.call(rbind, parallel::mclapply(scenarios, function(s) {
  rowMeans(replicate(300, one_replicate(s)))
}, mc.cores = 4))
results <- data.frame(scenario = scenarios, rates)
saveRDS(results, "data-raw/validation_node.rds", version = 2)
print(results, digits = 3)
