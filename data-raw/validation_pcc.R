# Checks pcc_test() against Table 1 of Peres-Neto & Kembel (2015): rejection
# rates at the 5% level with 100 sites, 50 species, one Yule tree, incidence
# data, 99 permutations and 1000 replicates. Their tree is not available, so
# power and the species test in the phylogeny-only scenario are expected to
# differ. Run from the package root; writes data-raw/validation_pcc.rds.
library(comdynamite)

RNGkind("L'Ecuyer-CMRG")
set.seed(2015)
tree <- ape::rphylo(50, birth = 1, death = 0)
design <- expand.grid(scenario = c("neither", "environment", "phylogeny", "both"),
                      tolerance = c(5, 10), stringsAsFactors = FALSE)

one_replicate <- function(scenario, tolerance) {
  traits <- species_traits(tree, lambda = as.numeric(scenario %in% c("phylogeny", "both")),
                           optima_range = c(-1, 101), sigma_mean = tolerance, sigma_sd = 10,
                           rmax_range = c(0, 30))
  env <- runif(100, 0, 100)
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  used <- if (scenario %in% c("environment", "both")) env else runif(100, 0, 100)
  keep <- rowSums(comm) > 0
  out <- pcc_test(comm[keep, ], tree, used[keep], n_perm = 99)
  c(max = out$p_max, site = out$p_site, taxa = out$p_species)
}

one_cell <- function(i) {
  p <- replicate(1000, one_replicate(design$scenario[i], design$tolerance[i]))
  c(rowMeans(p <= 0.05))
}

rates <- do.call(rbind, parallel::mclapply(seq_len(nrow(design)), one_cell, mc.cores = 8))
results <- cbind(design, rates)
saveRDS(results, "data-raw/validation_pcc.rds", version = 2)
print(results, digits = 3)
