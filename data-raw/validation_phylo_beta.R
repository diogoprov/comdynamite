# Checks phylo_beta_test() against Table 1 of Peres-Neto & Kembel (2015), with
# pcc_test() on the same data: rejection rates at the 5% level with 100 sites, 50
# species, one Yule tree, incidence data, one gradient, 99 permutations and 1000
# replicates. Their tree is not available, so power and the species test in the
# phylogeny-only scenario are expected to differ. Run from the package root;
# writes data-raw/validation_phylo_beta.rds.
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
  beta <- phylo_beta_test(comm[keep, ], tree, used[keep], n_perm = 99)
  pcc <- pcc_test(comm[keep, ], tree, used[keep], n_perm = 99)
  c(beta_max = beta$p_max, beta_site = beta$p_site, beta_taxa = beta$p_species,
    pcc_max = pcc$p_max, pcc_site = pcc$p_site, pcc_taxa = pcc$p_species, beta_r2 = beta$r2)
}

one_cell <- function(i) {
  p <- replicate(1000, one_replicate(design$scenario[i], design$tolerance[i]))
  rejected <- p[rownames(p) != "beta_r2", ] <= 0.05
  c(rowMeans(rejected), beta_r2 = mean(p["beta_r2", ]))
}

rates <- do.call(rbind, parallel::mclapply(seq_len(nrow(design)), one_cell, mc.cores = 8))
results <- cbind(design, rates)
saveRDS(results, "data-raw/validation_phylo_beta.rds", version = 2)
print(results, digits = 3)
