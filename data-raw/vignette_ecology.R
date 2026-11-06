# Simulations of vignettes/phylogenetic-structure.Rmd. Run from the package root
# (about 5 minutes); writes vignettes/phylogenetic-structure.rds, which the
# vignette reads so that it builds quickly. Needs picante and vegan.
library(comdynamite)
library(picante)
library(vegan)
set.seed(2026)

n_rep <- 30
scenarios <- expand.grid(dispersal = c(0.001, 0.01, 0.05, 0.1, 0.5),
                         lambda = c(0, 1), competition = c("none", "phylo"),
                         stringsAsFactors = FALSE)

# community statistics of one simulation; the null model draws from all the
# species in the tree, the regional pool
summarise_sim <- function(sim, pool_dist) {
  pa <- (sim$abundance > 0) * 1
  pa <- pa[rowSums(pa) > 1, , drop = FALSE]
  z <- ses.mpd(pa, pool_dist, null.model = "taxa.labels", runs = 99)$mpd.obs.z
  c(alpha = mean(rowSums(pa)), gamma = sum(colSums(pa) > 0),
    ses_mpd = mean(z[is.finite(z)]),
    beta = mean(vegdist(pa[, colSums(pa) > 0], "jaccard", binary = TRUE)))
}

results <- do.call(rbind, lapply(seq_len(n_rep), function(i) {
  # tree, landscape and traits are shared by the scenarios of a replicate
  tree <- ape::rcoal(30)
  pool_dist <- cophenetic(tree)
  land <- landscape_implicit(25)
  traits <- list(`0` = species_traits(tree, lambda = 0, sigma_mean = 8, sigma_sd = 2),
                 `1` = species_traits(tree, lambda = 1, sigma_mean = 8, sigma_sd = 2))
  alpha <- list(none = competition_matrix(tree, "none"),
                phylo = competition_matrix(tree, "phylo", inter_max = 0.5))
  stats <- t(mapply(function(d, l, comp) {
    sim <- simulate_metacommunity(land, traits[[as.character(l)]], alpha[[comp]],
                                  dispersal = d, n_steps = 1000)
    summarise_sim(sim, pool_dist)
  }, scenarios$dispersal, scenarios$lambda, scenarios$competition))
  cbind(scenarios, replicate = i, stats)
}))

saveRDS(results, "vignettes/phylogenetic-structure.rds", version = 2)
