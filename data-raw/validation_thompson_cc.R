# Competition-colonisation scenario of Thompson et al. (2020), to compare with
# their Figure 3g and h. As in their code, the first 30% of species are dominant:
# they disperse ten times less, and each has an effect from U(1, 1.5) on the species
# after it, whereas all other effects are from U(0, 1). Each population is lost with
# probability 0.002 per step. Their code also halves the number of immigrants in
# this scenario only, which the article does not mention and which we reproduce by
# halving the dispersal matrix. Uses the settings of validation_thompson.R. Their
# code removes populations before immigrants arrive and comdynamite after, as
# mcomsimr. Run from the package root; writes data-raw/validation_thompson_cc.rds.
pilot <- TRUE
source("data-raw/validation_thompson.R")

set.seed(2021)
n_dominant <- 0.3 * n_species
cc_breadths <- c(0.5, 10)

# effect[i, j] is the effect of i on j, as in their code; alpha is its transpose
trade_off <- function() {
  effect <- matrix(runif(n_species^2), n_species)
  effect[seq_len(n_dominant), ] <- runif(n_dominant * n_species, 1, 1.5)
  lower <- lower.tri(effect, diag = TRUE)
  effect[lower] <- runif(sum(lower))
  diag(effect) <- 1
  t(effect) * 0.05
}

one_landscape_cc <- function(replicate) {
  xy <- matrix(runif(2 * n_sites, 1, 100), n_sites)
  kernel <- exp(-0.1 * torus_distances(xy))
  diag(kernel) <- 0
  kernel <- kernel / rowSums(kernel)
  env <- environment_field(xy)
  optima <- runif(n_species)
  alpha <- trade_off()
  design <- expand.grid(dispersal = dispersal_rates, niche = cc_breadths)
  out <- t(vapply(seq_len(nrow(design)), function(i) {
    traits <- data.frame(optimum = optima, sigma = sqrt(2) * design$niche[i], r_max = 5)
    rates <- rep(design$dispersal[i], n_species)
    rates[seq_len(n_dominant)] <- rates[seq_len(n_dominant)] / 10
    sim <- mc_dynamics(function(step) growth_rates(env[, step], traits), alpha, rates,
                       n_steps = n_steps, seed_steps = seeds, disp_matrix = 0.5 * kernel,
                       record_steps = record, extirpation = 0.002)
    statistics(sim$history)
  }, numeric(4)))
  cbind(replicate = replicate, design, out)
}

results <- do.call(rbind, parallel::mclapply(1:30, one_landscape_cc, mc.cores = 15))
saveRDS(results, "data-raw/validation_thompson_cc.rds", version = 2)
print(aggregate(cbind(alpha, gamma, beta_spatial, beta_temporal) ~ dispersal + niche,
                results, median), digits = 3)
