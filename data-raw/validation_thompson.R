# Reproduces the diversity patterns of Thompson et al. (2020), Figure 3 and the
# hypotheses of their Box 2. Settings follow their "Simulation details": 100
# patches on a torus, 50 species, r_max = 5, 15 dispersal rates, 2200 steps with
# the environment held constant for the first 200, and statistics from step
# 1000 on, every 20 steps. The competition-colonisation scenario is left out
# because comdynamite has no stochastic extirpation.
#
# Differences from the original: RandomFields is no longer on CRAN, so the
# environment is a Gaussian field with exponential spatial covariance evolving
# as an AR(1) process, rescaled to [0, 1]. Their Table 1 was not available to
# us, so the definitions of beta are ours: spatial beta is gamma / alpha at each
# step, and temporal beta is the species seen in a patch over the period divided
# by its mean richness. Run from the package root; writes
# data-raw/validation_thompson.rds.
library(comdynamite)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args)) as.integer(args[1]) else 30
RNGkind("L'Ecuyer-CMRG")
set.seed(2020)

n_sites <- 100
n_species <- 50
n_steps <- 2200
record <- seq(1000, n_steps, 20)
dispersal_rates <- 10^seq(-4, log10(0.464), length.out = 15)
# their sigma, ours is sqrt(2) times this: 0.5 and 10 are the columns of their
# Figure 3, and the narrower values are needed for their hypothesis H4
breadths <- c("0.05" = 0.05, "0.1" = 0.1, "0.2" = 0.2, "0.5" = 0.5, "10" = 10)

torus_distances <- function(xy, size = 100) {
  dx <- abs(outer(xy[, 1], xy[, 1], `-`))
  dy <- abs(outer(xy[, 2], xy[, 2], `-`))
  sqrt(pmin(dx, size - dx)^2 + pmin(dy, size - dy)^2)
}

environment_field <- function(distances, spatial_scale = 50, temporal_scale = 500) {
  # torus distances do not give a positive definite matrix, so negative
  # eigenvalues are set to zero
  e <- eigen(exp(-distances / spatial_scale), symmetric = TRUE)
  root <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
  rho <- exp(-1 / temporal_scale)
  z <- matrix(0, n_sites, n_steps)
  z[, 200] <- root %*% rnorm(n_sites)
  for (t in 201:n_steps) {
    z[, t] <- rho * z[, t - 1] + sqrt(1 - rho^2) * (root %*% rnorm(n_sites))
  }
  z[, 1:199] <- z[, 200]
  (z - min(z)) / (max(z) - min(z))
}

competition <- function(type) {
  a <- switch(type,
              equal = matrix(1, n_species, n_species),
              stabilising = matrix(runif(n_species^2, 0, 0.5), n_species),
              mixed = matrix(runif(n_species^2, 0, 1.5), n_species))
  diag(a) <- 1
  a * 0.05
}

statistics <- function(history) {
  pa <- history > 0
  local <- apply(pa, c(1, 3), sum)
  regional <- apply(apply(pa, c(2, 3), any), 2, sum)
  mean_local <- rowMeans(local)
  seen <- rowSums(apply(pa, c(1, 2), any))
  occupied <- mean_local > 0
  c(alpha = mean(local), gamma = mean(regional),
    beta_spatial = mean((regional / colMeans(local))[colMeans(local) > 0]),
    beta_temporal = mean(seen[occupied] / mean_local[occupied]))
}

one_landscape <- function(replicate) {
  xy <- matrix(runif(2 * n_sites, 1, 100), n_sites)
  distances <- torus_distances(xy)
  kernel <- exp(-0.1 * distances)
  diag(kernel) <- 0
  kernel <- kernel / rowSums(kernel)
  env <- environment_field(distances)
  optima <- runif(n_species)
  alphas <- lapply(c(equal = "equal", stabilising = "stabilising", mixed = "mixed"), competition)
  design <- expand.grid(dispersal = dispersal_rates, niche = names(breadths),
                        competition = names(alphas), stringsAsFactors = FALSE)
  out <- t(vapply(seq_len(nrow(design)), function(i) {
    traits <- data.frame(optimum = optima, sigma = sqrt(2) * breadths[[design$niche[i]]],
                         r_max = 5)
    sim <- mc_dynamics(function(step) growth_rates(env[, step], traits),
                       alphas[[design$competition[i]]], design$dispersal[i],
                       n_steps = n_steps, disp_matrix = kernel, record_steps = record)
    statistics(sim$history)
  }, numeric(4)))
  cbind(replicate = replicate, design, out)
}

if (!exists("pilot")) {
  results <- do.call(rbind, parallel::mclapply(seq_len(n_rep), one_landscape, mc.cores = 15))
  saveRDS(results, "data-raw/validation_thompson.rds", version = 2)
}
