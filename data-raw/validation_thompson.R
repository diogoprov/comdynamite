# Reproduces the diversity patterns of Thompson et al. (2020), Figures 3 and 4 and
# the hypotheses of their Box 2: 100 patches on a torus, 50 species, r_max = 5, 15
# dispersal rates, 2200 steps with the environment held constant for the first 200,
# and statistics from step 1000 on, every 20 steps. The competition-colonisation
# scenario is in validation_thompson_cc.R.
#
# Where the article and the authors' code (Meta_com_framework v1.0, Zenodo
# 10.5281/zenodo.3708350) differ, we follow the code, which produced the figures:
# - beta is gamma / alpha, not the gamma - alpha of their Table 1;
# - alpha, gamma and temporal gamma leave out empty patches and time points;
# - the environment is converted to ranks over all patches and steps, so it is
#   uniform on [0, 1], which the article does not mention;
# - dispersal starts at 0.00001, not at 0.0001;
# - sites are seeded every 20 steps up to step 100, not every 10.
#
# Differences from the original: RandomFields is no longer on CRAN, so the
# environment before the conversion to ranks is a Gaussian field with exponential
# spatial covariance evolving as an AR(1) process, whereas theirs has one covariance
# over space and time. Their code draws emigrants from a Poisson distribution and
# the article from a binomial, which we use. We run 6 of their 13 niche breadths.
# Run from the package root; writes data-raw/validation_thompson.rds.
library(comdynamite)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args)) as.integer(args[1]) else 30
RNGkind("L'Ecuyer-CMRG")
set.seed(2020)

n_sites <- 100
n_species <- 50
n_steps <- 2200
record <- seq(1000, n_steps, 20)
seeds <- seq(20, 100, 20)
dispersal_rates <- exp(seq(log(1e-5), 0, length.out = 16))[1:15]
# their sigma, ours is sqrt(2) times this: 0.5 and 10 are the columns of their
# Figure 3, and the others are from the grid of their Figure 4
breadths <- c("0.0464" = 0.0464, "0.1" = 0.1, "0.215" = 0.215, "0.5" = 0.5, "1" = 1, "10" = 10)

torus_distances <- function(xy, size = 100) {
  dx <- abs(outer(xy[, 1], xy[, 1], `-`))
  dy <- abs(outer(xy[, 2], xy[, 2], `-`))
  sqrt(pmin(dx, size - dx)^2 + pmin(dy, size - dy)^2)
}

environment_field <- function(xy, spatial_scale = 50, temporal_scale = 500) {
  # as in their code, the field itself is not on a torus
  root <- t(chol(exp(-as.matrix(dist(xy)) / spatial_scale)))
  rho <- exp(-1 / temporal_scale)
  repeat {
    z <- matrix(0, n_sites, n_steps)
    z[, 200] <- root %*% rnorm(n_sites)
    for (t in 201:n_steps) {
      z[, t] <- rho * z[, t - 1] + sqrt(1 - rho^2) * (root %*% rnorm(n_sites))
    }
    z[, 1:199] <- z[, 200]
    ranks <- matrix(ecdf(z[, 200:n_steps])(z), n_sites)
    # their condition, applied before the conversion to ranks
    first <- (z[, 200] - min(z)) / (max(z) - min(z))
    if (max(first) - min(first) > 0.6) return(ranks)
  }
}

competition <- function(type) {
  a <- switch(type,
              equal = matrix(1, n_species, n_species),
              stabilising = matrix(runif(n_species^2, 0, 0.5), n_species),
              mixed = matrix(runif(n_species^2, 0, 1.5), n_species))
  diag(a) <- 1
  a * 0.05
}

# as in their process_julia_outputs.R, empty patches and time points do not count
statistics <- function(history) {
  pa <- history > 0
  drop_zero <- function(x) replace(x, x == 0, NA)
  local <- drop_zero(apply(pa, c(1, 3), sum))
  regional <- drop_zero(apply(apply(pa, c(2, 3), any), 2, sum))
  seen <- drop_zero(rowSums(apply(pa, c(1, 2), any)))
  alpha <- mean(rowMeans(local, na.rm = TRUE), na.rm = TRUE)
  gamma <- mean(regional, na.rm = TRUE)
  gamma_time <- mean(seen, na.rm = TRUE)
  c(alpha = alpha, gamma = gamma, beta_spatial = gamma / alpha,
    beta_temporal = gamma_time / alpha)
}

one_landscape <- function(replicate) {
  xy <- matrix(runif(2 * n_sites, 1, 100), n_sites)
  kernel <- exp(-0.1 * torus_distances(xy))
  diag(kernel) <- 0
  kernel <- kernel / rowSums(kernel)
  env <- environment_field(xy)
  optima <- runif(n_species)
  alphas <- lapply(c(equal = "equal", stabilising = "stabilising", mixed = "mixed"), competition)
  design <- expand.grid(dispersal = dispersal_rates, niche = names(breadths),
                        competition = names(alphas), stringsAsFactors = FALSE)
  out <- t(vapply(seq_len(nrow(design)), function(i) {
    traits <- data.frame(optimum = optima, sigma = sqrt(2) * breadths[[design$niche[i]]],
                         r_max = 5)
    sim <- mc_dynamics(function(step) growth_rates(env[, step], traits),
                       alphas[[design$competition[i]]], design$dispersal[i],
                       n_steps = n_steps, seed_steps = seeds, disp_matrix = kernel,
                       record_steps = record)
    statistics(sim$history)
  }, numeric(4)))
  cbind(replicate = replicate, design, out)
}

if (!exists("pilot")) {
  results <- do.call(rbind, parallel::mclapply(seq_len(n_rep), one_landscape, mc.cores = 15))
  saveRDS(results, "data-raw/validation_thompson.rds", version = 2)
}
