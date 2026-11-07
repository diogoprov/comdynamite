# Simulations of vignettes/evaluating-methods.Rmd. Run from the package root;
# writes vignettes/evaluating-methods.rds, which the vignette reads so that it
# builds quickly. Needs vegan.
library(comdynamite)
library(vegan)
set.seed(2026)

n_rep <- 200
n_perm <- 199
side <- 10
torus <- torus_permutations(side)

# the method under evaluation: proportion of the variation in Hellinger-
# transformed abundances explained by one predictor, as in a redundancy analysis
r2_env <- function(y, x) {
  y <- scale(y, scale = FALSE)
  x <- x - mean(x)
  sum((crossprod(x, y) / sum(x^2))^2) * sum(x^2) / sum(y^2)
}
p_value <- function(y, x, perms) {
  obs <- r2_env(y, x)
  null <- apply(perms, 1, function(p) r2_env(y, x[p]))
  (sum(null >= obs) + 1) / (length(null) + 1)
}

one_replicate <- function(i) {
  tree <- treeshapesim::sim_tree_shape(30)
  land <- landscape_grid(side, env_bandwidth = 2, disp_scale = 2)
  # autocorrelated like the environment, but unrelated to the community
  unrelated <- landscape_grid(side, env_bandwidth = 2)$env
  alpha <- competition_matrix(tree, "none")
  free <- t(replicate(n_perm, sample(side^2)))
  restricted <- torus[sample(nrow(torus), n_perm), ]
  do.call(rbind, lapply(c(0, 1), function(lambda) {
    traits <- species_traits(tree, lambda = lambda, sigma_mean = 8, sigma_sd = 2)
    sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.05, n_steps = 500)
    y <- decostand(sim$abundance[, colSums(sim$abundance) > 0], "hellinger")
    data.frame(replicate = i, lambda = lambda,
               predictor = rep(c("unrelated", "environment"), each = 2),
               null_model = rep(c("free", "torus"), 2),
               p = c(p_value(y, unrelated, free), p_value(y, unrelated, restricted),
                     p_value(y, land$env, free), p_value(y, land$env, restricted)))
  }))
}

results <- do.call(rbind, lapply(seq_len(n_rep), one_replicate))
saveRDS(results, "vignettes/evaluating-methods.rds", version = 2)
