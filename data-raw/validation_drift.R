# Validation of drift_test() on simulated clades. Run from the package root; writes
# data-raw/validation_drift.rds.
#
# A clade evolves in simulate_evolving_metacommunity() with six traits, speciation by fission
# and a known G whose eigenvalues span an order of magnitude. Under "neutral", every niche axis
# has sigma = Inf, so the species diverge by drift alone and B is expected to be proportional
# to G (slope 1). Under "selection", the environment varies among sites along the eigenvector
# of G with the smallest eigenvalue and the niche is finite on every axis, so divergent
# selection spreads the species along the minor axis of G (slope below 1). Each run applies
# drift_test() to the species means with W = G (known) and with W = the pooled phenotypic
# covariance within species, which adds the environmental variance E (not proportional to G)
# and the variance among sites within species.
library(comdynamite)

n_rep <- 60
q <- 6
axes <- paste0("t", seq_len(q))
land <- landscape_grid(6, disp_scale = 0.5)
set.seed(20261112)
rotation <- qr.Q(qr(matrix(rnorm(q * q), q)))
G <- rotation %*% diag(2^(3:-2)) %*% t(rotation)
dimnames(G) <- list(axes, axes)
minor <- rotation[, q]
gradient <- as.vector(scale(land$env))

traits_for <- function(sigma) {
  stats::setNames(lapply(axes, function(a) {
    data.frame(optimum = 0, sigma = sigma, r_max = 3, row.names = "sp1")
  }), axes)
}
scenarios <- list(
  neutral = list(traits = traits_for(Inf), env = matrix(0, length(gradient), q)),
  selection = list(traits = traits_for(4), env = outer(gradient, minor) * 6)
)

run_one <- function(scenario, seed) {
  set.seed(seed)
  s <- scenarios[[scenario]]
  env <- s$env
  colnames(env) <- axes
  sim <- simulate_evolving_metacommunity(land, s$traits, speciation_control("fission", rate = 0.01),
                                         evolution = evolution_control(G, env_var = diag(q),
                                                                       ne_ratio = 0.05),
                                         env = env, dispersal = 0.002, n_steps = 400)
  ind <- sample_individuals(sim, n = 5, source = "evolution")
  if (length(unique(ind$species)) < 3) return(NULL)
  rbind(
    data.frame(scenario = scenario, seed = seed, W = "G",
               unclass(drift_test(ind, "species", traits = axes, W = G))[c("slope", "p_value",
                                                                          "n_groups")]),
    data.frame(scenario = scenario, seed = seed, W = "pooled",
               unclass(drift_test(ind, "species", traits = axes))[c("slope", "p_value",
                                                                   "n_groups")])
  )
}

runs <- do.call(rbind, lapply(names(scenarios), function(sc) {
  do.call(rbind, lapply(seq_len(n_rep), function(i) run_one(sc, i)))
}))
summary_table <- do.call(rbind, lapply(split(runs, list(runs$scenario, runs$W)), function(d) {
  data.frame(scenario = d$scenario[1], W = d$W[1], runs = nrow(d),
             median_species = stats::median(d$n_groups), median_slope = stats::median(d$slope),
             rejected = mean(d$p_value < 0.05))
}))
rownames(summary_table) <- NULL
print(summary_table)
saveRDS(list(G = G, runs = runs, summary = summary_table), "data-raw/validation_drift.rds")
