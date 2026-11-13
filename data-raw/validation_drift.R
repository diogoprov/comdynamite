# Validation of drift_test() on simulated clades. Run from the package root; writes
# data-raw/validation_drift.rds.
#
# A clade evolves in simulate_evolving_metacommunity() with eight traits, speciation by fission
# and a known G whose eigenvalues span two orders of magnitude. Under "neutral", every niche
# axis has sigma = Inf, so the species diverge by drift alone and B is expected to be
# proportional to G (slope 1). Under "selection", the environment varies among sites along the
# eigenvector of G with the smallest eigenvalue and the niche is finite on every axis, so
# divergent selection spreads the species along the minor axis of G (slope below 1 with W = G).
# Each run applies drift_test() (with the sampling correction) to the species means with three
# choices of W: the known G; "population", the pooled phenotypic covariance within populations
# (species x site), which estimates P = G + E with E not proportional to G; and "species", the
# pooled covariance within species, which also holds the divergence among the populations of a
# species (local adaptation under "selection").
#
# Results (60 runs per scenario; runs with fewer than three species, or in which divergence on
# some component was below the sampling error of the means, were not tested):
# - neutral, W = G: median slope 0.94 (IQR 0.85-1.07), 3.6% rejected at 5%: the test holds its
#   size when G is known.
# - selection, W = G: median slope -0.22, every slope below 1, 36% rejected: the direction is
#   always right, but with eight traits (6 df) the test has little power.
# - neutral, W = P within populations: median slope 1.94, 57% rejected. With E = I and the
#   eigenvalues of G from 8 to 0.08, heritability is low on the minor axes, so P is far from
#   proportional to G and the test rejects drift that is true (as Proa et al. 2013 found for
#   dissimilar G and P).
# - selection, W = P within populations: every slope below 1, 10% rejected.
# - W pooled within species: slope 1.27 under drift (20% rejected) and 1.15 under selection
#   (none rejected): the divergence among the populations of a species enters W, and under
#   divergent selection it hides the selection among species. W must come from within
#   populations, not from species that span the gradient.
library(comdynamite)

n_rep <- 60
q <- 8
axes <- paste0("t", seq_len(q))
land <- landscape_grid(6, disp_scale = 0.5)
set.seed(20261112)
rotation <- qr.Q(qr(matrix(rnorm(q * q), q)))
G <- rotation %*% diag(2^seq(3, -3.6, length.out = q)) %*% t(rotation)
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
  # the pooled covariance within populations (species x site)
  pop <- interaction(ind$species, ind$site, drop = TRUE)
  x <- as.matrix(ind[axes])
  centred <- x - (rowsum(x, pop) / as.vector(table(pop)))[as.integer(pop), ]
  P <- crossprod(centred) / (nrow(x) - nlevels(pop))
  one <- function(label, W) {
    t <- tryCatch(drift_test(ind, "species", traits = axes, W = W), error = function(e) NULL)
    if (is.null(t)) {
      return(data.frame(scenario = scenario, seed = seed, W = label, slope = NA_real_,
                        p_value = NA_real_, n_groups = length(unique(ind$species))))
    }
    data.frame(scenario = scenario, seed = seed, W = label, slope = t$slope,
               p_value = t$p_value, n_groups = t$n_groups)
  }
  rbind(one("G", G), one("population", P), one("species", NULL))
}

runs <- do.call(rbind, lapply(names(scenarios), function(sc) {
  do.call(rbind, lapply(seq_len(n_rep), function(i) run_one(sc, i)))
}))
summary_table <- do.call(rbind, lapply(split(runs, list(runs$scenario, runs$W)), function(d) {
  ok <- !is.na(d$slope)
  data.frame(scenario = d$scenario[1], W = d$W[1], runs = nrow(d), tested = sum(ok),
             median_species = stats::median(d$n_groups),
             median_slope = stats::median(d$slope[ok]),
             slope_q25 = unname(stats::quantile(d$slope[ok], 0.25)),
             slope_q75 = unname(stats::quantile(d$slope[ok], 0.75)),
             below_1 = mean(d$slope[ok] < 1), rejected = mean(d$p_value[ok] < 0.05))
}))
rownames(summary_table) <- NULL
print(summary_table)
saveRDS(list(G = G, runs = runs, summary = summary_table), "data-raw/validation_drift.rds")
