# Branch-specific rates of speciation from ClaDS (Maliet et al. 2019) and tip rates DR (Jetz et
# al. 2012) as the edge covariates of fit_damocles() and damocles_samples(), with the example
# fit of ClaDS2 to the Caprimulgidae in RPANDA 2.5 (GPL, only used here; archived on CRAN:
# install the archived TESS and then remotes::install_github("hmorlon/PANDA")).
# 1. dr_rates() against the equal-splits measure of picante::evol.distinct() (DR = 1 / ES).
# 2. Posterior samples of the branch rates of ClaDS from the chains, transformed as in
#    RPANDA::getMAPS_ClaDS() (log rate = l_0 + alpha * alpha_effect + the value of the edge),
#    checked against its MAPs; columns 5 to ncol - 2 of the chains are the edges, in the order of
#    tree$edge.
# 3. Twenty communities simulated under DAMOCLES with colonisation proportional to the square of
#    the MAP branch rates of ClaDS (beta = 2; the rates vary little in this clade), then
#    damocles_samples() over 100 posterior samples of the rates, and fits with the MAP rates and
#    with DR per edge as covariates, compared by AIC with constant rates. The Caprimulgidae tree holds 61% of the species (sample_fraction), so it is
#    not a complete regional pool; the communities here are simulated on it.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(RPANDA)
})
cores <- as.integer(Sys.getenv("CORES", "6"))
RNGkind("L'Ecuyer-CMRG")
set.seed(2019)

data(Caprimulgidae_ClaDS2)
tree <- Caprimulgidae_ClaDS2$tree
sampler <- Caprimulgidae_ClaDS2$sampler

# 1. DR ------------------------------------------------------------------------------------------
es <- picante::evol.distinct(tree, type = "equal.splits")
dr <- dr_rates(tree)
cat("1. max |DR - 1 / ES|:", format(max(abs(dr[es$Species] - 1 / es$w)), digits = 3), "\n")

# 2. ClaDS rates -----------------------------------------------------------------------------
n_par <- ncol(sampler$chains[[1]]) - 2
stopifnot(n_par - 4 == nrow(tree$edge))
kept <- seq(floor(nrow(sampler$chains[[1]]) / 2), nrow(sampler$chains[[1]]))
draws <- do.call(rbind, lapply(sampler$chains, function(chain) chain[kept, , drop = FALSE]))
branch_rates <- t(apply(draws, 1, function(d) {
  exp(d[4] + d[2] * sampler$alpha_effect + d[5:n_par])
}))
maps <- getMAPS_ClaDS(sampler, thin = 1)[-(1:4)]
cat("2. correlation of log posterior medians with log MAPs of the branch rates:",
    round(stats::cor(log(apply(branch_rates, 2, stats::median)), log(maps)), 4), "\n")
samples <- branch_rates[sort(sample.int(nrow(branch_rates), 100)), ]
covariates <- lapply(seq_len(nrow(samples)), function(i) samples[i, ])

# 3. DAMOCLES with ClaDS and DR --------------------------------------------------------------
height <- max(ape::node.depth.edgelength(tree))
z <- log(maps) - mean(log(maps))
comm <- simulate_damocles(tree, exp(2 * z) / height, 1 / height, n = 20)
constant <- fit_damocles(comm, tree, shared = TRUE)
with_map <- fit_damocles(comm, tree, shared = TRUE, edge_covariate = maps, vary = "gamma")
with_dr <- fit_damocles(comm, tree, shared = TRUE, edge_covariate = dr_rates(tree, "edge"),
                        vary = "gamma")
posterior <- damocles_samples(comm, tree, edge_covariate = covariates, vary = "gamma",
                              cores = cores)
cat("3. correlation of log DR per edge with log MAP rates of ClaDS:",
    round(stats::cor(log(dr_rates(tree, "edge")), log(maps)), 3), "\n")
comparison <- rbind(map = damocles_lrt(constant, with_map)[2, ],
                    dr = damocles_lrt(constant, with_dr)[2, ])
comparison$delta_aic <- comparison$AIC - damocles_lrt(constant, with_map)$AIC[1]
comparison$beta_gamma <- c(with_map$beta_gamma[1], with_dr$beta_gamma[1])
print(comparison[, c("loglik", "delta_aic", "statistic", "p_value", "beta_gamma")],
      digits = 4)
cat("   100 posterior samples of the ClaDS rates\n")
print(summary(posterior)$comparison, row.names = FALSE)
print(summary(posterior)$parameters, digits = 3)
saveRDS(list(dr_difference = max(abs(dr[es$Species] - 1 / es$w)), comparison = comparison,
             posterior = posterior$samples), "data-raw/damocles_clads.rds")
