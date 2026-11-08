# T-statistics of Violle et al. (2012) under local adaptation and character displacement
# (scenario 10 of dev/notas-artigo-taxas.md). Individuals are drawn from simulated
# metacommunities with sample_individuals() and the three variance ratios computed with
# tstats(): observed values for 20 replicates, and standardised effect sizes from the null
# models in base R for one replicate.
#
# Expectations (see the section "Local adaptation and character displacement" in ?tstats):
# local adaptation pulls the populations of a site toward its environment, so co-occurring
# species become more alike: T_PC.PR and T_IC.IR fall and T_IP.IC rises. Character
# displacement pushes co-occurring populations apart: T_IP.IC falls and T_PC.PR rises.
# Divergence of populations that does not track the environment (population_sd) adds variance
# within sites and in the region alike, so it should change the ratios little.
#
# Run from the root of the package: Rscript data-raw/demo_tstats.R

devtools::load_all(quiet = TRUE)

n_rep <- 20
ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
scenarios <- rbind(
  expand.grid(local_adaptation = c(0, 0.25, 0.5, 0.75, 0.9), displacement = c(0, 2, 10),
              population_sd = 0),
  data.frame(local_adaptation = 0, displacement = 0, population_sd = 5)
)

simulate_one <- function(rep) {
  set.seed(rep)
  tree <- treeshapesim::sim_tree_shape(30)
  simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                         competition_matrix(tree, "none"), dispersal = 0.01,
                         n_steps = 300, seed_steps = 1:10)
}

draw <- function(sim, i) {
  sample_individuals(sim, n = 10, intraspecific_sd = 2,
                     local_adaptation = scenarios$local_adaptation[i],
                     displacement = scenarios$displacement[i],
                     population_sd = scenarios$population_sd[i])
}

one_replicate <- function(rep) {
  sim <- simulate_one(rep)
  rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    t_site <- tstats(draw(sim, i), nperm = 0)
    cbind(rep = rep, scenarios[i, ], as.data.frame(t(colMeans(t_site[ratios], na.rm = TRUE))))
  })
  do.call(rbind, rows)
}

replicates <- do.call(rbind, lapply(seq_len(n_rep), one_replicate))
summary_table <- aggregate(cbind(T_IP.IC, T_IC.IR, T_PC.PR) ~
                             local_adaptation + displacement + population_sd,
                           replicates, function(v) round(mean(v), 3))
print(summary_table)

# null models in base R for one replicate: mean standardised effect size over sites and the
# share of sites below the null at 5% (alternative = "less")
sim <- simulate_one(1)
set.seed(1)
null_scenarios <- c(1, 5, 11, 15)
ses <- do.call(rbind, lapply(null_scenarios, function(i) {
  res <- tstats(draw(sim, i), nperm = 199)
  cbind(scenarios[i, ],
        as.data.frame(t(colMeans(res[paste0("ses_", ratios)], na.rm = TRUE))),
        as.data.frame(t(colMeans(res[paste0("p_", ratios)] < 0.05, na.rm = TRUE))))
}))
print(ses, digits = 2)

saveRDS(list(replicates = replicates, summary = summary_table, ses = ses),
        "data-raw/demo_tstats.rds")
