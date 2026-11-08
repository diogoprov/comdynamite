# T-statistics of Violle et al. (2012) with and without local adaptation (scenario 10 of
# dev/notas-artigo-taxas.md). Individuals are drawn from simulated metacommunities with
# sample_individuals() and the three variance ratios computed with tstats(), in base R for the
# replicates and with cati::Tstats() (null models) for one replicate per scenario.
#
# Expectation (see the section "Local adaptation" in ?tstats): local adaptation pulls the
# populations of a site toward its environment, so co-occurring species become more alike.
# T_PC.PR and T_IC.IR should fall and T_IP.IC should rise as local_adaptation grows.
# Divergence of populations that does not track the environment (population_sd) adds variance
# within sites and in the region alike, so it should change the ratios little.
#
# Run from the root of the package: Rscript data-raw/demo_tstats.R

devtools::load_all(quiet = TRUE)

n_rep <- 20
adaptation <- c(0, 0.25, 0.5, 0.75, 1)
drift_sd <- 5

one_replicate <- function(rep) {
  set.seed(rep)
  tree <- treeshapesim::sim_tree_shape(30)
  sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                                competition_matrix(tree, "none"), dispersal = 0.01,
                                n_steps = 300, seed_steps = 1:10)
  scenarios <- rbind(data.frame(local_adaptation = adaptation, population_sd = 0),
                     data.frame(local_adaptation = 0, population_sd = drift_sd))
  rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    ind <- sample_individuals(sim, n = 10, intraspecific_sd = 2,
                              local_adaptation = scenarios$local_adaptation[i],
                              population_sd = scenarios$population_sd[i])
    t_site <- tstats(ind)
    cbind(rep = rep, scenarios[i, ],
          as.data.frame(t(colMeans(t_site[c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE))))
  })
  do.call(rbind, rows)
}

replicates <- do.call(rbind, lapply(seq_len(n_rep), one_replicate))
summary_table <- aggregate(cbind(T_IP.IC, T_IC.IR, T_PC.PR) ~ local_adaptation + population_sd,
                           replicates, function(v) round(mean(v), 3))
print(summary_table)

# null models of cati for one replicate, without and with strong local adaptation
ses <- NULL
if (requireNamespace("cati", quietly = TRUE)) {
  set.seed(1)
  tree <- treeshapesim::sim_tree_shape(30)
  sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                                competition_matrix(tree, "none"), dispersal = 0.01,
                                n_steps = 300, seed_steps = 1:10)
  ses <- lapply(c(none = 0, strong = 0.8), function(h) {
    res <- tstats(sample_individuals(sim, n = 10, intraspecific_sd = 2, local_adaptation = h),
                  method = "cati", nperm = 99)
    # ses.Tstats() gives T_IP.IC, T_IC.IR and T_PC.PR in this order
    s <- cati::ses.Tstats(res)
    stats::setNames(vapply(s, function(e) mean(e$ses, na.rm = TRUE), numeric(1)),
                    c("T_IP.IC", "T_IC.IR", "T_PC.PR"))
  })
  print(round(do.call(rbind, ses), 2))
}

saveRDS(list(replicates = replicates, summary = summary_table, ses = ses),
        "data-raw/demo_tstats.rds")
