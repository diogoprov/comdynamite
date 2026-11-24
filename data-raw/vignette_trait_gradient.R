# Simulated individuals and Bayesian fits for the article vignettes/articles/trait-gradient.Rmd.
# The same seed gives a different metacommunity on different platforms (the site is built on
# Linux, and macOS gave other species and values), so the article reads the individuals from
# here instead of simulating them, and every reader gets the same numbers. The fits of
# trait_gradient_bayes() need cmdstanr, which the site build does not have, so their output is
# saved here too.
#
#   Rscript data-raw/vignette_trait_gradient.R
# Writes vignettes/articles/trait-gradient.rds and vignettes/articles/trait_gradient_ppc.png.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)

set.seed(1)
tree <- treeshapesim::sim_tree_shape(30)
sim <- simulate_metacommunity(landscape_implicit(25), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)

# as in the article
tg <- trait_gradient(ind, abundance = "abundance")
set.seed(2)
ph <- trait_gradient_phylo(tg, tree, n_boot = 199)
tg_occ <- trait_gradient(ind, abundance = "abundance", species_weights = "occurrence")
bayes <- trait_gradient_bayes(tg_occ, tree, parallel_chains = 4, seed = 3)
a <- ph$models$a[ph$models$response == "beta" & ph$models$model == "OU"]
ou <- trait_gradient_bayes(tg_occ, tree, ou_alpha = a, parallel_chains = 4, seed = 4)
comparison <- loo::loo_compare(bayes$fit$loo(variables = "log_lik_p"),
                               ou$fit$loo(variables = "log_lik_p"))
ppc <- bayesplot::pp_check(bayes, resp = "sitemean")
ggplot2::ggsave("vignettes/articles/trait_gradient_ppc.png", ppc, width = 6, height = 4,
                dpi = 100, bg = "white")

saveRDS(list(tree = tree, individuals = ind, ou_alpha = a,
             bayes = list(diagnostics = bayes$diagnostics, signal = bayes$signal),
             ou = list(diagnostics = ou$diagnostics, signal = ou$signal),
             loo = as.data.frame(comparison[, c("elpd_diff", "se_diff")]),
             session = list(platform = R.version$platform,
                            cmdstan = cmdstanr::cmdstan_version(),
                            date = Sys.Date())),
        "vignettes/articles/trait-gradient.rds")
