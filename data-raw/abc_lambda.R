# Can the phylogenetic signal of the niche (lambda) be recovered by ABC with random forests from
# statistics that can be observed? Item (i) of block C of dev/ROADMAP.md. In
# data-raw/abc_recovery.R, lambda reached R2 0.55, almost only through Blomberg's K of the
# simulated optima (trait_k), which are not observed in the field; with incidence alone, as in
# the anurans of vignette("abc"), it was not recovered. Here niche_k and niche_dist_cor of
# summary_stats() (Blomberg's K of the realised niche positions of the species present, and the
# Mantel correlation of phylogenetic distance with their difference) are tested in the design of
# data-raw/abc_recovery.R: the same tree of 30 species, grid of 25 sites, priors and simulation,
# one survey at the last step, with 10,000 simulations and 200 pseudo-observations.
# Views of the statistics (all from one survey):
#   community        statistics of abundances and phylogeny (alpha to sad_sd, ses_mpd, ses_mntd,
#                    decay), none of them about the niche;
#   niche            community plus niche_k and niche_dist_cor, from abundances;
#   optima           community plus trait_env, trait_sd and trait_k of the simulated optima, the
#                    statistics of data-raw/abc_recovery.R, as a reference that the field lacks;
#   incidence        alpha, gamma, beta, occupancy, occupancy_sd, ses_mpd, ses_mntd and decay;
#   incidence_niche  incidence plus niche_k and niche_dist_cor from presences.
# Criterion set before running, as in data-raw/abc_recovery.R (revised 2026-11-21): recoverable
# with R2 >= 0.5 and coverage of the 90% interval of at least 85%; `conservative` flags coverage
# above 95%.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_lambda.R. N_REF,
# N_TEST, CORES and NTREE change the size; OUTPUT the file. SPECIES and SIDE change the tree and
# the grid: the sensitivity run with a larger design used SPECIES=100 SIDE=10
# OUTPUT=data-raw/abc_lambda_large.rds.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_lambda.rds")

set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(as.integer(Sys.getenv("SPECIES", "30"))),
              landscape = landscape_grid(side = as.integer(Sys.getenv("SIDE", "5"))))
prior <- list(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
              lambda = c(0, 1))
one_run <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                         competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max),
                         dispersal = 10^scenario$log_dispersal, n_steps = 1000)
}
summarise_run <- function(sim, shared) {
  s <- summary_stats(sim, shared$tree, coords = shared$landscape$coords, n_null = 49)
  s <- s[setdiff(names(s), c("cv_population", "cv_site", "cv_regional", "turnover"))]
  presence <- (sim$abundance > 0) * 1
  incidence <- niche_stats(presence, sim_niche(sim)$env, shared$tree)
  c(s, stats::setNames(incidence, paste0("i_", names(incidence))))
}

community <- c("alpha", "gamma", "beta", "occupancy", "occupancy_sd", "abundance", "evenness",
               "evenness_regional", "sad_sd", "ses_mpd", "ses_mntd", "decay")
incidence <- c("alpha", "gamma", "beta", "occupancy", "occupancy_sd", "ses_mpd", "ses_mntd",
               "decay")
views <- list(community = community,
              niche = c(community, "niche_k", "niche_dist_cor"),
              optima = c(community, "trait_env", "trait_sd", "trait_k"),
              incidence = incidence,
              incidence_niche = c(incidence, "i_niche_k", "i_niche_dist_cor"))

if (sys.nframe() == 0) {
  future::plan(future::multicore, workers = cores)
  started <- Sys.time()
  table <- abc_table(prior, one_run, n = n_ref + n_test, setup = function() fixed,
                     summarise = summarise_run, seed = 2031)
  future::plan(future::sequential)
  minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  results <- do.call(rbind, lapply(names(views), function(view) {
    rec <- abc_recovery(table, statistics = views[[view]], n_test = n_test, ntree = ntree,
                        ncores = cores, seed = 2032)
    cbind(view = view, rec[, c("parameter", "r2", "bias", "coverage", "width")],
          recoverable = rec$r2 >= 0.5 & rec$coverage >= 0.85,
          conservative = rec$coverage > 0.95)
  }))
  fit <- abc_fit(table, statistics = views$niche, ntree = ntree, ncores = cores)
  importance <- abc_importance(fit)
  saveRDS(list(results = results, importance = importance, minutes = minutes_sim,
               views = views, species = length(fixed$tree$tip.label),
               sites = nrow(fixed$landscape$coords), settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree),
               date = Sys.Date()),
          output, compress = "xz")
  print(results, digits = 3)
  print(importance[importance$parameter == "lambda", ][1:6, ], digits = 3)
  cat("simulation", round(minutes_sim, 1), "min\n")
}
