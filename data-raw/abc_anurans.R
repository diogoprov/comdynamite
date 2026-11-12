# ABC with random forests for the anurans of urban ponds (data/anurans.rda), the worked example of
# vignette("abc"). The ponds keep their coordinates (landscape_points()), the environment is the
# hydrology score, and each city has its own pool of species, so that the two cities, about 900 km
# apart, are separate metacommunities with shared parameters. The observed data are survey counts,
# whose totals depend on effort and detection, so the statistics use incidence only.
# Steps: (i) recovery of the parameters of the model with phylogenetic competition from
# pseudo-observations; (ii) choice between that model and one without competition; (iii) the
# posterior given the observed data, with abc_diagnose(); (iv) a posterior predictive check of the
# phylogenetic statistics, left out of a fit on the others; (v) cross-validation with spatial
# blocks, two in each city; (vi) model choice and posterior without total richness and beta
# diversity, since the observed total richness of each city equals its pool by construction.
# Run from the root of the package. N_REF, N_TEST, CORES and NTREE change the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_anurans.rds")

load("data/anurans.rda")
city <- c("Campo Grande" = "cg", "São José dos Campos" = "sjc")[as.character(anurans$city)]
pool <- t(vapply(city, function(k) colSums(anurans$comm[city == k, , drop = FALSE]) > 0,
                 logical(ncol(anurans$comm)))) * 1L
rownames(pool) <- paste0("site", seq_len(nrow(pool)))
fixed <- list(tree = anurans$tree, coords = anurans$coords, env = anurans$env[, "hydrology"],
              pool = pool, city = unname(city))
prior <- list(log_dispersal = c(-4, log10(0.3)), log_kernel = log10(c(0.5, 20)),
              inter_max = c(0, 1), sigma_mean = c(5, 100), lambda = c(0, 1))
simulate <- function(scenario, shared) {
  landscape <- landscape_points(shared$coords, shared$env, 10^scenario$log_kernel, lonlat = TRUE)
  alpha <- if (is.null(scenario$inter_max)) competition_matrix(shared$tree, "none") else
    competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max)
  simulate_metacommunity(landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean,
                                        sigma_sd = scenario$sigma_mean / 2),
                         alpha, dispersal = 10^scenario$log_dispersal, n_steps = 500,
                         pool = shared$pool)
}

# incidence statistics of a set of sites: those of summary_stats() that need no abundances,
# traits or series, plus the correlation of richness with hydrology and that of Jaccard
# dissimilarity with the difference in hydrology among occupied ponds
incidence <- c("alpha", "gamma", "beta", "occupancy", "occupancy_sd", "ses_mpd", "ses_mntd",
               "decay")
incidence_stats <- function(data, tree) {
  s <- summary_stats(data, tree, n_null = 19)[incidence]
  present <- data$comm > 0
  richness <- rowSums(present)
  env <- data$env[[1]]
  occupied <- richness > 0
  link <- function(x, y) {
    if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) NA_real_ else stats::cor(x, y)
  }
  composition <- NA_real_
  if (sum(occupied) >= 3) {
    p <- present[occupied, , drop = FALSE]
    shared_species <- tcrossprod(p)
    union <- outer(rowSums(p), rowSums(p), "+") - shared_species
    jaccard <- 1 - shared_species / union
    lower <- lower.tri(jaccard)
    composition <- link(jaccard[lower], abs(outer(env[occupied], env[occupied], "-"))[lower])
  }
  c(s, richness_env = link(richness, env), composition_env = composition)
}
by_city <- function(data, shared) {
  unlist(lapply(c("cg", "sjc"), function(k) {
    s <- incidence_stats(data[shared$city == k, ], shared$tree)
    stats::setNames(s, paste0(k, "_", names(s)))
  }))
}
# two spatial blocks in each city, for cross-validation; the statistics of each fold pool the
# sites of both cities
folds <- integer(length(city))
for (k in c("cg", "sjc")) {
  i <- city == k
  folds[i] <- spatial_folds(anurans$coords[i, ], k = 2) + if (k == "sjc") 2L else 0L
}
by_fold <- fold_summarise(folds, function(x, shared) incidence_stats(x, shared$tree))
summarise <- function(sim, shared) {
  data <- as_mc_data(sim)
  c(by_city(data, shared), by_fold(data, shared))
}
observed_data <- mc_data(anurans$comm, tree = anurans$tree,
                         env = data.frame(env = fixed$env), coords = anurans$coords)
observed <- summarise(observed_data, fixed)
city_stats <- grep("^(cg|sjc)_", names(observed), value = TRUE)
phylogenetic <- grep("ses_", city_stats, value = TRUE)

future::plan(future::multicore, workers = cores)
started <- Sys.time()
competition <- abc_table(prior, simulate, n = n_ref + n_test, setup = function() fixed,
                         summarise = summarise, seed = 2040)
none <- abc_table(prior[names(prior) != "inter_max"], simulate, n = n_ref,
                  setup = function() fixed, summarise = summarise, seed = 2041)
future::plan(future::sequential)
minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

keep <- function(t, rows) {
  out <- t[rows, ]
  attributes(out)[c("prior", "parameters", "statistics", "class")] <-
    attributes(t)[c("prior", "parameters", "statistics", "class")]
  out
}
reference <- keep(competition, seq_len(n_ref))
started <- Sys.time()
# (i) recovery, with the pseudo-observations left out of the reference table
recovery <- abc_recovery(competition, statistics = city_stats, n_test = n_test, ntree = ntree,
                         ncores = cores, seed = 2042)
# (ii) model choice
choice <- abc_choose(list(phylogenetic = reference, none = none), statistics = city_stats,
                     ntree = ntree, ncores = cores)
chosen <- predict(choice, observed, ntree = ntree, ncores = cores)
# (iii) posterior and diagnostics
fit <- abc_fit(reference, statistics = city_stats, ntree = ntree, ncores = cores)
posterior <- predict(fit, observed)
importance <- abc_importance(fit)
diagnose <- abc_diagnose(fit, observed)
# (iv) posterior predictive check of the phylogenetic statistics
fit_rest <- abc_fit(reference, statistics = setdiff(city_stats, phylogenetic), ntree = ntree,
                    ncores = cores)
ppc_phylo <- abc_ppc(fit_rest, observed, reference, statistics = phylogenetic, n = 500,
                     seed = 2043)
ppc_all <- abc_ppc(fit, observed, reference, statistics = city_stats, n = 500, seed = 2044)
# (v) cross-validation with spatial blocks
cv <- abc_cv(as.data.frame(reference), observed, parameters = names(prior), ntree = ntree,
             ncores = cores, seed = 2045)
# (vi) sensitivity to the statistics bound by the pool
reduced <- grep("_(gamma|beta)$", city_stats, value = TRUE, invert = TRUE)
choice_reduced <- abc_choose(list(phylogenetic = reference, none = none), statistics = reduced,
                             ntree = ntree, ncores = cores)
chosen_reduced <- predict(choice_reduced, observed, ntree = ntree, ncores = cores)
fit_reduced <- abc_fit(reference, statistics = reduced, ntree = ntree, ncores = cores)
posterior_reduced <- predict(fit_reduced, observed)
recovery_reduced <- abc_recovery(competition, statistics = reduced, n_test = n_test,
                                 ntree = ntree, ncores = cores, seed = 2042)
minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

saveRDS(list(observed = observed, folds = folds, prior = prior, recovery = recovery,
             estimates = attr(recovery, "estimates"),
             choice = choice[c("error", "confusion", "statistics", "models")], chosen = chosen,
             posterior = posterior, importance = importance, diagnose = diagnose,
             draws = attr(ppc_all, "draws")[c(names(prior), city_stats)],
             ppc_phylo = ppc_phylo,
             reduced = list(statistics = reduced, recovery = recovery_reduced,
                            choice = choice_reduced[c("error", "confusion")],
                            chosen = chosen_reduced, posterior = posterior_reduced), ppc_all = ppc_all, cv = cv, cv_summary = summary(cv),
             minutes = c(simulation = minutes_sim, forests = minutes_fit),
             settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree)), output)
print(recovery, digits = 3)
cat("model choice: out-of-bag error", round(choice$error, 3), "\n")
print(choice$confusion)
print(chosen, digits = 3)
print(posterior, digits = 3)
print(diagnose, digits = 3)
print(ppc_phylo, digits = 3)
print(as.data.frame(ppc_all)[c("statistic", "observed", "p")], digits = 3)
print(summary(cv), digits = 3)
cat("without gamma and beta: out-of-bag error", round(choice_reduced$error, 3), "\n")
print(chosen_reduced, digits = 3)
print(posterior_reduced, digits = 3)
print(recovery_reduced, digits = 3)
cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
