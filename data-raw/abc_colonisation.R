# Posterior predictive checks and model choice under limited colonisation from the region, item
# (iii) of block C of dev/ROADMAP.md. In vignette("abc"), model choice favoured phylogenetic
# competition for the anurans of urban ponds mostly through the low richness of the ponds, and
# the vignette suggested that colonisation from a region where most species are rare could be
# read as competition; a filter conserved in the phylogeny did not remove the misfit of the
# standardised MPD (data-raw/abc_conserved_filter.R). In the simulator of data-raw/abc_anurans.R
# (sourced), every species of the pool gets Poisson(0.5) seeds in every pond at each seeding
# step. Here the seeds of species s have mean 10^log_seed * w_s, with log w_s = sigma_region times
# a standardised trait drawn on the tree by Brownian motion with Pagel's lambda (lambda_region),
# and w scaled to mean 1, so that a few species are common in the region and, with lambda_region
# near 1, related species are similarly common. Priors: log_seed U(-2, 0), sigma_region U(0, 3),
# lambda_region U(0, 1); with sigma_region = 0 and log_seed = log10(0.5), the model of the
# vignette. Four models, each with a reference table of N_REF + N_TEST simulations:
#   random_competition    the model of the vignette (equal seeds, phylogenetic competition);
#   random_none           equal seeds, no competition (inter_max = 0);
#   limited_competition   limited colonisation and phylogenetic competition;
#   limited_none          limited colonisation, no competition.
# (a) Does the check see the error, and is limited colonisation read as competition? The N_TEST
#     pseudo-observations of limited_none, and as a control those of random_competition, are
#     fitted with the model of the vignette (forests on the statistics without the phylogenetic
#     ones, trained on random_competition) and checked on ses_mpd and ses_mntd of each city with
#     abc_ppc(); the vignette's choice between random_competition and random_none, and the choice
#     among the four models, are predicted for the same pseudo-observations.
# (b) The observed anurans: choice among the four models, posterior of each limited model with
#     the posterior predictive check of the phylogenetic statistics, and recovery of the
#     parameters of both limited models.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_colonisation.R.
# N_REF, N_TEST, CORES, NTREE and N_PPC change the size; OUTPUT the file; TABLES, if set, a file
# for the reference tables (large, kept out of git).
source("data-raw/abc_anurans.R")
options(comdynamite.cite = FALSE)
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
n_ppc <- as.integer(Sys.getenv("N_PPC", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_colonisation.rds")

# mean seeds of each species of the pool, from a regional abundance drawn on the tree
regional_seeds <- function(tree, species, log_seed, sigma_region, lambda_region) {
  z <- trait_bm(tree, lambda_region)[species]
  w <- exp(sigma_region * (z - mean(z)) / stats::sd(z))
  10^log_seed * w / mean(w)
}
simulate_limited <- function(scenario, shared) {
  landscape <- landscape_points(shared$coords, shared$env, 10^scenario$log_kernel, lonlat = TRUE)
  alpha <- if (is.null(scenario$inter_max)) competition_matrix(shared$tree, "none") else
    competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max)
  traits <- species_traits(shared$tree, scenario$lambda, sigma_mean = scenario$sigma_mean,
                           sigma_sd = scenario$sigma_mean / 2)
  simulate_metacommunity(landscape, traits, alpha, dispersal = 10^scenario$log_dispersal,
                         n_steps = 500,
                         seed_mean = regional_seeds(shared$tree, rownames(traits),
                                                    scenario$log_seed, scenario$sigma_region,
                                                    scenario$lambda_region),
                         pool = urban_pool(shared$pool, shared$urban, scenario$q_tol))
}
no_competition <- prior[names(prior) != "inter_max"]
limited <- list(log_seed = c(-2, 0), sigma_region = c(0, 3), lambda_region = c(0, 1))
models <- list(
  random_competition = list(prior = prior, simulate = simulate),
  random_none = list(prior = no_competition, simulate = simulate),
  limited_competition = list(prior = c(prior, limited), simulate = simulate_limited),
  limited_none = list(prior = c(no_competition, limited), simulate = simulate_limited))

# rows of a reference table, keeping its attributes
keep_rows <- function(t, rows) {
  out <- t[rows, ]
  attributes(out)[c("prior", "parameters", "statistics", "class")] <-
    attributes(t)[c("prior", "parameters", "statistics", "class")]
  out
}
# posterior predictive p-values of the phylogenetic statistics for each pseudo-observation
check_rows <- function(fit, table, rows, reference) {
  do.call(rbind, lapply(rows, function(i) {
    obs <- unlist(as.data.frame(table)[i, city_stats])
    out <- abc_ppc(fit, obs, reference, statistics = phylogenetic, n = n_ppc, seed = i)
    data.frame(row = i, statistic = out$statistic, observed = out$observed, p = out$p)
  }))
}

if (sys.nframe() == 0) {
  future::plan(future::multicore, workers = cores)
  started <- Sys.time()
  tables <- lapply(seq_along(models), function(m) {
    abc_table(models[[m]]$prior, models[[m]]$simulate, n = n_ref + n_test,
              setup = function() fixed, summarise = summarise, seed = 5100 + m)
  })
  names(tables) <- names(models)
  if (nzchar(Sys.getenv("TABLES"))) saveRDS(tables, Sys.getenv("TABLES"))
  future::plan(future::sequential)
  minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  started <- Sys.time()
  train <- lapply(tables, keep_rows, rows = seq_len(n_ref))
  # pseudo-observations with every statistic (with few colonists, a city may hold too few
  # species for some of them)
  complete <- function(t) {
    rows <- n_ref + seq_len(n_test)
    rows[stats::complete.cases(as.data.frame(t)[rows, city_stats])]
  }
  test_rows <- complete(tables$limited_none)
  control_rows <- complete(tables$random_competition)
  # (a) the model of the vignette against pseudo-observations of each truth
  fit_rest <- abc_fit(train$random_competition, statistics = setdiff(city_stats, phylogenetic),
                      ntree = ntree, ncores = cores)
  ppc_limited <- check_rows(fit_rest, tables$limited_none, test_rows,
                              train$random_competition)
  ppc_limited[c("log_seed", "sigma_region", "lambda_region")] <-
    as.data.frame(tables$limited_none)[ppc_limited$row, c("log_seed", "sigma_region", "lambda_region")]
  ppc_control <- check_rows(fit_rest, tables$random_competition, control_rows,
                            train$random_competition)
  choice2 <- abc_choose(train[c("random_competition", "random_none")], statistics = city_stats,
                        ntree = ntree, ncores = cores)
  choice4 <- abc_choose(train, statistics = city_stats, ntree = ntree, ncores = cores)
  pseudo <- as.data.frame(tables$limited_none)[test_rows, city_stats]
  chosen2_pseudo <- predict(choice2, pseudo, ntree = ntree, ncores = cores)
  chosen4_pseudo <- predict(choice4, pseudo, ntree = ntree, ncores = cores)
  truth <- as.data.frame(tables$limited_none)[test_rows, c("log_seed", "sigma_region",
                                                          "lambda_region", "q_tol")]
  chosen2_pseudo <- cbind(chosen2_pseudo, truth)
  chosen4_pseudo <- cbind(chosen4_pseudo, truth)
  # (b) the observed anurans
  chosen4 <- predict(choice4, observed, ntree = ntree, ncores = cores)
  posterior <- lapply(train[c("limited_competition", "limited_none")], function(t) {
    fit <- abc_fit(t, statistics = city_stats, ntree = ntree, ncores = cores)
    rest <- abc_fit(t, statistics = setdiff(city_stats, phylogenetic), ntree = ntree,
                    ncores = cores)
    list(posterior = predict(fit, observed), diagnose = abc_diagnose(fit, observed),
         importance = abc_importance(fit),
         ppc_phylo = abc_ppc(rest, observed, t, statistics = phylogenetic, n = n_ppc,
                             seed = 2043))
  })
  recovery <- lapply(tables[c("limited_competition", "limited_none")], abc_recovery,
                     statistics = city_stats, n_test = n_test, ntree = ntree, ncores = cores,
                     seed = 2042)
  minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  saveRDS(list(ppc_limited = ppc_limited, ppc_control = ppc_control,
               choice2 = choice2[c("error", "confusion")],
               choice4 = choice4[c("error", "confusion")],
               chosen2_pseudo = chosen2_pseudo, chosen4_pseudo = chosen4_pseudo,
               chosen4 = chosen4, posterior = posterior, recovery = recovery,
               minutes = c(simulation = minutes_sim, forests = minutes_fit),
               settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree, n_ppc = n_ppc),
               complete = c(limited_none = length(test_rows),
                            random_competition = length(control_rows)),
               date = Sys.Date()),
          output, compress = "xz")
  flag <- function(d) aggregate(p ~ statistic, d, function(p) mean(p < 0.05))
  cat("(a) share of pseudo-observations with p < 0.05\n-- limited colonisation, no competition\n")
  print(flag(ppc_limited), digits = 3)
  cat("-- control: the model of the vignette\n")
  print(flag(ppc_control), digits = 3)
  cat("vignette choice on limited pseudo-observations:\n")
  print(table(chosen2_pseudo$model))
  cat("four-model choice on limited pseudo-observations:\n")
  print(table(chosen4_pseudo$model))
  print(choice4$confusion)
  cat("(b) observed\n")
  print(chosen4, digits = 3)
  for (m in names(posterior)) {
    cat("--", m, "\n")
    print(posterior[[m]]$posterior, digits = 3)
    print(posterior[[m]]$ppc_phylo, digits = 3)
  }
  print(recovery, digits = 3)
  cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
}
