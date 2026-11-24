# Posterior predictive checks under a model that is wrong in a known way: an urban filter
# conserved in the phylogeny, against competition among relatives. Item (ii) of block C of
# dev/ROADMAP.md. In vignette("abc"), the model with phylogenetic competition was chosen with a
# probability near one for the anurans of urban ponds, but its posterior predictive check failed
# on the standardised MPD, and the vignette suggested, as a hypothesis, a filter that favours some
# clades. Here the simulator of data-raw/abc_anurans.R (sourced, with its landscape, pools and
# statistics by city) gains a filter whose thresholds of tolerance to urban cover are conserved:
# a trait drawn on the tree by Brownian motion with Pagel's lambda (lambda_filter, U(0, 1)) ranks
# the species, the share q_tol with the highest values tolerates any cover, and the others get
# thresholds from 0 to 100% in the order of the trait, so that each threshold has the same
# uniform distribution as in the random filter, which lambda_filter = 0 reproduces in law.
# Four models, each with a reference table of N_REF + N_TEST simulations:
#   random_competition     the model of the vignette (random filter, phylogenetic competition);
#   random_none            random filter, no competition (inter_max = 0);
#   conserved_competition  conserved filter and phylogenetic competition;
#   conserved_none         conserved filter, no competition.
# (a) Does the check see the error? The N_TEST pseudo-observations of conserved_none, and as a
#     control those of random_competition, are fitted with the model of the vignette (forests on
#     the statistics without the phylogenetic ones, trained on random_competition) and checked on
#     ses_mpd and ses_mntd of each city with abc_ppc(); the vignette's choice between
#     random_competition and random_none, and the choice among the four models, are predicted for
#     the same pseudo-observations.
# (b) The observed anurans: choice among the four models, posterior of the model chosen, and the
#     posterior predictive check of the phylogenetic statistics under each conserved model.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_conserved_filter.R.
# N_REF, N_TEST, CORES, NTREE and N_PPC change the size; OUTPUT the file.
source("data-raw/abc_anurans.R")
options(comdynamite.cite = FALSE)
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
n_ppc <- as.integer(Sys.getenv("N_PPC", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_conserved_filter.rds")

conserved_pool <- function(pool, urban, q_tol, tree, lambda_filter) {
  k <- ncol(pool)
  z <- trait_bm(tree, lambda_filter)[colnames(pool)]
  rank_z <- rank(z, ties.method = "random")
  tolerant <- rank_z > k - round(q_tol * k)
  # thresholds in the order of the trait, from the uniform quantiles of the random filter
  threshold <- sort(stats::runif(k, 0, 100))[rank_z]
  threshold[tolerant] <- Inf
  filtered <- pool * outer(urban, threshold, "<=")
  lost <- colSums(filtered) == 0
  filtered[, lost] <- pool[, lost]
  filtered
}
simulate_conserved <- function(scenario, shared) {
  landscape <- landscape_points(shared$coords, shared$env, 10^scenario$log_kernel, lonlat = TRUE)
  alpha <- if (is.null(scenario$inter_max)) competition_matrix(shared$tree, "none") else
    competition_matrix(shared$tree, "phylo", inter_max = scenario$inter_max)
  simulate_metacommunity(landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean,
                                        sigma_sd = scenario$sigma_mean / 2),
                         alpha, dispersal = 10^scenario$log_dispersal, n_steps = 500,
                         pool = conserved_pool(shared$pool, shared$urban, scenario$q_tol,
                                               shared$tree, scenario$lambda_filter))
}
no_competition <- prior[names(prior) != "inter_max"]
models <- list(
  random_competition = list(prior = prior, simulate = simulate),
  random_none = list(prior = no_competition, simulate = simulate),
  conserved_competition = list(prior = c(prior, list(lambda_filter = c(0, 1))),
                               simulate = simulate_conserved),
  conserved_none = list(prior = c(no_competition, list(lambda_filter = c(0, 1))),
                        simulate = simulate_conserved))

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
              setup = function() fixed, summarise = summarise, seed = 5000 + m)
  })
  names(tables) <- names(models)
  future::plan(future::sequential)
  minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  started <- Sys.time()
  train <- lapply(tables, keep_rows, rows = seq_len(n_ref))
  test_rows <- n_ref + seq_len(n_test)
  # (a) the model of the vignette against pseudo-observations of each truth
  fit_rest <- abc_fit(train$random_competition, statistics = setdiff(city_stats, phylogenetic),
                      ntree = ntree, ncores = cores)
  ppc_conserved <- check_rows(fit_rest, tables$conserved_none, test_rows,
                              train$random_competition)
  ppc_conserved$lambda_filter <- tables$conserved_none$lambda_filter[ppc_conserved$row]
  ppc_control <- check_rows(fit_rest, tables$random_competition, test_rows,
                            train$random_competition)
  choice2 <- abc_choose(train[c("random_competition", "random_none")], statistics = city_stats,
                        ntree = ntree, ncores = cores)
  choice4 <- abc_choose(train, statistics = city_stats, ntree = ntree, ncores = cores)
  pseudo <- as.data.frame(tables$conserved_none)[test_rows, city_stats]
  chosen2_pseudo <- predict(choice2, pseudo, ntree = ntree, ncores = cores)
  chosen4_pseudo <- predict(choice4, pseudo, ntree = ntree, ncores = cores)
  chosen2_pseudo$lambda_filter <- tables$conserved_none$lambda_filter[test_rows]
  chosen4_pseudo$lambda_filter <- chosen2_pseudo$lambda_filter
  # (b) the observed anurans
  chosen4 <- predict(choice4, observed, ntree = ntree, ncores = cores)
  posterior <- lapply(train[c("conserved_competition", "conserved_none")], function(t) {
    fit <- abc_fit(t, statistics = city_stats, ntree = ntree, ncores = cores)
    rest <- abc_fit(t, statistics = setdiff(city_stats, phylogenetic), ntree = ntree,
                    ncores = cores)
    list(posterior = predict(fit, observed), diagnose = abc_diagnose(fit, observed),
         importance = abc_importance(fit),
         ppc_phylo = abc_ppc(rest, observed, t, statistics = phylogenetic, n = n_ppc,
                             seed = 2043))
  })
  recovery <- abc_recovery(tables$conserved_none, statistics = city_stats, n_test = n_test,
                           ntree = ntree, ncores = cores, seed = 2042)
  minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  saveRDS(list(ppc_conserved = ppc_conserved, ppc_control = ppc_control,
               choice2 = choice2[c("error", "confusion")],
               choice4 = choice4[c("error", "confusion")],
               chosen2_pseudo = chosen2_pseudo, chosen4_pseudo = chosen4_pseudo,
               chosen4 = chosen4, posterior = posterior, recovery = recovery,
               minutes = c(simulation = minutes_sim, forests = minutes_fit),
               settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree, n_ppc = n_ppc),
               date = Sys.Date()),
          output, compress = "xz")
  flag <- function(d) aggregate(p ~ statistic, d, function(p) mean(p < 0.05))
  cat("(a) share of pseudo-observations with p < 0.05\n-- conserved filter, no competition\n")
  print(flag(ppc_conserved), digits = 3)
  cat("-- control: the model of the vignette\n")
  print(flag(ppc_control), digits = 3)
  cat("vignette choice on conserved pseudo-observations:\n")
  print(table(chosen2_pseudo$model))
  cat("four-model choice on conserved pseudo-observations:\n")
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
