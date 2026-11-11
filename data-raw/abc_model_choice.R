# Model choice and posterior predictive checks by ABC with random forests (abc_choose(),
# abc_ppc()), in the design of data-raw/abc_recovery.R: the same tree of 30 species, grid of 25
# sites, priors and statistics. Two models: phylogenetic competition (the reference table of
# abc_recovery.rds, subsampled) and no competition (new simulations, inter_max fixed at 0 and
# the other three parameters from the same priors). Questions set before running:
# (i) how often the two models are told apart (out-of-bag error and confusion matrix), and
# (ii) whether a posterior predictive check on the temporal statistics, left out of the fit on
# the survey statistics, flags pseudo-observations simulated under the wrong model more often
# than those simulated under the right one.
# Run from the root of the package after abc_recovery.R. N_MODEL, N_TEST, CORES and NTREE
# change the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_model <- as.integer(Sys.getenv("N_MODEL", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "100"))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_model_choice.rds")

recovery <- readRDS("data-raw/abc_recovery.rds")
competition <- recovery$table
set.seed(20261107)
fixed <- list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
prior <- attr(competition, "prior")[c("log_dispersal", "sigma_mean", "lambda")]
no_competition <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                         competition_matrix(shared$tree, "none"),
                         dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                         record_steps = seq(800, 1000, by = 50))
}
future::plan(future::multicore, workers = cores)
started <- Sys.time()
neutral <- abc_table(prior, no_competition, n = n_model + n_test, setup = function() fixed,
                     summarise = function(sim, shared) summary_stats(sim, shared$tree, n_null = 49),
                     seed = 2029)
future::plan(future::sequential)
minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

set.seed(2030)
take <- sample.int(nrow(competition), n_model + n_test)
competition <- competition[take, ]
attributes(competition)[c("prior", "parameters", "statistics", "class")] <-
  attributes(recovery$table)[c("prior", "parameters", "statistics", "class")]
is_test <- function(t) seq_len(nrow(t)) > n_model
keep <- function(t, rows) {
  out <- t[rows, ]
  attributes(out)[c("prior", "parameters", "statistics", "class")] <-
    attributes(t)[c("prior", "parameters", "statistics", "class")]
  out
}
temporal <- c("cv_population", "cv_site", "cv_regional", "turnover")
survey <- setdiff(attr(competition, "statistics"), temporal)

# (i) model choice with all statistics, and with the survey only
started <- Sys.time()
train <- list(competition = keep(competition, !is_test(competition)),
              neutral = keep(neutral, !is_test(neutral)))
choice <- lapply(list(series = c(survey, temporal), survey = survey), function(s) {
  abc_choose(train, statistics = s, ntree = ntree, ncores = cores)
})
tests <- list(competition = keep(competition, is_test(competition)),
              neutral = keep(neutral, is_test(neutral)))
held_out <- do.call(rbind, lapply(names(choice), function(view) {
  do.call(rbind, lapply(names(tests), function(truth) {
    p <- predict(choice[[view]], as.data.frame(tests[[truth]]), ntree = ntree, ncores = cores)
    data.frame(view = view, truth = truth, chosen = p$model, probability = p$probability)
  }))
}))

# (ii) fit the competition model on the survey statistics; check the temporal ones
fit <- abc_fit(train$competition, statistics = survey, ntree = ntree, ncores = cores)
ppc <- do.call(rbind, lapply(names(tests), function(truth) {
  t <- as.data.frame(tests[[truth]])
  ok <- which(rowSums(!is.finite(as.matrix(t[c(survey, temporal)]))) == 0)
  do.call(rbind, lapply(ok, function(i) {
    observed <- unlist(t[i, c(survey, temporal)])
    check <- abc_ppc(fit, observed, train$competition, statistics = temporal, n = 500,
                     seed = 2031 + i)
    data.frame(truth = truth, test = i, statistic = check$statistic, p = check$p)
  }))
}))
minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

flagged <- aggregate(p ~ truth + test, ppc, function(p) any(p < 0.05 / length(temporal)))
saveRDS(list(choice = lapply(choice, `[`, c("error", "confusion", "statistics")),
             held_out = held_out, ppc = ppc, flagged = flagged, neutral = neutral,
             minutes = c(simulation = minutes_sim, forests = minutes_fit),
             settings = list(n_model = n_model, n_test = n_test, ntree = ntree)), output)
for (view in names(choice)) {
  cat(view, ": out-of-bag error", round(choice[[view]]$error, 3), "\n")
  print(choice[[view]]$confusion)
}
print(aggregate(cbind(correct = truth == chosen, probability) ~ view + truth, held_out, mean),
      digits = 3)
cat("PPC, fraction of pseudo-observations with a temporal statistic flagged (Bonferroni):\n")
print(aggregate(p ~ truth, flagged, mean), digits = 3)
print(aggregate(p ~ truth + statistic, ppc, function(p) mean(p < 0.05)), digits = 3)
cat("simulation", round(minutes_sim, 1), "min; forests and checks", round(minutes_fit, 1),
    "min\n")
