# Application of the eco-evolutionary engine to the pond amphibians of Brie (Jeliazkov et al.
# 2014; data from CESTES, see data-raw/jeliazkov2014/prepare.R): 135 ponds at their real
# coordinates, 11 species, and an environmental gradient from woodland to crops and towns (first
# axis of a PCA of the environmental variables). Three models are compared by ABC with random
# forests (abc_choose()):
#   neutral    species equivalent (infinite niche breadth);
#   sorting    Gaussian niches on the gradient, fixed in time;
#   adaptation the same niches, with the optimum of each population evolving (constant G).
# Parameters with uniform priors, shared: log10 emigration probability (-4 to -0.5), log10 mean
# dispersal distance in km (-0.5 to 0.7), ratio of inter- to intraspecific competition (0.2 to
# 1), log10 carrying capacity (2 to 4), log10 detection probability per individual (-3 to 0);
# niche models: niche breadth sigma (5 to 60); adaptation: log10 G (-1 to 1.5) and log10 V_E
# (-0.5 to 1.5). The optima of the 11 species are drawn anew in each run from the range of the
# gradient. Each run lasts 400 steps; observed counts are binomial samples of the abundances at
# the last step with the detection probability, and statistics are computed on the ponds with at
# least one individual counted, as in the data (summary_stats() with gradient and coordinates).
# Before the data are used, recovery is checked on pseudo-observations of each model (criterion
# as in data-raw/abc_recovery.R: R2 >= 0.5 and coverage 85-95%), and the model choice error
# out of bag. Then the observed statistics are classified, the parameters of the chosen model
# estimated and the fit checked with abc_ppc() on statistics held out of the fit.
# Run from the root of the package with abcrf installed. N_REF (per model), N_TEST, CORES and
# NTREE change the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "8000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/application_jeliazkov2014.rds")

dat <- readRDS("data-raw/jeliazkov2014/jeliazkov2014.rds")
species <- colnames(dat$counts)
n_sp <- length(species)
env <- dat$env
coords <- dat$coords

base <- list(log_dispersal = c(-4, -0.5), log_scale = c(-0.5, 0.7), competition = c(0.2, 1),
             log_K = c(2, 4), log_detect = c(-3, 0))
priors <- list(neutral = base,
               sorting = c(base, list(sigma = c(5, 60))),
               adaptation = c(base, list(sigma = c(5, 60), log_G = c(-1, 1.5),
                                         log_VE = c(-0.5, 1.5))))

one_run <- function(scenario, shared, model) {
  land <- landscape_points(coords, env, disp_scale = 10^scenario$log_scale)
  r_max <- 5
  intra <- (r_max - 1) / 10^scenario$log_K
  alpha <- matrix(intra * scenario$competition, n_sp, n_sp, dimnames = list(species, species))
  diag(alpha) <- intra
  sigma <- if (model == "neutral") Inf else scenario$sigma
  traits <- data.frame(optimum = stats::runif(n_sp, 0, 100), sigma = sigma, r_max = r_max,
                       row.names = species)
  evolution <- if (model == "adaptation") {
    evolution_control(10^scenario$log_G, ne_ratio = 0.5, env_var = 10^scenario$log_VE)
  }
  sim <- simulate_metacommunity(land, traits, alpha, dispersal = 10^scenario$log_dispersal,
                                n_steps = 400, record_steps = 400, evolution = evolution)
  counts <- matrix(stats::rbinom(length(sim$abundance), sim$abundance, 10^scenario$log_detect),
                   nrow(sim$abundance), dimnames = dimnames(sim$abundance))
  counts
}
stats_of <- function(counts) {
  keep <- rowSums(counts) > 0
  if (sum(keep) < 10 || sum(colSums(counts) > 0) < 2) return(c(n_ponds = sum(keep)))
  c(summary_stats(counts[keep, colSums(counts) > 0, drop = FALSE], env = env[keep],
                  coords = coords[keep, ], n_null = 19),
    n_ponds = sum(keep), log_total = log10(sum(counts)),
    max_share = max(colSums(counts)) / sum(counts))
}

observed <- stats_of(dat$counts)
print(round(observed, 3))

future::plan(future::multicore, workers = cores)
started <- Sys.time()
tables <- lapply(names(priors), function(m) {
  abc_table(priors[[m]], function(scenario, shared) one_run(scenario, shared, m),
            n = n_ref + n_test, summarise = function(sim, shared) stats_of(sim),
            seed = 3000 + match(m, names(priors)))
})
names(tables) <- names(priors)
future::plan(future::sequential)
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
cat("Simulation:", round(minutes, 1), "min\n")

statistics <- Reduce(intersect, lapply(tables, function(t) attr(t, "statistics")))
statistics <- statistics[statistics %in% names(observed)[is.finite(observed)]]
held_out <- intersect(c("decay", "max_share"), statistics)   # for the posterior predictive check
fit_stats <- setdiff(statistics, held_out)

choice <- abc_choose(tables, statistics = fit_stats, ntree = ntree, ncores = cores)
cat("\nModel choice, out-of-bag error:", round(choice$error, 3), "\n")
print(round(prop.table(choice$confusion, 1), 3))
chosen <- predict(choice, observed[fit_stats])
print(chosen)

recovery <- lapply(c("sorting", "adaptation"), function(m) {
  abc_recovery(tables[[m]], statistics = fit_stats, n_test = n_test, ntree = ntree,
               ncores = cores, seed = 7)
})
names(recovery) <- c("sorting", "adaptation")
for (m in names(recovery)) {
  cat("\nRecovery,", m, "\n")
  print(recovery[[m]], digits = 3, row.names = FALSE)
}

fits <- lapply(c("sorting", "adaptation"), function(m) {
  fit <- abc_fit(tables[[m]], statistics = fit_stats, ntree = ntree, ncores = cores)
  list(fit = fit, estimate = predict(fit, observed[fit_stats], level = 0.9),
       ppc = abc_ppc(fit, observed, tables[[m]], statistics = c(held_out, fit_stats), n = 500))
})
names(fits) <- c("sorting", "adaptation")
for (m in names(fits)) {
  cat("\nEstimates,", m, "\n")
  print(fits[[m]]$estimate, digits = 3, row.names = FALSE)
  cat("Posterior predictive check (held out:", paste(held_out, collapse = ", "), ")\n")
  print(fits[[m]]$ppc, digits = 3)
}

saveRDS(list(observed = observed, tables = tables, choice = choice[c("error", "confusion",
                                                                     "statistics", "models")],
             chosen = chosen, recovery = recovery,
             estimates = lapply(fits, `[[`, "estimate"), ppc = lapply(fits, `[[`, "ppc"),
             priors = priors, held_out = held_out, minutes = minutes, n_ref = n_ref,
             date = Sys.Date()), output, compress = "xz")
