# Can the additive genetic variance G of the eco-evolutionary engine be recovered by ABC with
# random forests (abc_fit(), abc_recovery()) from data that a field study could collect? Eight
# species with ancestral optima spread along an environmental gradient on a 10 x 10 grid, with
# five parameters drawn from uniform priors: log10 G, log10 dispersal, the niche breadth sigma,
# log10 Ne/N, which sets drift, the ratio of inter- to intraspecific competition and log10 of
# the non-genetic phenotypic variance V_E, three nuisance parameters that could be confounded
# with G. V_E is unknown, as in a field study, so the phenotypic variance P = G + V_E of measured
# individuals does not reveal G by itself (a first version fixed V_E = 2, which let the spread of
# measured means within populations carry G directly). Two environments: static for 400 steps, or static for 300 steps and then warming by
# 0.1 per step for 100 steps (+10 units, about two thirds of a mean niche breadth). Every run is
# summarised at its last step in three views of increasing detail:
#   community   abundances of species in sites, with the environment and coordinates
#               (summary_stats() without traits);
#   species     plus the mean trait of each species, the abundance-weighted mean of its
#               populations, as from a trait database;
#   populations plus the mean trait of each population, measured on 10 individuals, from which
#               the slope and R2 of trait on environment within species and the spread of
#               population means are computed (abundance-weighted means over species).
# Criterion set before running, the same as data-raw/abc_recovery.R: a parameter is recoverable
# with R2 >= 0.5 and coverage of the 90% interval between 85% and 95%, in 200
# pseudo-observations left out of the table.
# Run from the root of the package, with abcrf installed. N_REF, N_TEST, CORES and NTREE change
# the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "5000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_genetic_variance.rds")

set.seed(20261118)
side <- 10
land <- landscape_grid(side, env_bandwidth = 3, disp_scale = 1.5)
n_sp <- 8
species <- paste0("sp", seq_len(n_sp))
n_ind <- 10
prior <- list(log_G = c(-1, 1.5), log_dispersal = c(-3, -0.5), sigma = c(8, 25),
              log_ne_ratio = c(-1.3, 0), competition = c(0.2, 1), log_VE = c(-0.5, 1.5))
environments <- list(
  static = matrix(land$env, length(land$env), 400),
  warming = sapply(1:400, function(t) land$env + 0.1 * max(0, t - 300))
)

one_run <- function(scenario, shared) {
  traits <- data.frame(optimum = seq(20, 80, length.out = n_sp), sigma = scenario$sigma,
                       r_max = 5, row.names = species)
  alpha <- matrix(0.02 * scenario$competition, n_sp, n_sp, dimnames = list(species, species))
  diag(alpha) <- 0.02
  simulate_metacommunity(land, traits, alpha, dispersal = 10^scenario$log_dispersal,
                         n_steps = 400, env_time = shared$env_time, record_steps = 400,
                         evolution = evolution_control(10^scenario$log_G, ne_ratio = 10^scenario$log_ne_ratio,
                                                       env_var = 10^scenario$log_VE))
}

# statistics of the three views; traits are measured on phenotypes, with sampling error
summarise_views <- function(sim, shared) {
  comm <- sim$abundance
  env <- sim$evolution$environment
  coords <- land$coords
  community <- summary_stats(comm, env = env, coords = coords, n_null = 19)
  names(community) <- paste0("c_", names(community))
  g <- sim$evolution$settings$genetic_var
  env_var <- sim$evolution$settings$env_var
  pheno <- sim$evolution$phenotype
  measured <- pheno + matrix(stats::rnorm(length(pheno), 0, sqrt((g + env_var) / n_ind)),
                             nrow(pheno))
  measured[comm == 0] <- NA
  total <- colSums(comm)
  sp_mean <- vapply(seq_len(n_sp), function(j) {
    if (total[j] == 0) NA_real_ else sum(measured[, j] * comm[, j], na.rm = TRUE) / total[j]
  }, 0)
  present <- total > 0
  species_view <- summary_stats(comm[, present, drop = FALSE],
                                trait = stats::setNames(sp_mean[present], species[present]),
                                env = env, coords = coords, n_null = 19)
  species_view <- species_view[startsWith(names(species_view), "trait")]
  names(species_view) <- paste0("s_", names(species_view))
  within <- t(vapply(seq_len(n_sp), function(j) {
    ok <- comm[, j] > 0
    if (sum(ok) < 3 || stats::sd(env[ok]) == 0) return(c(NA, NA, NA, NA))
    fit <- stats::lm(measured[ok, j] ~ env[ok])
    c(stats::coef(fit)[2], summary(fit)$r.squared, stats::sd(measured[ok, j]),
      stats::sd(stats::residuals(fit)))
  }, numeric(4)))
  w <- ifelse(is.na(within[, 1]), 0, total)
  wmean <- function(v) if (sum(w) == 0) NA_real_ else sum(v * w, na.rm = TRUE) / sum(w)
  populations <- c(p_slope = wmean(within[, 1]), p_r2 = wmean(within[, 2]),
                   p_sd = wmean(within[, 3]), p_resid_sd = wmean(within[, 4]),
                   p_slope_min = suppressWarnings(min(within[, 1], na.rm = TRUE)),
                   p_slope_max = suppressWarnings(max(within[, 1], na.rm = TRUE)))
  populations[!is.finite(populations)] <- NA
  c(community, species_view, populations)
}

future::plan(future::multicore, workers = cores)
results <- list()
for (e in names(environments)) {
  started <- Sys.time()
  shared <- list(env_time = environments[[e]])
  table <- abc_table(prior, one_run, n = n_ref + n_test, setup = function() shared,
                     summarise = summarise_views, seed = 2026 + match(e, names(environments)))
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  stats <- attr(table, "statistics")
  views <- list(community = stats[startsWith(stats, "c_")],
                species = stats[startsWith(stats, "c_") | startsWith(stats, "s_")],
                populations = stats)
  recovery <- lapply(views, function(v) {
    abc_recovery(table, statistics = v, n_test = n_test, ntree = ntree, ncores = cores,
                 seed = 7)
  })
  results[[e]] <- list(table = table, recovery = recovery, minutes = minutes)
  cat("\n==", e, "(", round(minutes, 1), "min of simulation )\n")
  for (v in names(recovery)) {
    cat("--", v, "\n")
    print(recovery[[v]], digits = 3, row.names = FALSE)
  }
}
future::plan(future::sequential)

summary_table <- do.call(rbind, lapply(names(results), function(e) {
  do.call(rbind, lapply(names(results[[e]]$recovery), function(v) {
    data.frame(environment = e, view = v, results[[e]]$recovery[[v]][, c(
      "parameter", "r2", "bias", "coverage", "width")])
  }))
}))
summary_table$recoverable <- summary_table$r2 >= 0.5 & summary_table$coverage >= 0.85 &
  summary_table$coverage <= 0.95
print(summary_table, digits = 3, row.names = FALSE)
saveRDS(list(summary = summary_table,
             recovery = lapply(results, `[[`, "recovery"),
             tables = lapply(results, `[[`, "table"),
             minutes = vapply(results, `[[`, 0, "minutes"),
             prior = prior, n_ref = n_ref, n_test = n_test, date = Sys.Date()),
        output, compress = "xz")
