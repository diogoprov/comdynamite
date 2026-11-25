# Does a temporal statistic of drift between consecutive generations make the effective size
# identifiable by ABC with random forests? Item (x) of block C of dev/ROADMAP.md. In
# data-raw/abc_identifiability.R (item viii) the two surveys were 50 steps apart, so selection
# and gene flow moved the means as much as drift did, and Ne/N gained little (R2 at most +0.10).
# In data-raw/abc_g_ne.R (item ix), abundances and trait means informed G / (Ne/N) well (R2 0.92
# to 0.96) but not Ne/N (0.29 to 0.61). Here population_stats() gains temporal_step, the mean
# squared change of population means between consecutive surveys, corrected for sampling and
# divided by the interval and the variance within populations, and temporal_step_slope, its
# slope on 1 / abundance.
# Simulator, statistics of the community, priors and landscape of data-raw/abc_identifiability.R
# (sourced through data-raw/abc_g_ne.R, which adds the derived parameters), with unknown ancestral optima, a static or warming environment and G constant or set
# by mutation (four cases). Sampling designs, all with the complete community survey:
#   base           10 individuals per population at every occupied site, at step 400;
#   surveys2       the same at steps 350 and 400 (the design of item viii);
#   consecutive10  10 individuals at each of the steps 391 to 400;
#   consecutive30  30 individuals at each of the steps 391 to 400.
# Prediction written before running: drift per generation is also proportional to G / Ne, so the
# temporal data should again inform G / Ne, and Ne/N only as far as G is known; G / Ne should get
# more precise and Ne/N gain little.
# Criterion set before running: Ne/N recoverable (R2 >= 0.5 and coverage of the 90% interval of
# at least 85%, in 200 pseudo-observations) and a gain of at least 0.10 in its R2 over surveys2.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_temporal_drift.R.
# N_REF, N_TEST, CORES and NTREE change the size; OUTPUT the file; TABLES_DIR, if set, keeps the
# reference tables.
source("data-raw/abc_g_ne.R")  # also sources data-raw/abc_identifiability.R
output <- Sys.getenv("OUTPUT", "data-raw/abc_temporal_drift.rds")
tables_dir <- Sys.getenv("TABLES_DIR")
steps <- c(350, 391:400)

run_series <- function(scenario, shared) {
  centre <- scenario$opt_centre %||% 50
  spread <- scenario$opt_spread %||% 60
  traits <- data.frame(optimum = centre + spread * seq(-0.5, 0.5, length.out = n_sp),
                       sigma = scenario$sigma, r_max = 5, row.names = species)
  alpha <- matrix(0.02 * scenario$competition, n_sp, n_sp, dimnames = list(species, species))
  diag(alpha) <- 0.02
  evolution <- if (is.null(scenario$log_Vm)) {
    evolution_control(10^scenario$log_G, ne_ratio = 10^scenario$log_ne_ratio,
                      env_var = 10^scenario$log_VE)
  } else {
    evolution_control(1, ne_ratio = 10^scenario$log_ne_ratio, env_var = 10^scenario$log_VE,
                      mutation_var = 10^scenario$log_Vm)
  }
  simulate_metacommunity(land, traits, alpha, dispersal = 10^scenario$log_dispersal,
                         n_steps = 400, env_time = shared$env_time, record_steps = steps,
                         evolution = evolution)
}

# individuals measured in every population at every recorded step (sites x species x step x
# individual); with mutation, G of the last step is used at every step, a small approximation
measure_series <- function(sim) {
  ev <- sim$evolution
  g <- ev$genetic_variance
  if (is.null(g)) g <- matrix(ev$settings$genetic_var, side^2, n_sp)
  sd_ij <- sqrt(g + ev$settings$env_var)
  sd_ij[!is.finite(sd_ij)] <- 0
  abund <- sim$history[, , paste0("step", steps), drop = FALSE]
  out <- array(NA_real_, c(side^2, n_sp, length(steps), max_ind))
  for (k in seq_along(steps)) {
    means <- ev$history[, , k]
    for (r in seq_len(max_ind)) {
      x <- means + stats::rnorm(length(means), 0, sd_ij)
      x[abund[, , k] < r] <- NA
      out[, , k, r] <- x
    }
  }
  list(ind = out, abund = abund)
}

series_designs <- list(base = list(n = 10, steps = 400),
                       surveys2 = list(n = 10, steps = c(350, 400)),
                       consecutive10 = list(n = 10, steps = 391:400),
                       consecutive30 = list(n = 30, steps = 391:400))

series_stats <- function(m, env, n_ind, at) {
  k <- match(at, steps)
  x <- m$ind[, , k, seq_len(n_ind), drop = FALSE]
  idx <- arrayInd(seq_along(x), dim(x))
  ok <- is.finite(x)
  ind <- data.frame(step = at[idx[ok, 3]], site = site_names[idx[ok, 1]],
                    species = species[idx[ok, 2]],
                    abundance = m$abund[, , k, drop = FALSE][idx[ok, 1:3, drop = FALSE]],
                    env = env[idx[ok, 1]], trait = x[ok])
  out <- population_stats(ind, abundance = "abundance", step = if (length(at) > 1) "step",
                          coords = site_coords)
  out[!is.finite(out)] <- NA
  out
}

summarise_series <- function(sim, shared) {
  env <- sim$evolution$environment
  community <- summary_stats(sim$abundance, env = env, coords = land$coords, n_null = 19)
  names(community) <- paste0("c_", names(community))
  m <- measure_series(sim)
  pops <- unlist(lapply(names(series_designs), function(nm) {
    d <- series_designs[[nm]]
    stats::setNames(series_stats(m, env, d$n, d$steps),
                    paste0(nm, "_", population_stat_names_all))
  }))
  c(community, pops)
}
population_stat_names_all <- c(population_stat_names, "slope_min", "slope_max")

series_views <- function(stats) {
  community <- stats[startsWith(stats, "c_")]
  out <- list()
  for (nm in names(series_designs)) {
    s <- stats[startsWith(stats, paste0(nm, "_"))]
    out[[nm]] <- c(community, s)
  }
  out
}
series_cases <- cases[cases$optima == "unknown", ]

if (sys.nframe() == 0) {
  future::plan(future::multicore, workers = cores)
  results <- list()
  for (i in seq_len(nrow(series_cases))) {
    cs <- series_cases[i, ]
    key <- paste(cs, collapse = "_")
    started <- Sys.time()
    shared <- list(env_time = environments[[cs$environment]])
    prior <- make_prior(cs$optima, cs$genetics)
    table <- abc_table(prior, run_series, n = n_ref + n_test, setup = function() shared,
                       summarise = summarise_series, seed = 8000 + i)
    minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    if (nzchar(tables_dir)) {
      dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(table, file.path(tables_dir, paste0(key, ".rds")), compress = "xz")
    }
    views <- series_views(attr(table, "statistics"))
    # statistics that are NA in every simulation of a design (the temporal ones of base) are
    # left out of its forest
    views <- lapply(views, function(v) v[vapply(v, function(s) any(is.finite(table[[s]])), TRUE)])
    table <- add_derived(table, cs$genetics)
    parameters <- c(names(prior), derived_of(cs$genetics))
    recovery <- lapply(views, function(v) {
      abc_recovery(table, parameters = parameters, statistics = v, n_test = n_test,
                   ntree = ntree, ncores = cores, seed = 7)
    })
    importance <- abc_importance(abc_fit(table, parameters = "log_ne_ratio",
                                         statistics = views$consecutive30, ntree = ntree,
                                         ncores = cores))
    results[[key]] <- list(case = cs, recovery = recovery, importance = importance,
                           minutes = minutes)
    cat("\n==", key, "(", round(minutes, 1), "min of simulation )\n")
    for (v in names(recovery)) {
      r <- recovery[[v]]
      cat("--", v, ":", paste(sprintf("%s %.2f/%.2f", r$parameter, r$r2, r$coverage),
                              collapse = "; "), "\n")
    }
  }
  future::plan(future::sequential)
  summary_table <- do.call(rbind, lapply(results, function(x) {
    do.call(rbind, lapply(names(x$recovery), function(v) {
      data.frame(x$case, view = v,
                 x$recovery[[v]][, c("parameter", "r2", "bias", "coverage", "width")],
                 row.names = NULL)
    }))
  }))
  summary_table$recoverable <- summary_table$r2 >= 0.5 & summary_table$coverage >= 0.85
  rownames(summary_table) <- NULL
  saveRDS(list(summary = summary_table,
               importance = lapply(results, function(x) head(x$importance, 10)),
               minutes = vapply(results, `[[`, 0, "minutes"), cases = series_cases,
               designs = series_designs, n_ref = n_ref, n_test = n_test, date = Sys.Date()),
          output, compress = "xz")
  message("saved ", output)
}
