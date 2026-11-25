# Which statistics make the genetic variance G and the ratio of effective to census size Ne/N
# identifiable by ABC with random forests when the ancestral optima of the species are unknown,
# and which sampling design does a field study need? Item (viii) of block C of dev/ROADMAP.md,
# from the fourth limitation of the manuscript: in data-raw/abc_genetic_variance.R the ancestral
# optima were the same in every run, so the departure of trait means from them carried
# information that a study of real species would lack, and Ne/N was never recovered.
#
# Design as in data-raw/abc_genetic_variance.R (eight species on a 10 x 10 grid, uniform priors
# on log10 G, log10 dispersal, niche breadth sigma, log10 Ne/N, competition and log10 V_E; static
# environment for 400 steps, or static for 300 and warming by 0.1 per step for 100), crossed
# with:
#   optima    known: ancestral optima evenly spaced from 20 to 80 in every run, as before;
#             unknown: centre (35 to 65) and spread (20 to 70) of the evenly spaced optima drawn
#             from uniform priors in each run, two further parameters;
#   genetics  constant: one G for every population (the engine of the manuscript);
#             mutation: G of each population set by mutation, drift and selection
#             (evolution_control(mutation_var = )), starting at 1, with log10 V_m (-3 to -1)
#             in place of log10 G, so that the variance within populations depends on Ne.
#
# Individuals are measured with error-free phenotypes drawn around the population mean with
# variance G_ij + V_E (G_ij of the population at the last step; the mutation runs do not record
# G at the first survey, so the second survey's G is used for both, a small approximation over
# 50 steps). Statistics of the populations of each species, from population_stats() (with at
# least three sampled populations, then averaged over species weighted by abundance; the run
# saved in data-raw/abc_identifiability.rds computed them with a copy of the same code inside
# this script, checked to give identical values):
#   old   slope and R2 of population means on the environment, SD of the means and of their
#         residuals, minimum and maximum slope (the statistics of the earlier study);
#   new   within    log of the mean within-population variance, and its slope on log abundance;
#         qst       residual variance among population means relative to it plus twice the
#                   within variance, a QST-like divergence beyond the environment;
#         distance  semivariance of residuals between neighbours (distance <= 1.5) and beyond
#                   (> 3), relative to the residual variance, and the correlation of squared
#                   differences with distance (drift raises the first, gene flow lowers it);
#         size      correlation of absolute residuals with log abundance (drift is stronger in
#                   small populations);
#         temporal  (two surveys, steps 350 and 400) change of sample means between surveys,
#                   minus its sampling variance, relative to the within variance, and its slope
#                   on 1 / abundance (the temporal method for Ne).
# Sampling designs computed from the same simulations: 10 individuals per population at every
# occupied site in one survey (base), 5 or 30 individuals, 10 individuals at 20 fixed sites, two
# surveys, and 30 individuals in two surveys at every site; the community survey (summary_stats()
# of abundances, environment and coordinates) is always complete.
#
# Criterion set before running, as in data-raw/abc_recovery.R (revised 2026-11-21): a parameter is
# recoverable with R2 >= 0.5 and coverage of the 90% interval of at least 85%, in 200
# pseudo-observations left out of the table; `conservative` flags coverage above 95%.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_identifiability.R.
# N_REF, N_TEST, CORES and NTREE change the size, PARTS (for example "1,2") the tables run, and
# OUTPUT the file. Each table is also saved in OUTPUT_DIR as it finishes.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_identifiability.rds")
output_dir <- Sys.getenv("OUTPUT_DIR", file.path(tempdir(), "abc_identifiability"))

set.seed(20261123)
side <- 10
land <- landscape_grid(side, env_bandwidth = 3, disp_scale = 1.5)
n_sp <- 8
species <- paste0("sp", seq_len(n_sp))
sites20 <- sort(sample(side^2, 20))
site_names <- paste0("site", seq_len(side^2))
site_coords <- land$coords
rownames(site_coords) <- site_names
surveys <- c(350, 400)
max_ind <- 30
environments <- list(
  static = matrix(land$env, length(land$env), 400),
  warming = sapply(1:400, function(t) land$env + 0.1 * max(0, t - 300))
)
base_prior <- list(log_dispersal = c(-3, -0.5), sigma = c(8, 25), log_ne_ratio = c(-1.3, 0),
                   competition = c(0.2, 1), log_VE = c(-0.5, 1.5))
make_prior <- function(optima, genetics) {
  genetic <- if (genetics == "constant") list(log_G = c(-1, 1.5)) else list(log_Vm = c(-3, -1))
  optimum <- if (optima == "unknown") list(opt_centre = c(35, 65), opt_spread = c(20, 70))
  c(genetic, base_prior, optimum)
}

one_run <- function(scenario, shared) {
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
                         n_steps = 400, env_time = shared$env_time, record_steps = surveys,
                         evolution = evolution)
}

# individuals measured in every population at both surveys: an array of sites x species x
# survey x individual, NA where the population is absent or smaller than the individual number
measure <- function(sim) {
  ev <- sim$evolution
  g <- ev$genetic_variance
  if (is.null(g)) g <- matrix(ev$settings$genetic_var, side^2, n_sp)
  sd_ij <- sqrt(g + ev$settings$env_var)
  sd_ij[!is.finite(sd_ij)] <- 0                 # absent populations, masked below
  abund <- sim$history[, , paste0("step", surveys), drop = FALSE]
  out <- array(NA_real_, c(side^2, n_sp, 2, max_ind))
  for (k in 1:2) {
    means <- ev$history[, , k]
    for (r in seq_len(max_ind)) {
      x <- means + stats::rnorm(length(means), 0, sd_ij)
      x[abund[, , k] < r] <- NA
      out[, , k, r] <- x
    }
  }
  list(ind = out, abund = abund)
}

# individuals of one sampling design as a data frame, as sample_individuals() gives them
individuals_of <- function(m, env, n_ind, sites, two) {
  k <- if (two) 1:2 else 2
  x <- m$ind[sites, , k, seq_len(n_ind), drop = FALSE]
  idx <- arrayInd(seq_along(x), dim(x))
  ok <- is.finite(x)
  data.frame(step = surveys[k][idx[ok, 3]], site = site_names[sites][idx[ok, 1]],
             species = species[idx[ok, 2]],
             abundance = m$abund[sites, , k, drop = FALSE][idx[ok, 1:3, drop = FALSE]],
             env = env[sites][idx[ok, 1]], trait = x[ok])
}

# statistics of one sampling design, from population_stats()
design_stats <- function(m, env, n_ind, sites, two) {
  out <- population_stats(individuals_of(m, env, n_ind, sites, two), abundance = "abundance",
                          step = if (two) "step", coords = site_coords)
  # temporal_step (added in item x) is left out, so that the tables of this study stay the same
  out <- out[setdiff(names(out), c("temporal_step", "temporal_step_slope"))]
  if (!two) out <- out[setdiff(names(out), c("temporal", "temporal_slope"))]
  out[!is.finite(out)] <- NA
  out
}

designs <- list(base = list(n = 10, sites = "all", two = FALSE),
                n5 = list(n = 5, sites = "all", two = FALSE),
                n30 = list(n = 30, sites = "all", two = FALSE),
                sites20 = list(n = 10, sites = "20", two = FALSE),
                surveys2 = list(n = 10, sites = "all", two = TRUE),
                best = list(n = 30, sites = "all", two = TRUE))
old_stats <- c("slope", "r2", "sd", "resid_sd", "slope_min", "slope_max")

summarise_run <- function(sim, shared) {
  env <- sim$evolution$environment
  community <- summary_stats(sim$abundance, env = env, coords = land$coords, n_null = 19)
  names(community) <- paste0("c_", names(community))
  m <- measure(sim)
  pops <- unlist(lapply(names(designs), function(nm) {
    d <- designs[[nm]]
    s <- design_stats(m, env, d$n, if (d$sites == "all") seq_len(side^2) else sites20, d$two)
    stats::setNames(s, paste0(nm, "_", names(s)))
  }))
  c(community, pops)
}

views_of <- function(stats) {
  community <- stats[startsWith(stats, "c_")]
  out <- list(community = community,
              old = c(community, paste0("base_", old_stats)))
  for (nm in names(designs)) out[[nm]] <- c(community, stats[startsWith(stats, paste0(nm, "_"))])
  out
}

cases <- expand.grid(environment = names(environments), genetics = c("constant", "mutation"),
                     optima = c("known", "unknown"), stringsAsFactors = FALSE)
parts <- as.integer(strsplit(Sys.getenv("PARTS", paste(seq_len(nrow(cases)), collapse = ",")),
                             ",")[[1]])

if (sys.nframe() == 0) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  future::plan(future::multicore, workers = cores)
  results <- list()
  for (i in parts) {
    cs <- cases[i, ]
    key <- paste(cs, collapse = "_")
    started <- Sys.time()
    shared <- list(env_time = environments[[cs$environment]])
    table <- abc_table(make_prior(cs$optima, cs$genetics), one_run, n = n_ref + n_test,
                       setup = function() shared, summarise = summarise_run, seed = 3000 + i)
    minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    recovery <- lapply(views_of(attr(table, "statistics")), function(v) {
      abc_recovery(table, statistics = v, n_test = n_test, ntree = ntree, ncores = cores,
                   seed = 7)
    })
    results[[key]] <- list(case = cs, table = table, recovery = recovery, minutes = minutes)
    saveRDS(results[[key]], file.path(output_dir, paste0(key, ".rds")), compress = "xz")
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
                 x$recovery[[v]][, c("parameter", "r2", "bias", "coverage", "width")])
    }))
  }))
  summary_table$recoverable <- summary_table$r2 >= 0.5 & summary_table$coverage >= 0.85
  summary_table$conservative <- summary_table$coverage > 0.95
  rownames(summary_table) <- NULL
  # the reference tables (about 40 MB) stay in the files of OUTPUT_DIR
  saveRDS(list(summary = summary_table, recovery = lapply(results, `[[`, "recovery"),
               minutes = vapply(results, `[[`, 0, "minutes"), cases = cases,
               designs = designs, sites20 = sites20, n_ref = n_ref, n_test = n_test,
               date = Sys.Date()),
          output, compress = "xz")
  message("saved ", output)
}
