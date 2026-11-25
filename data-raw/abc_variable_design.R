# Recovery by ABC with random forests when the tree and the landscape vary among simulations,
# with and without imperfect detection. Item (v) of block C of dev/ROADMAP.md. The recovery
# studies (data-raw/abc_recovery.R, abc_detection.R, abc_lambda.R) used one tree of 30 species
# and one grid of 25 sites, drawn once, so their conclusions may hold only for that draw. Here,
# in the design of data-raw/abc_detection.R (same priors and model, series of five surveys, one
# per 50 steps from step 800, three visits with detection log-uniform between 0.02 and 0.5 and
# the largest count kept), four reference tables of N_REF + N_TEST simulations:
#   fixed_perfect     one tree and one landscape, perfect detection;
#   fixed_detection   one tree and one landscape, imperfect detection;
#   varying_perfect   a new tree (treeshapesim::sim_tree_shape(30)) and a new grid
#                     (landscape_grid(side = 5)) in every simulation, perfect detection;
#   varying_detection the same, imperfect detection.
# The two fixed tables are run again because summary_stats() has gained statistics since
# 2026-11-11. Statistics: summary_stats() with the tree, the environment and the coordinates.
# Descriptors of the design, which can be computed for any observed data set: the mean and the
# standard deviation of the phylogenetic distances among all species of the tree (tree_mpd,
# tree_sd_pd) and Moran's I of the environment on the grid with weights 1 / distance
# (env_moran).
# Questions and criterion set before running (R2 >= 0.5 and coverage of the 90% interval of at
# least 85%, as in data-raw/abc_recovery.R after 2026-11-21):
# (i) when the tree and the landscape vary, are the parameters recovered with a fixed design
#     still recovered, and is detection?
# (ii) do the descriptors of the design, added to the statistics, restore what was lost?
# (iii) in the varying tables, does the absolute error of the estimates change with the
#     descriptors (Spearman correlation over the pseudo-observations)? This says whether the
#     results of one tree are typical.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_variable_design.R.
# N_REF, N_TEST, CORES and NTREE change the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
n_ref <- as.integer(Sys.getenv("N_REF", "20000"))
n_test <- as.integer(Sys.getenv("N_TEST", "200"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_variable_design.rds")

draw_design <- function() {
  list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
}
set.seed(20261107)
fixed <- draw_design()
prior <- list(log_dispersal = log10(c(0.001, 0.5)), inter_max = c(0, 1), sigma_mean = c(2, 30),
              lambda = c(0, 1))
prior_detection <- c(prior, list(log_detection = log10(c(0.02, 0.5))))

# one simulation; with `vary`, the tree and the landscape are drawn with the seed of the
# simulation, which abc_table() sets for every row
make_simulate <- function(vary, detection) {
  function(scenario, shared) {
    design <- if (vary) draw_design() else shared
    sim <- simulate_metacommunity(design$landscape,
                                  species_traits(design$tree, scenario$lambda,
                                                 sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                                  competition_matrix(design$tree, "phylo",
                                                     inter_max = scenario$inter_max),
                                  dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                                  record_steps = seq(800, 1000, by = 50))
    if (detection) sim <- observe_detection(sim, 10^scenario$log_detection, visits = 3)
    list(sim = sim, design = design)
  }
}
moran <- function(x, d) {
  w <- 1 / d
  diag(w) <- 0
  z <- x - mean(x)
  length(x) / sum(w) * sum(w * outer(z, z)) / sum(z^2)
}
summarise <- function(out, shared) {
  design <- out$design
  pd <- stats::cophenetic(design$tree)
  pd <- pd[upper.tri(pd)]
  c(summary_stats(out$sim, design$tree, env = design$landscape$env,
                  coords = design$landscape$coords, n_null = 49),
    tree_mpd = mean(pd), tree_sd_pd = stats::sd(pd),
    env_moran = moran(design$landscape$env, design$landscape$dist))
}
cases <- list(
  fixed_perfect = list(prior = prior, vary = FALSE, detection = FALSE),
  fixed_detection = list(prior = prior_detection, vary = FALSE, detection = TRUE),
  varying_perfect = list(prior = prior, vary = TRUE, detection = FALSE),
  varying_detection = list(prior = prior_detection, vary = TRUE, detection = TRUE))
descriptors <- c("tree_mpd", "tree_sd_pd", "env_moran")

if (sys.nframe() == 0) {
  future::plan(future::multicore, workers = cores)
  started <- Sys.time()
  tables <- lapply(seq_along(cases), function(m) {
    abc_table(cases[[m]]$prior, make_simulate(cases[[m]]$vary, cases[[m]]$detection),
              n = n_ref + n_test, setup = function() fixed, summarise = summarise,
              seed = 6000 + m)
  })
  names(tables) <- names(cases)
  future::plan(future::sequential)
  minutes_sim <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  started <- Sys.time()
  views <- function(table) {
    s <- setdiff(attr(table, "statistics"), descriptors)
    # statistics constant in a fixed design carry nothing and are left out
    s <- s[vapply(s, function(v) stats::var(table[[v]], na.rm = TRUE) > 0, logical(1))]
    list(statistics = s, with_design = c(s, descriptors))
  }
  recovery <- list()
  results <- do.call(rbind, lapply(names(tables), function(case) {
    v <- views(tables[[case]])
    if (!cases[[case]]$vary) v <- v["statistics"]
    do.call(rbind, lapply(names(v), function(view) {
      rec <- abc_recovery(tables[[case]], statistics = v[[view]], n_test = n_test,
                          ntree = ntree, ncores = cores, seed = 6100)
      recovery[[paste(case, view)]] <<- rec
      cbind(case = case, view = view,
            rec[, c("parameter", "n", "r2", "bias", "coverage", "width")],
            recoverable = rec$r2 >= 0.5 & rec$coverage >= 0.85,
            conservative = rec$coverage > 0.95)
    }))
  }))
  # (iii) absolute error against the descriptors, in the varying tables, without them in the
  # forest
  errors <- do.call(rbind, lapply(c("varying_perfect", "varying_detection"), function(case) {
    est <- attr(recovery[[paste(case, "statistics")]], "estimates")
    est$error <- abs(est$estimate - est$truth)
    d <- as.data.frame(tables[[case]])
    do.call(rbind, lapply(unique(est$parameter), function(p) {
      e <- est[est$parameter == p, ]
      do.call(rbind, lapply(descriptors, function(x) {
        r <- stats::cor.test(e$error, d[[x]][as.integer(e$name)], method = "spearman",
                             exact = FALSE)
        data.frame(case = case, parameter = p, descriptor = x, rho = unname(r$estimate),
                   p = r$p.value)
      }))
    }))
  }))
  minutes_fit <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  saveRDS(list(results = results, errors = errors,
               design_range = lapply(tables[c("varying_perfect", "varying_detection")],
                                     function(t) summary(as.data.frame(t)[descriptors])),
               minutes = c(simulation = minutes_sim, forests = minutes_fit),
               settings = list(n_ref = n_ref, n_test = n_test, ntree = ntree),
               date = Sys.Date()),
          output, compress = "xz")
  if (nzchar(Sys.getenv("TABLES"))) saveRDS(tables, Sys.getenv("TABLES"))
  print(results, digits = 3, row.names = FALSE)
  print(errors[errors$p < 0.05, ], digits = 3, row.names = FALSE)
  cat("simulation", round(minutes_sim, 1), "min; forests", round(minutes_fit, 1), "min\n")
}
