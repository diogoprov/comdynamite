# Qualitative check of evolution_control() against the regional neutrality of Leibold et al.
# (2019, PNAS 116: 2612-2617): when populations adapt to their patches faster than gene flow
# swamps them, whichever species arrives first in a patch adapts to it and holds it, so species
# that started with different niches end up with similar regional distributions. Without
# evolution, each species keeps the patches near its own optimum.
# Design (ours, not a replication of their individual-based model): two species with optima 30
# and 70 on an environment U(0, 100), equal competition within and between species, 100 sites
# on a spatially implicit landscape, a grid of dispersal and patch extirpation, with and
# without evolution. Sorting index: difference between the two species in the mean environment
# of the sites they occupy, over the difference of their initial optima (1, sorted by niche; 0,
# the same regional distribution). Run from the root of the package; REPS changes the size.
# Result (10 replicates, 2026-11-12): regional neutrality appears only with patch extirpation,
# low dispersal and fast adaptation. With G = 50, sorting fell from 0.99 without evolution to
# 0.04 (dispersal 0.001, extirpation 0.01) and to 0.22 (dispersal 0.01, extirpation 0.01), while
# the control with the same G and no evolutionary change kept it at 0.90 and 0.73. Without
# extirpation sorting stayed near 1 with any G: the mechanism needs empty patches, where the
# first species to arrive adapts before the other comes. With G = 200 the wider response of
# growth alone lowers sorting (control 0.54-0.56 with extirpation 0.002), so G = 50 is the cleaner
# contrast. Dispersal 0.1 removes sorting with or without evolution (mass effects).
suppressMessages(devtools::load_all(".", quiet = TRUE))
reps <- as.integer(Sys.getenv("REPS", "10"))
traits <- data.frame(optimum = c(30, 70), sigma = 15, r_max = 5, row.names = c("sp1", "sp2"))
alpha <- matrix(0.05, 2, 2, dimnames = list(rownames(traits), rownames(traits)))
# genetic_var 0 is the model without evolution; larger values make adaptation faster
# (a = G / (sigma^2 + G) per generation: 0.04, 0.18 and 0.47 for G = 9, 50 and 200)
grid <- expand.grid(dispersal = c(1e-4, 1e-3, 1e-2, 1e-1), extirpation = c(0, 0.002, 0.01),
                    genetic_var = c(0, 9, 50, 200), rep = seq_len(reps))
grid$response <- TRUE
# control: the same G, and so the same widened response of growth, with no evolutionary change
# (selection and drift divided by 1e9 steps per generation)
control <- expand.grid(dispersal = c(1e-4, 1e-3, 1e-2, 1e-1), extirpation = c(0, 0.002, 0.01),
                       genetic_var = c(50, 200), rep = seq_len(reps), response = FALSE)
grid <- rbind(grid, control)
grid$evolution <- grid$genetic_var > 0
run <- function(i) {
  g <- grid[i, ]
  set.seed(1000 + i)
  ls <- landscape_implicit(100, env_range = c(0, 100))
  ev <- if (g$evolution) evolution_control(genetic_var = g$genetic_var, ne_ratio = 0.5,
                                          steps_per_generation = if (g$response) 1 else 1e9)
  sim <- simulate_metacommunity(ls, traits, alpha, dispersal = g$dispersal,
                                extirpation = g$extirpation, n_steps = 2000, evolution = ev)
  occ <- sim$abundance > 0
  mean_env <- apply(occ, 2, function(x) if (any(x)) mean(ls$env[x]) else NA_real_)
  out <- data.frame(g, sorting = unname(diff(mean_env)) / 40,
                    occupancy_1 = mean(occ[, 1]), occupancy_2 = mean(occ[, 2]))
  if (g$evolution) {
    z <- sim$evolution$genetic
    out$trait_env_cor <- mean(vapply(1:2, function(j) {
      ok <- !is.na(z[, j])
      if (sum(ok) > 3) stats::cor(ls$env[ok], z[ok, j]) else NA_real_
    }, numeric(1)), na.rm = TRUE)
  } else {
    out$trait_env_cor <- NA_real_
  }
  out
}
started <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run,
                                         mc.cores = as.integer(Sys.getenv("CORES", "6"))))
summary_table <- stats::aggregate(cbind(sorting, occupancy_1, occupancy_2, trait_env_cor) ~
                                    dispersal + extirpation + genetic_var + response, data = res,
                                  FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
saveRDS(list(results = res, summary = summary_table, reps = reps,
             minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))),
        "data-raw/evolution_leibold.rds")
print(summary_table, digits = 2)
