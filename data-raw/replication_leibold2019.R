# Replication of the individual-based model of Leibold et al. (2019, PNAS 116: 2612-2617,
# doi:10.1073/pnas.1808615116) with the quantitative-genetic engine of comdynamite.
#
# Their model (SI Appendix, Table S2; bold values): a 10 x 10 grid of patches with two patch
# types (environment 0.4 or 0.6) placed at random, half each; two species, each pre-adapted
# to one type; every patch filled at the start with its pre-adapted species; Gaussian
# stabilising selection of width 0.2; carrying capacity 400 individuals per patch (pooled
# over species, so competition is equal within and between species); exponential dispersal
# kernel with rho = 2 and about Nm = 10 dispersing offspring per patch and generation; patch
# extinction 0.01 per generation; 10,000 generations; 10 replicates. Evolution comes from
# 10 diploid additive loci with mutation 1e-4 and no standing variation at the start.
#
# Translation used here (environment and traits multiplied by 100):
# * patch types 40 and 60, niche width sigma = 20, non-genetic SD 2 (VE = 4);
# * Beverton-Holt growth with r_max = 5 and alpha = 0.01 for all pairs, so K = 400;
# * one step = one generation; dispersal Nm / 400 (0.025 for Nm = 10 emigrants per patch and
#   generation);
#   kernel exp(-distance / 0.5), which is rho = 2 normalised as in their model, no torus;
# * disturbance (whole-patch extinction) 0.01 per step (their Table S2 gives P_ext per
#   generation);
# * evolution of the optimum by Lande's equation with a constant G (no mutation or loss of
#   variance), so the speed of adaptation is set by G. G = 0 is the run without evolution.
# * Control with a single patch type (50) and two identical species (their Fig. S2).
#
# What is not translated: loci, mutation and the build-up of variance from zero; sexual
# reproduction and overlapping generations; microsite variation in the environment.
#
# Response variables, as in their Figs. 3, 5, 6 and S1-S2:
# * proportion of species 1 in the metacommunity and persistence of both species;
# * mean genetic optimum of each species in each patch type (abundance weighted);
# * fitness difference between species within each patch type: mean growth rate of the
#   individuals of species 1 minus that of species 2;
# * variation partitioning (vegan::varpart) of Hellinger-transformed composition, and of
#   the community mean optimum, into patch type and space (PCNM, vegan::pcnm), with the
#   edge patches left out as in their analysis; we keep the adjusted R2 of patch type and of
#   space each alone ([a + b] and [b + c]), since the pure fractions exceed 1 when the
#   adjusted R2 of space is negative;
# * classes of patches at generation 5000: dominated (>= 90% of individuals) by the
#   pre-adapted species, dominated by the other species ("switched"), or shared.

devtools::load_all(quiet = TRUE)

side <- 10
xy <- expand.grid(x = seq_len(side) - 1, y = seq_len(side) - 1)
edge <- xy$x %in% c(0, side - 1) | xy$y %in% c(0, side - 1)
species <- c("sp1", "sp2")
alpha <- matrix(0.01, 2, 2, dimnames = list(species, species))
pcnm_space <- vegan::scores(vegan::pcnm(dist(xy[!edge, ])))

run_leibold <- function(G, landscape_type = c("two", "one"), seed, n_steps = 10000,
                        every = 100, nm = 10, immigrant_selection = FALSE,
                        ne_ratio = 0.5, founder = TRUE, mutation_var = NULL,
                        recruitment = "beverton_holt", mortality = 1,
                        variance_selection = FALSE) {
  landscape_type <- match.arg(landscape_type)
  set.seed(seed)
  env <- if (landscape_type == "two") sample(rep(c(40, 60), side^2 / 2)) else rep(50, side^2)
  optima <- if (landscape_type == "two") c(40, 60) else c(50, 50)
  ls <- landscape_points(xy, env, disp_scale = 0.5)
  traits <- data.frame(optimum = optima, sigma = 20, r_max = 5, row.names = species)
  initial <- if (landscape_type == "two") {
    cbind(sp1 = (env == 40) * 400L, sp2 = (env == 60) * 400L)
  } else {
    first <- sample(rep(c(TRUE, FALSE), side^2 / 2))
    cbind(sp1 = first * 400L, sp2 = (!first) * 400L)
  }
  # with overlapping generations (lottery with mortality < 1), one generation is 1 / mortality
  # steps: n_steps and every stay in generations, and patch extinction is per generation
  per <- 1 / mortality
  steps <- seq(every, n_steps, every)
  sim <- simulate_metacommunity(ls, traits, alpha, dispersal = nm / 400,
                                n_steps = round(n_steps * per),
                                seed_steps = integer(0), initial = initial,
                                record_steps = round(steps * per),
                                disturbance = 1 - (1 - 0.01)^(1 / per),
                                recruitment = recruitment, microsites = 400,
                                mortality = mortality,
                                evolution = evolution_control(
                                  genetic_var = G, env_var = 4, ne_ratio = ne_ratio,
                                  immigrant_selection = immigrant_selection,
                                  founder = founder, mutation_var = mutation_var,
                                  variance_selection = variance_selection))
  list(env = env, preadapted = ifelse(initial[, 1] > 0, 1L, 2L), steps = steps,
       abundance = sim$history, genetic = sim$evolution$history, G = G,
       final = sim$abundance, variance = sim$evolution$genetic_variance,
       landscape_type = landscape_type, seed = seed)
}

# summaries of one run at each recorded step
summarise_leibold <- function(run, var_total = 400 + run$G + 4) {
  types <- sort(unique(run$env))
  rows <- lapply(seq_along(run$steps), function(s) {
    n <- run$abundance[, , s]
    z <- run$genetic[, , s]
    tot <- sum(n)
    out <- data.frame(step = run$steps[s], prop_sp1 = if (tot > 0) sum(n[, 1]) / tot else NA,
                      both = all(colSums(n) > 0), occupied = mean(rowSums(n) > 0))
    for (ty in types) {
      in_ty <- run$env == ty
      for (j in 1:2) {
        w <- n[in_ty, j]
        out[[sprintf("z_sp%d_e%g", j, ty)]] <- if (sum(w) > 0) {
          sum(w * z[in_ty, j], na.rm = TRUE) / sum(w)
        } else NA
      }
      # fitness: growth rate of the mean phenotype in each population, averaged over
      # individuals (the variance term is the same for both species)
      r <- 5 * exp(-(ty - z[in_ty, , drop = FALSE])^2 / (2 * var_total))
      fit <- vapply(1:2, function(j) {
        w <- n[in_ty, j]
        if (sum(w) > 0) sum(w * r[, j], na.rm = TRUE) / sum(w) else NA
      }, numeric(1))
      out[[sprintf("fitdiff_e%g", ty)]] <- fit[1] - fit[2]
    }
    inner <- !edge & rowSums(n) > 0
    out$env_comp <- out$space_comp <- out$env_trait <- out$space_trait <- NA
    if (length(types) == 2 && sum(inner) > 10 && all(colSums(n[inner, ]) > 0)) {
      comp <- vegan::decostand(n[inner, ], "hellinger")
      env_f <- data.frame(type = factor(run$env[inner]))
      sp <- pcnm_space[inner[!edge], , drop = FALSE]
      vp <- vegan::varpart(comp, env_f, sp)$part$fract$Adj.R.squared
      out$env_comp <- vp[1]
      out$space_comp <- vp[2]
      cwm <- rowSums(n[inner, ] * ifelse(is.na(z[inner, ]), 0, z[inner, ])) /
        rowSums(n[inner, ])
      vt <- vegan::varpart(cwm, env_f, sp)$part$fract$Adj.R.squared
      out$env_trait <- vt[1]
      out$space_trait <- vt[2]
    }
    out
  })
  do.call(rbind, rows)
}

classify_leibold <- function(run, step = 5000) {
  s <- match(step, run$steps)
  n <- run$abundance[, , s]
  occ <- rowSums(n) > 0
  share <- n[occ, , drop = FALSE] / rowSums(n[occ, , drop = FALSE])
  dom <- max.col(share, ties.method = "first")
  monop <- apply(share, 1, max) >= 0.9
  cls <- ifelse(!monop, "shared", ifelse(dom == run$preadapted[occ], "preadapted", "switched"))
  table(factor(cls, c("preadapted", "switched", "shared"))) / sum(occ)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_steps <- if (length(args)) as.integer(args[1]) else 10000
  n_rep <- if (length(args) > 1) as.integer(args[2]) else 10
  # G sets the speed of adaptation (G / (sigma^2 + G + VE) of the mismatch per generation:
  # 0, 5%, 11%, 20%, 33% and 50%); Nm as in their Table S2. The same seed (landscape) for a
  # replicate across settings. Control with one patch type and identical species (Fig. S2).
  design <- rbind(expand.grid(G = c(0, 20, 50, 100, 200, 400), nm = c(1, 10, 100),
                              landscape_type = "two", rep = seq_len(n_rep),
                              stringsAsFactors = FALSE),
                  expand.grid(G = c(0, 200), nm = 10, landscape_type = "one",
                              rep = seq_len(n_rep), stringsAsFactors = FALSE))
  t0 <- Sys.time()
  results <- parallel::mclapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    run <- run_leibold(d$G, d$landscape_type, seed = 1000 + d$rep, n_steps = n_steps,
                       nm = d$nm)
    list(G = d$G, nm = d$nm, landscape_type = d$landscape_type, rep = d$rep,
         summary = summarise_leibold(run),
         classes = if (n_steps >= 5000) classify_leibold(run) else NULL)
  }, mc.cores = max(1, parallel::detectCores() - 2), mc.preschedule = FALSE)
  message("time: ", format(Sys.time() - t0))
  out_file <- if (n_steps == 10000 && n_rep == 10) {
    "data-raw/replication_leibold2019.rds"
  } else {
    file.path(tempdir(), "replication_leibold2019_pilot.rds")
  }
  bad <- vapply(results, inherits, logical(1), "try-error")
  if (any(bad)) stop("runs failed: ", results[[which(bad)[1]]])
  saveRDS(results, out_file)
  message("saved ", out_file)
}
