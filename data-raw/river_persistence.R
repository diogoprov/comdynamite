# Block A, item 7d: long-term persistence in riverine metacommunities.
# Full factorial, N_REP replicates per cell:
# - landscape: dispersal along the channel with equal scales downstream and upstream
#   ("river_symmetric"), with a downstream bias ("river_downstream", upstream scale a quarter of
#   the downstream one), or with the same kernel on straight-line distances among the same sites
#   ("euclidean", the lattice-like control of Carrara et al. 2012);
# - floods: none, washout harsher in large rivers (size_effect = 1) or in headwaters
#   (size_effect = -1), in half of the years on average;
# - network shape: little (branching = 0.2) or much branched (0.6).
# Within a replicate, the network, the environment and the flood years are shared by all cells
# (common random numbers). The environment is a tail-up field along the network (env_river) with
# a gradient from the outlet to the headwaters, slow change among years and a seasonal cycle of
# larger amplitude in headwaters (env_seasonal). Twenty species with optima along the gradient
# compete (Beverton-Holt) for 3000 steps (150 years of 20 steps); the environment repeats its
# first year during the first 1000 steps of assembly.
# Responses: regional richness every year, species persistence over the last 50 years
# (persistence(): occupancy and tau), and local richness and mean Jaccard dissimilarity of
# headwaters and of the other sites at the end.
# Run from the package root with the package installed or loadable; N_REP (default 100) and
# CORES (default 8). Writes data-raw/river_persistence.rds.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_rep <- as.integer(Sys.getenv("N_REP", "100"))
cores <- as.integer(Sys.getenv("CORES", "8"))
n_sites <- 50
n_species <- 20
period <- 20
burn <- 1000
n_steps <- 3000
scale_down <- 2

traits <- data.frame(optimum = seq(25, 95, length.out = n_species), sigma = 15, r_max = 3,
                     row.names = sprintf("sp%02d", seq_len(n_species)))
alpha <- matrix(0.03, n_species, n_species, dimnames = list(rownames(traits), rownames(traits)))
diag(alpha) <- 0.1
cells <- expand.grid(landscape = c("euclidean", "river_symmetric", "river_downstream"),
                     floods = c("none", "mainstem", "headwaters"),
                     branching = c(0.2, 0.6), stringsAsFactors = FALSE)

jaccard_mean <- function(pa) {
  pa <- pa[rowSums(pa) > 0, , drop = FALSE]
  if (nrow(pa) < 2) return(NA_real_)
  mean(stats::dist(pa, method = "binary"))
}

one_replicate <- function(rep) {
  out <- vector("list", 2)
  for (b in seq_along(c(0.2, 0.6))) {
    branching <- c(0.2, 0.6)[b]
    set.seed(1000 * rep + b)
    net <- river_network(n_sites, branching = branching)
    river <- landscape_river(net, scale_down = scale_down)
    rel <- net$distance / max(net$distance)
    base <- env_river(river, n_steps = n_steps - burn, tail_up = c(25, 5), nugget = 2,
                      temporal_scale = 100, mean = 30 + 60 * rel)
    head_amp <- 2 + 6 * (1 - net$upstream_length / max(net$upstream_length))
    year <- env_seasonal(base, period = period, amplitude = head_amp, year_sd = 2)
    first <- year[, seq_len(period)]
    env <- cbind(first[, rep(seq_len(period), burn / period)], year)
    flood_years <- (burn / period) + which(stats::runif((n_steps - burn) / period) < 0.5)
    floods <- function(size_effect) {
      # as in river_floods(): the most exposed site has a washout probability of 0.3
      size <- (net$upstream_length / max(net$upstream_length))^size_effect
      p <- matrix(0, n_sites, n_steps)
      p[, (flood_years - 1) * period + 1] <- 0.3 * size / max(size)
      p
    }
    flood_mats <- list(none = 0, mainstem = floods(1), headwaters = floods(-1))
    lands <- list(
      euclidean = landscape_points(net[, c("x", "y")], env = env[, 1], disp_scale = scale_down),
      river_symmetric = river,
      river_downstream = landscape_river(net, scale_down = scale_down,
                                         scale_up = scale_down / 4))
    rows <- list()
    for (l in names(lands)) for (f in names(flood_mats)) {
      set.seed(1000 * rep + b + 7)
      sim <- simulate_metacommunity(lands[[l]], traits, alpha, dispersal = 0.05,
                                    n_steps = n_steps, env_time = env,
                                    disturbance = flood_mats[[f]],
                                    record_steps = seq(period, n_steps, by = period))
      pa <- sim$abundance > 0
      regional <- apply(sim$history > 0, 3, function(x) sum(colSums(x) > 0))
      pers <- persistence(sim, window = 50)
      ok <- is.finite(pers$tau)
      rows[[length(rows) + 1]] <- data.frame(
        rep = rep, branching = branching, landscape = l, floods = f,
        regional_end = sum(colSums(pa) > 0),
        regional_min = min(regional[-seq_len(burn / period)]),
        alpha_head = mean(rowSums(pa[net$headwater, , drop = FALSE])),
        alpha_other = mean(rowSums(pa[!net$headwater, , drop = FALSE])),
        beta_head = jaccard_mean(pa[net$headwater, , drop = FALSE]),
        beta_other = jaccard_mean(pa[!net$headwater, , drop = FALSE]),
        occupancy = mean(pers$occupancy[ok]), tau = mean(pers$tau[ok]),
        n_headwaters = sum(net$headwater))
      rows[[length(rows)]]$regional_series <- list(regional)
    }
    out[[b]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(n_rep), one_replicate, mc.cores = cores)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("replicates failed: ", paste(which(bad), collapse = ", "))
res <- do.call(rbind, res)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
summary_tab <- stats::aggregate(
  cbind(regional_end, regional_min, alpha_head, alpha_other, beta_head, beta_other, occupancy,
        tau) ~ branching + landscape + floods, res, mean)
saveRDS(list(results = res, summary = summary_tab, cells = cells,
             settings = list(n_rep = n_rep, n_sites = n_sites, n_species = n_species,
                             period = period, burn = burn, n_steps = n_steps,
                             scale_down = scale_down, minutes = elapsed)),
        "data-raw/river_persistence.rds", version = 2)
print(summary_tab, digits = 3)
cat("minutes:", round(elapsed, 1), "\n")
