# Block A, item 7d: flood-pulse scenarios inspired by the Pantanal (upper Paraguay basin),
# replicated. The set-up (zones, hydrological year, regimes) is described in
# data-raw/river_pantanal_setup.R. Factorial: regime (normal, drought, regulated) x dispersal
# along the channel (symmetric, downstream-biased), N_REP replicates (default 100) on CORES
# cores (default 8). Responses: regional richness every December, local richness and mean
# Jaccard dissimilarity of floodplain and plateau sites at high water (July) and low water
# (November) in years 11-20 (before the drought) and 51-60 (the end), and species persistence
# (persistence(): occupancy and tau) over the last 10 years.
# Run from the package root. Writes data-raw/river_pantanal.rds.
source("data-raw/river_pantanal_setup.R")
n_rep <- as.integer(Sys.getenv("N_REP", "100"))
cores <- as.integer(Sys.getenv("CORES", "8"))
cells <- expand.grid(regime = c("normal", "drought", "regulated"),
                     direction = c("symmetric", "downstream"), stringsAsFactors = FALSE)

jaccard_mean <- function(pa) {
  pa <- pa[rowSums(pa) > 0, , drop = FALSE]
  if (nrow(pa) < 2) return(NA_real_)
  mean(stats::dist(pa, method = "binary"))
}

one_replicate <- function(rep) {
  s <- pantanal_setup(rep)
  rec_year <- rep(seq_len(years), each = period)
  rec_month <- rep(seq_len(period), years)
  rows <- list()
  for (d in c("symmetric", "downstream")) for (g in c("normal", "drought", "regulated")) {
    sim <- pantanal_run(s, d, g, record_steps = burn_years * period + seq_len(years * period))
    h <- sim$history > 0
    regional <- vapply(which(rec_month == 12), function(t) sum(colSums(h[, , t]) > 0),
                       numeric(1))
    pre <- mean(regional[11:20])
    after <- which(regional >= 0.95 * pre & seq_len(years) > max(drought_years))
    summ <- function(sites, m, yrs, f) {
      at <- which(rec_month == m & rec_year %in% yrs)
      mean(vapply(at, function(t) f(h[sites, , t]), numeric(1)), na.rm = TRUE)
    }
    alpha_f <- function(pa) mean(rowSums(pa))
    pers <- persistence(sim, window = 10 * period)
    ok <- is.finite(pers$tau)
    out <- data.frame(rep = rep, regime = g, direction = d,
                      regional_pre = pre, regional_drought_min = min(regional[drought_years]),
                      regional_end = mean(regional[51:60]),
                      recovery_years = if (length(after)) min(after) - max(drought_years)
                                       else NA_real_,
                      n_years_below = sum(regional[33:years] < 0.95 * pre))
    for (w in c("pre", "end")) {
      yrs <- if (w == "pre") 11:20 else 51:60
      for (m in c(high = 7, low = 11)) {
        tag <- paste0(if (m == 7) "high" else "low", "_", w)
        out[[paste0("alpha_plain_", tag)]] <- summ(s$plain, m, yrs, alpha_f)
        out[[paste0("beta_plain_", tag)]] <- summ(s$plain, m, yrs, jaccard_mean)
        out[[paste0("alpha_plateau_", tag)]] <- summ(!s$plain, m, yrs, alpha_f)
        out[[paste0("beta_plateau_", tag)]] <- summ(!s$plain, m, yrs, jaccard_mean)
      }
    }
    out$occupancy <- mean(pers$occupancy[ok])
    out$tau <- mean(pers$tau[ok])
    out$regional_series <- list(regional)
    # mean number of species of each guild per floodplain lake in November (low water), the
    # scale of lagoon surveys
    out$lake_guild_series <- list(vapply(levels(guild), function(k) {
      vapply(which(rec_month == 11), function(t) mean(rowSums(h[s$lakes, guild == k, t])),
             numeric(1))
    }, numeric(years)))
    out$guild_series <- list(vapply(levels(guild), function(k) {
      vapply(which(rec_month == 12), function(t) sum(colSums(h[, guild == k, t]) > 0),
             numeric(1))
    }, numeric(years)))
    rows[[length(rows) + 1]] <- out
  }
  do.call(rbind, rows)
}

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(n_rep), function(r) try(one_replicate(r)), mc.cores = cores)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("replicates failed: ", paste(which(bad), collapse = ", "), "\n", res[bad][[1]])
res <- do.call(rbind, res)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
num <- setdiff(names(res)[vapply(res, is.numeric, logical(1))], "rep")
summary_tab <- stats::aggregate(res[num], res[c("regime", "direction")], mean, na.rm = TRUE)
saveRDS(list(results = res, summary = summary_tab, cells = cells,
             settings = list(n_rep = n_rep, n_sites = n_sites, n_species = n_species,
                             period = period, burn_years = burn_years, years = years,
                             drought_years = drought_years, scale_down = scale_down,
                             overland = overland, emigration = emigration, guild = guild,
                             guild_model = guild_model, sensitivity = sensitivity,
                             migrant_scale = migrant_scale,
                             traits = traits, alpha_inter = alpha[1, 2], high = high,
                             low = low, rain = rain, spate = spate, drying = drying,
                             minutes = elapsed)),
        if (guild_model == "full") "data-raw/river_pantanal.rds" else
          "data-raw/river_pantanal_emigration.rds", version = 2)
print(t(summary_tab), digits = 3)
cat("minutes:", round(elapsed, 1), "\n")
