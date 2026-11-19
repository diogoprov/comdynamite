# NAAMP pilot v4: dynamic occupancy with route, stop, year and observer random effects.
# Unit = listening stop (SiteID); season = survey year (unified protocol, 2001-2015);
# replicate = run within year. Stops of routes where the species was ever detected, with
# >= MIN_YEARS surveyed years.
suppressMessages({ library(cmdstanr) })
D <- "data-raw/naamp"
S <- "data-raw/persistence_proxy"
SPECIES <- strsplit(Sys.getenv("SPECIES", "Pseudacris crucifer;Anaxyrus fowleri"), ";")[[1]]
MIN_YEARS <- as.integer(Sys.getenv("MIN_YEARS", "5"))
MAX_ROUTES <- as.integer(Sys.getenv("MAX_ROUTES", "150"))
set.seed(2029)

runs <- read.csv(file.path(D, "Runs.csv")); names(runs) <- trimws(names(runs))
runs <- subset(runs, UnifiedProtocol == 1)
stops <- read.csv(file.path(D, "Stops.csv"))
stops <- merge(stops[stops$SkippedStop == 0, c("RunID", "StopNumber", "SiteID")],
               runs[, c("RunID", "RouteNumber", "SurveyYear", "RunNumber", "ObserverTrackingID", "State")])
counts <- read.csv(file.path(D, "Counts.csv"))
counts <- counts[counts$CallingIndex > 0, c("RunID", "StopNumber", "Species")]
yrs_site <- tapply(stops$SurveyYear, stops$SiteID, function(v) length(unique(v)))
stops <- stops[stops$SiteID %in% names(yrs_site)[yrs_site >= MIN_YEARS], ]
years <- sort(unique(stops$SurveyYear))

mod <- cmdstan_model(file.path(S, "dynocc_re4.stan"))

ppc_stats <- function(nvis, ndet, first) {
  occ <- ndet > 0; obs <- nvis > 0
  naive <- colSums(occ & obs) / pmax(colSums(obs), 1)
  sw <- 0; for (i in seq_len(nrow(occ))) {
    o <- occ[i, obs[i, ]]; if (length(o) > 1) sw <- sw + sum(diff(o) != 0) }
  share <- rowSums(occ & obs) / pmax(rowSums(obs), 1)
  c(mean_naive = mean(naive), switches = sw, ever = mean(rowSums(occ) > 0), sd_stop = sd(share))
}

fit_species <- function(sp) {
  st <- stops
  det <- merge(st, counts[counts$Species == sp, ], by = c("RunID", "StopNumber"))
  routes <- unique(det$RouteNumber)
  if (length(routes) > MAX_ROUTES) routes <- sample(routes, MAX_ROUTES)
  st <- st[st$RouteNumber %in% routes, ]
  st$det <- paste(st$RunID, st$StopNumber) %in% paste(det$RunID, det$StopNumber)
  sites <- sort(unique(st$SiteID))
  st$i <- match(st$SiteID, sites); st$t <- match(st$SurveyYear, years)
  st <- st[order(st$i, st$t, st$RunNumber), ]
  key <- paste(st$i, st$t); ukey <- unique(key)
  vstart <- match(ukey, key); vend <- c(vstart[-1] - 1L, nrow(st))
  ci <- st$i[vstart]; ct <- st$t[vstart]
  cell <- matrix(0L, length(sites), length(years)); cell[cbind(ci, ct)] <- seq_along(ukey)
  route_site <- st$RouteNumber[match(seq_along(sites), st$i)]
  rr <- match(route_site, sort(unique(route_site)))
  oo <- match(st$ObserverTrackingID, sort(unique(st$ObserverTrackingID)))
  gg <- match(st$State, sort(unique(st$State)))
  first <- apply(cell > 0, 1, function(v) which(v)[1])
  dat <- list(N = length(sites), T = length(years), R = max(rr), W = 4L, route = rr, first = first,
              cell = cell, K = length(ukey), vstart = vstart, vend = vend, V = nrow(st),
              y = as.integer(st$det), win = as.integer(st$RunNumber),
              O = max(oo), obs = oo, M = 500L,
              G = max(gg), gst = gg, vyr = st$t)
  nvis <- ndet <- matrix(0L, dat$N, dat$T)
  nvis[cbind(ci, ct)] <- vend - vstart + 1L
  ndet[cbind(ci, ct)] <- tapply(st$det, factor(key, levels = ukey), sum)
  cat(sp, ": stops", dat$N, " routes", dat$R, " stop-years", dat$K, " visits", dat$V, " observers", dat$O, " states", dat$G,
      " naive occupancy", round(mean(ndet[nvis > 0] > 0), 3), "\n")
  t0 <- Sys.time()
  fit <- mod$sample(data = dat, chains = 4, parallel_chains = 4, iter_warmup = as.integer(Sys.getenv("ITER", "1000")),
                    iter_sampling = as.integer(Sys.getenv("ITER", "1000")), refresh = 0, seed = 2029, show_messages = FALSE)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  pars <- c("a_psi", "a_col", "a_ext", "b_p", "s_psi_i", "s_col_r", "s_ext_r", "s_col_i", "s_ext_i",
            "s_col_t", "s_ext_t", "s_p_r", "s_p_i", "s_p_o", "s_sw", "s_yw", "col_m", "ext_m", "eps_m", "tau_m",
            "col_bar", "ext_bar", "eps_bar", "tau_bar")
  summ <- fit$summary(pars)
  diag <- fit$diagnostic_summary(quiet = TRUE)
  dr <- posterior::as_draws_matrix(fit$draws()); dr <- matrix(as.numeric(dr), nrow(dr), dimnames = list(NULL, colnames(dr)))
  obs_stat <- ppc_stats(nvis, ndet, first)
  g <- function(d, nm) dr[d, grep(paste0("^", nm, "\\["), colnames(dr))]
  sims <- t(vapply(sample(nrow(dr), 200), function(d) {
    lc_r <- dr[d, "a_col"] + dr[d, "s_col_r"] * g(d, "z_col_r")[rr] + dr[d, "s_col_i"] * g(d, "z_col_i")
    le_r <- dr[d, "a_ext"] + dr[d, "s_ext_r"] * g(d, "z_ext_r")[rr] + dr[d, "s_ext_i"] * g(d, "z_ext_i")
    lpsi <- dr[d, "a_psi"] + dr[d, "s_psi_i"] * g(d, "z_psi_i")
    lc_t <- dr[d, "s_col_t"] * g(d, "z_col_t"); le_t <- dr[d, "s_ext_t"] * g(d, "z_ext_t")
    offp <- dr[d, "s_p_r"] * g(d, "z_p_r")[rr] + dr[d, "s_p_i"] * g(d, "z_p_i")
    bp <- g(d, "b_p")
    z <- matrix(0L, dat$N, dat$T)
    for (s in seq_len(dat$N)) {
      zz <- rbinom(1, 1, plogis(lpsi[s]))
      for (tt in first[s]:dat$T) {
        if (tt > first[s]) zz <- if (zz == 1) rbinom(1, 1, 1 - plogis(le_r[s] + le_t[tt])) else rbinom(1, 1, plogis(lc_r[s] + lc_t[tt]))
        z[s, tt] <- zz } }
    yv <- rbinom(dat$V, 1, z[cbind(st$i, st$t)] * plogis(bp[dat$win] + offp[st$i] + dr[d, "s_p_o"] * g(d, "z_p_o")[oo] +
      g(d, "s_sw")[dat$win] * matrix(g(d, "z_sw"), dat$G)[cbind(gg, dat$win)] +
      dr[d, "s_yw"] * matrix(g(d, "z_yw"), dat$T)[cbind(st$t, dat$win)]))
    ys <- matrix(0L, dat$N, dat$T); ys[cbind(ci, ct)] <- tapply(yv, factor(key, levels = ukey), sum)
    ppc_stats(nvis, ys, first)
  }, numeric(4)))
  ppc <- rbind(observed = obs_stat, sim_q025 = apply(sims, 2, quantile, 0.025),
               sim_q975 = apply(sims, 2, quantile, 0.975),
               p_value = colMeans(sweep(sims, 2, obs_stat, ">=")))
  gq <- posterior::as_draws_df(fit$draws(c("col_m", "ext_m", "col_bar", "ext_bar", "eps_bar", "tau_bar")))
  list(species = sp, data = dat, summary = summ, gq = gq, iter = Sys.getenv("ITER", "1000"), diagnostics = diag, minutes = mins, ppc = ppc)
}

out <- lapply(SPECIES, fit_species)
for (o in out) {
  cat("\n====", o$species, sprintf("(%.1f min)", o$minutes), "\n")
  print(o$summary[, c("variable", "mean", "q5", "q95", "rhat", "ess_bulk", "ess_tail")], n = 30)
  cat("divergences:", sum(o$diagnostics$num_divergent), " max treedepth hits:",
      sum(o$diagnostics$num_max_treedepth), " E-BFMI:", round(o$diagnostics$ebfmi, 2), "\n")
  cat("PPC:\n"); print(round(o$ppc, 3))
}
saveRDS(lapply(out, function(o) o[names(o) != "data"]), file.path(S, "fits", paste0("pilot4_", gsub(" ", "_", SPECIES[1]), ".rds")))
