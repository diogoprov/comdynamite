# Shared set-up of the flood-pulse scenarios inspired by the Pantanal (block A, item 7d), sourced
# by data-raw/river_pantanal.R (the replicated study) and data-raw/animation_river_pantanal.R
# (the animation and the design figure).
# Monthly steps (period 12, step 1 = January). The network is split into a plateau (the half of
# the sites farthest from the outlet, Cerrado tributaries) and a floodplain (the half nearest the
# outlet). Within the floodplain, the sites of largest Shreve magnitude are the main channel; the
# others are floodplain lakes and side channels.
# Hydrological year (Marengo et al. 2021; Harris et al. 2005; Junk et al. 1989; Junk 2005):
# - plateau: rains from December to March bring short, unpredictable spates (polymodal pulse of
#   low-order rivers), each site and month washed out with probability `spate`;
# - floodplain: a monomodal, predictable pulse with high water from April to July, during which
#   floodplain sites are also linked overland (lateral connectivity) and non-migratory species
#   emigrate more (lateral moves into the flooded zone);
# - low water from September to November: floodplain lakes may dry out (lose all populations)
#   with probability `drying`; the main channel is a permanent refuge;
# - outside high water, floodplain lakes are cut off from the river: they neither send nor
#   receive migrants until the next flood (lakes isolated for months or years and connected only
#   at high water, Junk et al. 1989; a dead arm of the Miranda isolated in the dry season until
#   the next flood, Resende 2000). In a drought they stay cut off for 12 years.
# Regimes after 100 years of assembly under the normal pulse:
# - "normal": the pulse every year;
# - "drought": 12 years without the pulse (as in 1962-1973 at Ladario), with harsher drying of
#   floodplain lakes, then the pulse returns;
# - "regulated": the pulse is removed for good (no lateral connectivity, no dispersal pulse, no
#   drying), as downstream of a dam.
# Dispersal along the channel is symmetric or downstream-biased (upstream scale a quarter of the
# downstream one). Spates and drying come from uniforms drawn once per replicate, so that all
# scenarios of a replicate share them (common random numbers).
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_sites <- 50
n_species <- 20
period <- 12
burn_years <- 100
years <- 60
drought_years <- 21:32
n_steps <- (burn_years + years) * period
scale_down <- 2
overland <- 0.5
high <- 4:7
low <- 9:11
rain <- c(12, 1:3)
spate <- 0.05
drying <- c(normal = 0.03, drought = 0.15)

# Biota: a fish-like regional pool of 40 species (see the article for sources). Species optima
# are evenly spaced along the gradient from the floodplain (env about 30) to the plateau (about
# 90), with `skew` > 1 crowding them towards the floodplain. Niche breadth, the skew and
# interspecific competition were calibrated (scratchpad runs of 8 replicates, normal pulse) so
# that, in a December snapshot of year 20, about 17% of the species present occur only on the
# plateau and 17% only in the floodplain, against 17% and 14% for the fish of the upper Paraguay
# basin (Fernando et al. 2025). Species are split into migratory guilds in the shares reported
# for the Miranda River (15% long-distance migrants, 42% intermediate, 43% sedentary; Resende &
# Palmeira 1999, cited by Junk 2005), each spread evenly along the gradient.
n_species <- as.integer(Sys.getenv("N_SPECIES", "40"))
skew <- as.numeric(Sys.getenv("SKEW", "1"))
sigma <- as.numeric(Sys.getenv("SIGMA", "12"))
guild_share <- c(migrant = 0.15, intermediate = 0.42, sedentary = 0.43)
optima <- 20 + 80 * (seq_len(n_species) - 0.5)^skew / n_species^skew
traits <- data.frame(optimum = optima, sigma = sigma, r_max = 3,
                     row.names = sprintf("sp%02d", seq_len(n_species)))
alpha <- matrix(as.numeric(Sys.getenv("ALPHA", "0.005")), n_species, n_species,
                dimnames = list(rownames(traits), rownames(traits)))
diag(alpha) <- 0.1
counts <- round(guild_share * n_species)
counts[["sedentary"]] <- n_species - sum(counts[c("migrant", "intermediate")])
# each guild spread evenly along the gradient
keys <- unlist(lapply(counts, function(k) (seq_len(k) - 0.5) / k))
guild <- factor(rep(names(counts), counts)[order(keys)], names(guild_share))
# emigration per month: base rate, and the rate during the movements of each guild, the
# upstream spawning migration (piracema, October to April) for migrants and the lateral moves
# of high water for the others (Alho 2008); rates are assumptions, only their order and timing
# come from the literature
emigration <- rbind(base = c(migrant = 0.05, intermediate = 0.03, sedentary = 0.01),
                    moving = c(migrant = 0.3, intermediate = 0.1, sedentary = 0.03))
piracema <- c(10:12, 1:4)
# GUILD_MODEL = "full" (the default) also gives the guilds their own routes and tolerances:
# - migrants follow a channel kernel four times longer than the others, the same upstream and
#   downstream whatever `direction`, since the piracema is an active upstream migration;
#   the other guilds follow the kernel of `direction`;
# - the lateral pulse of high water is confined to floodplain sites (`dispersal_site`, emigration
#   three times higher there, for all guilds);
# - a disturbed site loses each population with the sensitivity of its guild: 1 for migrants,
#   0.85 for intermediate and 0.7 for sedentary species, so that, weighted by the guild shares,
#   about 19% of the species survive the drying of a lake, close to the 7 of 34 species (21%)
#   that survived years without flooding in a lagoon of the middle Sao Francisco, where
#   migratory species were the most likely to disappear and most survivors breathe air (Pompeu &
#   Godinho 2006). The same sensitivity applies to the spates of the plateau.
# GUILD_MODEL = "emigration" keeps the earlier model, in which guilds differ only in emigration
# through time (one kernel for all, all populations of a disturbed site lost, the lateral pulse
# in every site).
guild_model <- match.arg(Sys.getenv("GUILD_MODEL", "full"), c("full", "emigration"))
sensitivity <- c(migrant = 1, intermediate = 0.85, sedentary = 0.7)
migrant_scale <- 4 * scale_down

# network, zones, environment and shared uniforms of one replicate
pantanal_setup <- function(rep) {
  set.seed(1000 * rep)
  net <- river_network(n_sites)
  rel <- net$distance / max(net$distance)
  plain <- rel <= stats::median(rel)
  channel <- plain & net$magnitude >= stats::quantile(net$magnitude[plain], 0.75)
  river <- landscape_river(net, scale_down = scale_down)
  base <- env_river(river, n_steps = years * period, tail_up = c(25, 5), nugget = 2,
                    temporal_scale = 60, mean = 30 + 60 * rel)
  amp <- 2 + 6 * (1 - net$upstream_length / max(net$upstream_length))
  yr <- env_seasonal(base, period = period, amplitude = amp, year_sd = 2)
  steps <- seq_len(n_steps)
  list(rep = rep, net = net, plain = plain, channel = channel, lakes = plain & !channel,
       env = cbind(yr[, rep(seq_len(period), burn_years)], yr),
       u_spate = matrix(stats::runif(n_sites * n_steps), n_sites),
       u_dry = matrix(stats::runif(n_sites * n_steps), n_sites),
       month = (steps - 1) %% period + 1, year = (steps - 1) %/% period + 1 - burn_years)
}

# schedule of one regime: kernel index (1 low water, 2 high water), emigration multiplier and
# the sites x steps matrix of disturbances (1 = the site loses all its populations)
pantanal_schedule <- function(s, regime) {
  pulse <- s$year < 1 | regime == "normal" | (regime == "drought" & !s$year %in% drought_years)
  flood_on <- pulse & s$month %in% high
  p_dry <- ifelse(regime == "regulated" & s$year >= 1, 0,
                  ifelse(regime == "drought" & s$year %in% drought_years, drying[["drought"]],
                         drying[["normal"]]))
  p <- matrix(0, n_sites, n_steps)
  p[!s$plain, s$month %in% rain] <- spate
  dry_cols <- which(s$month %in% low)
  p[s$lakes, dry_cols] <- rep(p_dry[dry_cols], each = sum(s$lakes))
  hit <- (!s$plain & s$u_spate < p) | (s$lakes & s$u_dry < p)
  g <- as.character(guild)
  disp <- matrix(emigration["base", g], n_steps, n_species, byrow = TRUE)
  mig <- g == "migrant"
  disp[s$month %in% piracema, mig] <- emigration["moving", "migrant"]
  # lakes cut off outside high water: no emigrants from them (and no immigrants, through the
  # kernels of low water)
  site <- matrix(1, n_sites, n_steps)
  site[s$lakes, !flood_on] <- 0
  if (guild_model == "emigration") {
    disp[flood_on, !mig] <- rep(emigration["moving", g[!mig]], each = sum(flood_on))
    return(list(kernel_step = ifelse(flood_on, 2L, 1L), dispersal_time = disp,
                disturbance = hit * 1, flood_on = flood_on, dispersal_site = site))
  }
  site[s$plain, flood_on] <- 3
  list(kernel_step = cbind(resident = ifelse(flood_on, 2L, 1L),
                           migrant = ifelse(flood_on, 4L, 3L)),
       dispersal_time = disp, disturbance = hit * 1, flood_on = flood_on,
       dispersal_site = site)
}

# dispersal matrices at low water (channel only) and high water (channel plus overland links
# among floodplain sites)
pantanal_kernels <- function(s, direction) {
  land <- landscape_river(s$net, scale_down = scale_down,
                          scale_up = if (direction == "symmetric") scale_down else scale_down / 4)
  mig <- landscape_river(s$net, scale_down = migrant_scale)$dispersal
  straight <- dispersal_kernel(as.matrix(stats::dist(s$net[, c("x", "y")])), 1)
  straight[!s$plain, ] <- 0
  straight[, !s$plain] <- 0
  diag(straight) <- 0
  lat <- straight / pmax(rowSums(straight), 1e-300)
  k_high <- land$dispersal
  k_high[s$plain, ] <- (1 - overland) * land$dispersal[s$plain, ] + overland * lat[s$plain, ]
  # at low water, no site sends migrants to a lake (rows rescaled by the engine)
  cut_off <- function(k) {
    k[, s$lakes] <- 0
    k
  }
  if (guild_model == "emigration") {
    return(list(landscape = land, kernels = list(low = cut_off(land$dispersal), high = k_high)))
  }
  mig_high <- mig
  mig_high[s$plain, ] <- (1 - overland) * mig[s$plain, ] + overland * lat[s$plain, ]
  list(landscape = land, kernels = list(low = cut_off(land$dispersal), high = k_high,
                                        migrant_low = cut_off(mig), migrant_high = mig_high))
}

pantanal_run <- function(s, direction, regime, record_steps) {
  k <- pantanal_kernels(s, direction)
  sched <- pantanal_schedule(s, regime)
  set.seed(1000 * s$rep + 7)
  simulate_metacommunity(k$landscape, traits, alpha,
                         dispersal = unname(emigration["base", as.character(guild)]),
                         n_steps = n_steps,
                         env_time = s$env, disturbance = sched$disturbance,
                         dispersal_time = sched$dispersal_time, kernels = k$kernels,
                         kernel_step = sched$kernel_step,
                         kernel_group = if (guild_model == "full")
                           ifelse(guild == "migrant", "migrant", "resident"),
                         dispersal_site = sched$dispersal_site,
                         disturbance_sensitivity = if (guild_model == "full")
                           unname(sensitivity[as.character(guild)]) else 1,
                         record_steps = record_steps)
}
