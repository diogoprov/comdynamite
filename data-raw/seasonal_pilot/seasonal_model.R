# Prototype of a seasonal (monthly-step) metacommunity of pond-breeding anurans, for the
# identifiability pilot of issue #97. Not exported: if the pilot shows that monthly data
# separate the community models, the simulator moves into R/.
#
# One step is one month. Each species has a breeding window (a cyclic function of the month),
# a larval period in months and a niche on one pond axis (e.g. hydroperiod) that acts on
# larval survival. Larvae of all species compete in a pond through a Beverton-Holt term with
# the competition matrix of the community model; metamorphs join the adults of their pond or,
# with probability `disp`, emigrate to other ponds through an exponential kernel. Temporary
# ponds dry for some months each year, which kills their larvae and stops breeding there.
# Field data are monthly counts of tadpoles (dip-netting) and of calling males (only inside
# the breeding window), both with imperfect detection.
#
# The four community models differ only in the competition matrix. SS: interspecific weaker
# than intraspecific, U(0, inter); NT: equal competition and flat niches; PE (strict priority effects): every pair
# U(1, 1 + inter); PM (mixed priority effects): a share p_pri of pairs as in PE, the rest as
# in SS. Coefficients are relative to intraspecific competition (1) and scaled by `a0`.

library(comdynamite)

# competition matrix of a community model, relative to intraspecific competition
seasonal_alpha <- function(k, model, p) {
  a <- matrix(0, k, k)
  up <- upper.tri(a)
  a[up] <- switch(model,
                  SS = stats::runif(sum(up), 0, p$inter),
                  NT = 1,
                  PE = stats::runif(sum(up), 1, 1 + p$inter),
                  PM = ifelse(stats::runif(sum(up)) < p$p_pri,
                              stats::runif(sum(up), 1, 1 + p$inter),
                              stats::runif(sum(up), 0, p$inter_low)))
  a <- a + t(a)
  diag(a) <- 1
  a
}

# prior of a community model; shared parameters have the same prior in every model
seasonal_prior <- function(model) {
  list(model = model,
       disp = 10^stats::runif(1, -3, log10(0.3)),     # share of metamorphs that emigrate
       kern = 10^stats::runif(1, log10(0.1), log10(2)), # kernel scale (km)
       sigma = if (model == "NT") 1e4 else stats::runif(1, 5, 100),
       inter = switch(model, SS = stats::runif(1, 0, 0.8), NT = 1,
                      PE = stats::runif(1, 0.2, 2), PM = stats::runif(1, 0.2, 2)),
       p_pri = if (model == "PM") stats::runif(1, 0.15, 0.5) else NA_real_,
       inter_low = if (model == "PM") stats::runif(1, 0, 0.8) else NA_real_,
       a0 = 10^stats::runif(1, -3, -1.5),             # per-capita larval competition
       fec = 10^stats::runif(1, log10(20), log10(200)), # eggs per adult per breeding month
       s_larva = stats::runif(1, 0.5, 0.9),           # monthly larval survival without crowding
       s_meta = stats::runif(1, 0.05, 0.3),           # survival from metamorphosis to adult
       s_adult = stats::runif(1, 0.85, 0.97),         # monthly adult survival
       kappa = stats::runif(1, 0.5, 4),               # concentration of breeding windows
       q_larva = stats::runif(1, 0.01, 0.1),          # detection per tadpole
       q_adult = stats::runif(1, 0.2, 0.8),           # detection per calling male
       immig = 10^stats::runif(1, -3, -1))            # adults arriving from outside per pond
}

# design: ponds with coordinates (km), an environmental axis and months dry each year
seasonal_design <- function(n_ponds = 13, n_species = 15, temporary = 0.4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  xy <- cbind(stats::runif(n_ponds, 0, 3), stats::runif(n_ponds, 0, 3))
  env <- stats::runif(n_ponds, 0, 100)
  temp <- stats::runif(n_ponds) < temporary
  # temporary ponds dry from the end of the rainy season, for 1 to 4 months
  dry_len <- ifelse(temp, sample(1:4, n_ponds, replace = TRUE), 0)
  dry <- t(vapply(dry_len, function(l) seq_len(12) %in% (((6 + seq_len(l)) - 1) %% 12 + 1),
                  logical(12)))
  list(xy = xy, dist = as.matrix(stats::dist(xy)), env = env, dry = dry,
       n_species = n_species)
}

# species of a run: breeding peak, larval period and niche optimum
seasonal_species <- function(k) {
  data.frame(peak = sample(1:12, k, replace = TRUE),
             larval = sample(2:10, k, replace = TRUE),
             optimum = stats::runif(k, 0, 100))
}

# one run: burn-in years, then monthly sampling of `field_years` years
seasonal_run <- function(design, p, sp, burn_years = 15, field_years = 3) {
  m <- nrow(design$xy)
  k <- nrow(sp)
  amax <- max(sp$larval)
  alpha <- p$a0 * seasonal_alpha(k, p$model, p)
  kernel <- dispersal_kernel(design$dist, p$kern)
  niche <- exp(-outer(design$env, sp$optimum, "-")^2 / (2 * p$sigma^2))   # ponds x species
  window <- vapply(1:12, function(mo) {
    exp(p$kappa * (cos(2 * pi * (mo - sp$peak) / 12) - 1))
  }, numeric(k))                                                           # species x months
  larvae <- array(0L, c(m, k, amax))
  adults <- matrix(stats::rpois(m * k, 2), m, k)
  n_steps <- 12 * (burn_years + field_years)
  first <- 12 * burn_years + 1
  out_l <- out_a <- array(0L, c(m, k, n_steps - first + 1))
  mature <- outer(rep(1, m), sp$larval)                                   # ponds x species
  for (t in seq_len(n_steps)) {
    mo <- (t - 1) %% 12 + 1
    wet <- !design$dry[, mo]
    # breeding into wet ponds
    eggs <- matrix(stats::rpois(m * k, adults * outer(wet, window[, mo]) * p$fec), m, k)
    # crowding by all larvae, then survival with the niche of the pond
    tot <- apply(larvae, c(1, 2), sum) + eggs
    crowd <- 1 + tot %*% t(alpha)
    surv <- pmin(p$s_larva * niche / crowd, 1)
    # ageing: age 1 receives the eggs; survivors move up one month
    for (a in amax:2) {
      larvae[, , a] <- matrix(stats::rbinom(m * k, larvae[, , a - 1], surv), m, k)
    }
    larvae[, , 1] <- matrix(stats::rbinom(m * k, eggs, surv), m, k)
    # metamorphosis at the larval period of each species
    meta <- matrix(0L, m, k)
    for (i in seq_len(k)) {
      meta[, i] <- larvae[, i, sp$larval[i]]
      larvae[, i, sp$larval[i]:amax] <- 0L
    }
    recruits <- matrix(stats::rbinom(m * k, meta, p$s_meta), m, k)
    emig <- matrix(stats::rbinom(m * k, recruits, p$disp), m, k)
    arrive <- crossprod(kernel, emig)
    arrive <- matrix(stats::rpois(m * k, arrive), m, k)
    adults <- matrix(stats::rbinom(m * k, adults, p$s_adult), m, k) + recruits - emig + arrive +
      matrix(stats::rpois(m * k, p$immig / k), m, k)
    # drying kills the larvae of dry ponds
    if (any(!wet)) larvae[!wet, , ] <- 0L
    if (t >= first) {
      j <- t - first + 1
      out_l[, , j] <- matrix(stats::rbinom(m * k, apply(larvae, c(1, 2), sum), p$q_larva), m, k)
      calling <- matrix(stats::rbinom(m * k, adults, pmin(outer(wet, window[, mo]), 1)), m, k)
      out_a[, , j] <- matrix(stats::rbinom(m * k, calling, p$q_adult), m, k)
    }
  }
  list(larvae = out_l, adults = out_a)
}

# summary statistics of monthly counts (ponds x species x months), spatial and temporal
seasonal_stats <- function(x, design, prefix) {
  m <- dim(x)[1]
  n_t <- dim(x)[3]
  pres <- x > 0
  jac <- function(a, b) {
    u <- sum(a | b)
    if (u == 0) NA_real_ else 1 - sum(a & b) / u
  }
  pooled <- apply(pres, c(1, 2), any)                       # ponds x species over all months
  rich <- apply(pres, c(1, 3), sum)                         # ponds x months
  tot <- apply(x, c(1, 3), sum)
  turn <- function(lag) {
    if (n_t <= lag) return(NA_real_)
    mean(unlist(lapply(seq_len(m), function(i) {
      vapply(seq_len(n_t - lag), function(t) jac(pres[i, , t], pres[i, , t + lag]), 0)
    })), na.rm = TRUE)
  }
  before <- pres[, , -n_t, drop = FALSE]
  after <- pres[, , -1, drop = FALSE]
  residents <- array(rich[, -n_t], dim(before)[c(1, 3, 2)])
  residents <- aperm(residents, c(1, 3, 2))
  absent <- !before
  col_ev <- after[absent]
  res_abs <- residents[absent]
  dom <- apply(x, c(1, 3), function(v) if (sum(v) == 0) NA else which.max(v))
  dom_same <- dom[, -1] == dom[, -n_t]
  occ_pond <- apply(pres, c(1, 2), mean)                    # ponds x species
  env_cor <- apply(occ_pond, 2, function(v) {
    if (stats::sd(v) == 0) NA_real_ else abs(stats::cor(v, design$env))
  })
  dj <- as.matrix(vegan::vegdist(pooled * 1, "jaccard", binary = TRUE))
  ok <- rowSums(pooled) > 0
  dd <- if (sum(ok) > 3) suppressWarnings(stats::cor(dj[ok, ok][upper.tri(dj[ok, ok])],
                                                     design$dist[ok, ok][upper.tri(dj[ok, ok])]))
        else NA_real_
  cvs <- apply(tot, 1, function(v) if (mean(v) == 0) NA_real_ else stats::sd(v) / mean(v))
  sync <- suppressWarnings(stats::cor(t(tot)))
  out <- c(alpha = mean(rich), gamma = sum(colSums(pooled) > 0),
           beta = sum(colSums(pooled) > 0) / max(mean(rowSums(pooled)), 1e-9),
           occupancy = mean(pres),
           turn1 = turn(1), turn12 = turn(12),
           col = mean(after[absent]), ext = mean(!after[before]),
           col_res = if (stats::sd(res_abs) > 0 && stats::sd(col_ev) > 0)
             stats::cor(col_ev * 1, res_abs) else NA_real_,
           cv = mean(cvs, na.rm = TRUE), sync = mean(sync[upper.tri(sync)], na.rm = TRUE),
           env = mean(env_cor, na.rm = TRUE), decay = dd,
           dom_persist = mean(dom_same, na.rm = TRUE),
           log_n = log1p(mean(tot)))
  stats::setNames(out, paste0(prefix, "_", names(out)))
}

# statistics of one run for a design with the first `years` field years
seasonal_summarise <- function(run, design, years) {
  keep <- seq_len(12 * years)
  s <- c(seasonal_stats(run$larvae[, , keep, drop = FALSE], design, "tad"),
         seasonal_stats(run$adults[, , keep, drop = FALSE], design, "call"))
  if (years == 1) s <- s[!grepl("turn12", names(s))]
  s
}
