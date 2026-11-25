#' Statistics of trait variation within and among populations
#'
#' Summarises traits measured on individuals of several populations of each species, such as
#' the output of [sample_individuals()] or field data, into statistics that carry information
#' on the genetic variance, on drift and on gene flow, for inference by approximate Bayesian
#' computation ([abc_table()]) next to [summary_stats()]. Each statistic that lacks its inputs
#' is `NA`, so the result always has the same names.
#'
#' A population is a species in a site. The statistics of a species use its populations with at
#' least two individuals measured at the last survey, and need three such populations at sites
#' that differ in environment. `slope` and `r2` are the slope and \eqn{R^2}{R2} of a linear
#' regression of population means on the environment, `sd` the standard deviation of the means
#' and `resid_sd` that of the residuals. `within` is the logarithm of the mean variance within
#' populations (the phenotypic variance, genetic plus environmental) and `within_slope` the slope
#' of the variance within populations on the logarithm of abundance, positive when drift erodes
#' the genetic variance of small populations. `qst` is the residual variance among population
#' means, \eqn{V_B}{V_B}, relative to \eqn{V_B + 2 V_W}{V_B + 2 V_W}, where \eqn{V_W}{V_W} is
#' the mean variance within populations, a divergence among populations beyond what the
#' environment explains, in the form of \eqn{Q_{ST}}{QST} (Spitze 1993) but with phenotypes.
#' `near` and `far` are the semivariances of residuals between populations closer than `near`
#' and farther than `far`, relative to the residual variance, and `dist_cor` the correlation of
#' squared differences of residuals with distance: drift makes neighbours differ, gene flow
#' makes them alike. `size_cor` is the correlation of absolute residuals with the logarithm of
#' abundance, negative when small populations drift more.
#'
#' With two or more surveys (`step`), the first and the last are compared, in the populations
#' with at least two individuals measured in both. `temporal` is the mean squared change of
#' population means minus its sampling variance, relative to the variance within populations,
#' and `temporal_slope` the slope of that corrected change on the inverse of abundance, the
#' basis of the temporal method of estimating the effective size (Waples 1989). Over a long
#' interval, selection and gene flow also move the means, so `temporal_step` uses every pair of
#' consecutive surveys instead: the mean, over these pairs and their populations measured in
#' both, of the squared change minus its sampling variance, divided by the interval between the
#' surveys (in the units of `step`, so per generation when a step is a generation) and by the
#' variance within populations at the last survey. It estimates the variance of drift per
#' generation relative to the phenotypic variance, which is about
#' \eqn{G / (N_e (G + V_E))}{G / (Ne (G + V_E))} (Lande 1976). `temporal_step_slope` is the slope
#' of the corrected change per unit of time on the inverse of the abundance at the later survey.
#' Both need a numeric `step` and three changes in all; with two surveys, `temporal_step` is
#' `temporal` divided by the interval.
#'
#' The statistics of each species are averaged over species, weighted by the summed abundance of
#' their populations at the last survey (or by their number of populations without
#' `abundance`). `slope_min` and `slope_max` are the smallest and largest slope among species.
#' Means and variances are those of the measured individuals, so they include sampling error,
#' as field data would; [subsample_individuals()] changes the effort.
#'
#' @param individuals Data frame with one row per individual.
#' @param trait,site,species,env Names of the columns with the trait (numeric), the site, the
#'   species and the environment of the site (numeric), as in the output of
#'   [sample_individuals()].
#' @param abundance Name of a column with the abundance of the population of each individual,
#'   one value per population, or `NULL` (the default), which leaves `within_slope`, `size_cor`
#'   and `temporal_slope` as `NA`.
#' @param step Name of a column with the survey of each individual, or `NULL` (the default) for
#'   a single survey, which leaves `temporal`, `temporal_slope`, `temporal_step` and
#'   `temporal_step_slope` as `NA`. Without `step`, all
#'   the individuals of a population are pooled.
#' @param coords Matrix or data frame with the coordinates of the sites in two columns and the
#'   names of the sites as row names, or `NULL` (the default), which leaves `near`, `far` and
#'   `dist_cor` as `NA`.
#' @param near,far Distances (non-negative numbers, `near` no larger than `far`) that define
#'   neighbouring and distant pairs of populations.
#' @param by_species If `TRUE`, a data frame with the statistics of each species and its weight
#'   instead of the averages.
#' @return Named numeric vector, or with `by_species = TRUE` a data frame with one row per
#'   species.
#' @references Lande, R. (1976). Natural selection and random genetic drift in phenotypic
#'   evolution. *Evolution*, 30, 314--334. \doi{10.1111/j.1558-5646.1976.tb00911.x}
#'
#'   Spitze, K. (1993). Population structure in *Daphnia obtusa*: quantitative
#'   genetic and allozymic variation. *Genetics*, 135, 367--374.
#'   \doi{10.1093/genetics/135.2.367}
#'
#'   Waples, R. S. (1989). A generalized approach for estimating effective population size from
#'   temporal changes in allele frequency. *Genetics*, 121, 379--391.
#'   \doi{10.1093/genetics/121.2.379}
#' @family trait variation among individuals
#' @examples
#' set.seed(1)
#' land <- landscape_grid(side = 6, env_bandwidth = 2)
#' traits <- data.frame(optimum = c(30, 50, 70), sigma = 15, r_max = 5,
#'                      row.names = c("a", "b", "c"))
#' alpha <- matrix(0.01, 3, 3, dimnames = list(rownames(traits), rownames(traits)))
#' diag(alpha) <- 0.02
#' sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.01, n_steps = 120,
#'                               record_steps = c(100, 120),
#'                               evolution = evolution_control(4, ne_ratio = 0.2, env_var = 1))
#' ind <- sample_individuals(sim, n = 10, source = "evolution", step = c(100, 120))
#' coords <- land$coords
#' rownames(coords) <- rownames(sim$abundance)
#' round(population_stats(ind, abundance = "abundance", step = "step", coords = coords), 3)
#' @export
population_stats <- function(individuals, trait = "trait", site = "site", species = "species",
                             env = "env", abundance = NULL, step = NULL, coords = NULL,
                             near = 1.5, far = 3, by_species = FALSE) {
  for (a in c("trait", "site", "species", "env")) check_string(get(a), a)
  check_string(abundance, "abundance", null_ok = TRUE)
  check_string(step, "step", null_ok = TRUE)
  check_data_columns(individuals, c(trait, site, species, env, abundance, step))
  check_numeric_traits(individuals[[trait]], trait)
  if (!is.numeric(individuals[[env]])) stop("`env` must name a numeric column", call. = FALSE)
  check_number(near, "near", lower = 0)
  check_number(far, "far", lower = near)
  check_flag(by_species, "by_species")
  coords <- population_coords(coords)
  pops <- population_table(individuals, trait, site, species, env, abundance, step)
  series <- if (!is.null(step) && is.numeric(individuals[[step]])) {
    all_surveys <- population_table(individuals, trait, site, species, env, abundance, step,
                                    all = TRUE)
    split(all_surveys, all_surveys$species)
  }
  per <- lapply(split(pops, pops$species), function(p) {
    species_population_stats(p, coords, near, far, series[[p$species[1]]])
  })
  out <- data.frame(species = names(per), weight = vapply(per, attr, 0, "weight"),
                    do.call(rbind, per), row.names = NULL, check.names = FALSE)
  if (by_species) return(out)
  average_population_stats(out)
}

#' names of the statistics of population_stats(), in order
#' @noRd
population_stat_names <- c("slope", "r2", "sd", "resid_sd", "within", "within_slope", "qst",
                           "near", "far", "dist_cor", "size_cor", "temporal", "temporal_slope",
                           "temporal_step", "temporal_step_slope")

#' checks `coords` of population_stats(): NULL, or two numeric columns with site names
#' @noRd
population_coords <- function(coords) {
  if (is.null(coords)) return(NULL)
  coords <- as.matrix(coords)
  bad <- !is.numeric(coords) || ncol(coords) != 2 || is.null(rownames(coords)) ||
    anyNA(coords)
  if (bad) {
    stop("`coords` must have two numeric columns and the names of the sites as row names",
         call. = FALSE)
  }
  coords
}

#' one row per population and survey: mean, variance and number of the measured individuals,
#' abundance and environment; `survey` is 1 for the first survey and 2 for the last (2 also
#' with a single survey), or with `all = TRUE` the rank of every survey, and `time` its step
#' @noRd
population_table <- function(individuals, trait, site, species, env, abundance, step,
                             all = FALSE) {
  x <- individuals[[trait]]
  ok <- is.finite(x)
  survey <- rep(2, length(x))
  if (!is.null(step)) {
    s <- individuals[[step]]
    survey <- if (all) match(s, sort(unique(s))) else
      ifelse(s == max(s, na.rm = TRUE), 2, ifelse(s == min(s, na.rm = TRUE), 1, NA))
  }
  ok <- ok & !is.na(survey)
  key <- paste(individuals[[site]], individuals[[species]], survey, sep = "\r")[ok]
  x <- x[ok]
  first <- which(ok)[!duplicated(key)]
  n <- as.numeric(table(factor(key, unique(key))))
  mean_x <- as.numeric(tapply(x, factor(key, unique(key)), sum)) / n
  sq <- as.numeric(tapply(x^2, factor(key, unique(key)), sum))
  var_x <- ifelse(n >= 2, (sq - n * mean_x^2) / (n - 1), NA_real_)
  ab <- if (is.null(abundance)) rep(NA_real_, length(first)) else
    as.numeric(individuals[[abundance]][first])
  data.frame(site = as.character(individuals[[site]][first]),
             species = as.character(individuals[[species]][first]), survey = survey[first],
             n = n, mean = mean_x, var = var_x, abundance = ab,
             env = individuals[[env]][first],
             time = if (is.null(step)) NA else individuals[[step]][first],
             stringsAsFactors = FALSE)
}

#' statistics of one species from its rows of population_table(), with its weight as an
#' attribute; all NA unless three populations measured twice or more differ in environment
#' @noRd
species_population_stats <- function(p, coords, near, far, series = NULL) {
  last <- p[p$survey == 2 & p$n >= 2, , drop = FALSE]
  weight <- if (all(is.na(p$abundance))) nrow(p[p$survey == 2, ]) else
    sum(p$abundance[p$survey == 2])
  out <- stats::setNames(rep(NA_real_, length(population_stat_names)), population_stat_names)
  if (nrow(last) >= 3 && stats::sd(last$env) > 0) {
    out[] <- c(spatial_population_stats(last, coords, near, far),
               temporal_population_stats(p, last), step_population_stats(series, last))
  } else {
    weight <- 0
  }
  structure(out, weight = weight)
}

#' regression, within, QST-like, distance and size statistics of the last survey
#' @noRd
spatial_population_stats <- function(last, coords, near, far) {
  fit <- stats::lm(last$mean ~ last$env)
  res <- stats::residuals(fit)
  v_res <- stats::var(res)
  w <- mean(last$var)
  log_n <- log(last$abundance)
  c(slope = stats::coef(fit)[[2]], r2 = summary(fit)$r.squared, sd = stats::sd(last$mean),
    resid_sd = sqrt(v_res), within = log(w), within_slope = finite_slope(last$var, log_n),
    qst = v_res / (v_res + 2 * w), distance_population_stats(last$site, res, coords, near, far),
    size_cor = safe_cor(abs(res), log_n))
}

#' semivariances of residuals among near and far pairs, relative to their variance, and the
#' correlation of squared differences with distance
#' @noRd
distance_population_stats <- function(sites, res, coords, near, far) {
  if (is.null(coords) || !all(sites %in% rownames(coords))) {
    return(c(near = NA_real_, far = NA_real_, dist_cor = NA_real_))
  }
  d <- as.matrix(stats::dist(coords[sites, , drop = FALSE]))
  y <- outer(res, res, "-")^2
  low <- lower.tri(d)
  v_res <- stats::var(res)
  semi <- function(sel) if (any(sel)) mean(y[low][sel]) / 2 / v_res else NA_real_
  c(near = semi(d[low] <= near), far = semi(d[low] > far), dist_cor = safe_cor(y[low], d[low]))
}

#' change of population means between the first and the last survey
#' @noRd
temporal_population_stats <- function(p, last) {
  first <- p[p$survey == 1 & p$n >= 2, , drop = FALSE]
  both <- match(paste(last$site), paste(first$site))
  keep <- !is.na(both)
  if (sum(keep) < 3) return(c(temporal = NA_real_, temporal_slope = NA_real_))
  a <- last[keep, ]
  b <- first[both[keep], ]
  delta2 <- (a$mean - b$mean)^2
  noise <- b$var / b$n + a$var / a$n
  c(temporal = (mean(delta2) - mean(noise)) / mean(last$var),
    temporal_slope = finite_slope(delta2 - noise, 1 / a$abundance))
}

#' change of population means between consecutive surveys, per unit of time
#' @noRd
step_population_stats <- function(series, last) {
  none <- c(temporal_step = NA_real_, temporal_step_slope = NA_real_)
  if (is.null(series)) return(none)
  series <- series[series$n >= 2, , drop = FALSE]
  times <- sort(unique(series$time))
  if (length(times) < 2) return(none)
  pairs <- lapply(seq_len(length(times) - 1), function(k) {
    a <- series[series$time == times[k + 1], , drop = FALSE]
    b <- series[series$time == times[k], , drop = FALSE]
    both <- match(a$site, b$site)
    a <- a[!is.na(both), , drop = FALSE]
    b <- b[both[!is.na(both)], , drop = FALSE]
    dt <- times[k + 1] - times[k]
    data.frame(change = ((a$mean - b$mean)^2 - b$var / b$n - a$var / a$n) / dt,
               inverse = 1 / a$abundance)
  })
  pairs <- do.call(rbind, pairs)
  if (nrow(pairs) < 3) return(none)
  c(temporal_step = mean(pairs$change) / mean(last$var),
    temporal_step_slope = finite_slope(pairs$change, pairs$inverse))
}

#' slope of `y` on `x` over the finite pairs; NA with fewer than three or no spread in `x`
#' @noRd
finite_slope <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 3 || stats::sd(x[ok]) == 0) return(NA_real_)
  stats::coef(stats::lm(y[ok] ~ x[ok]))[[2]]
}

#' weighted means over species of the table of population_stats(by_species = TRUE), with the
#' smallest and largest slope
#' @noRd
average_population_stats <- function(out) {
  w <- out$weight
  means <- vapply(population_stat_names, function(s) {
    v <- out[[s]]
    ok <- is.finite(v) & w > 0
    if (!any(ok)) NA_real_ else sum(v[ok] * w[ok]) / sum(w[ok])
  }, 0)
  slopes <- out$slope[is.finite(out$slope)]
  c(means, slope_min = if (length(slopes)) min(slopes) else NA_real_,
    slope_max = if (length(slopes)) max(slopes) else NA_real_)
}
