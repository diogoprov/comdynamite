#' Maladaptation of populations in a simulation with evolution
#'
#' Measures how far the mean phenotype of each population is from the environment of its site
#' in a simulation run with `evolution` ([evolution_control()]), and the share of growth lost
#' because of that distance (the lag load).
#'
#' The phenotype of a population is its genetic mean plus its plastic response, \eqn{\bar p =
#' \bar z + b (E - E_0)}{pbar = zbar + b (E - E0)}, and its deviation is \eqn{d = E - \bar p}{d
#' = E - pbar}. Under the Gaussian niche of [growth_rates()], with niche breadth
#' \eqn{\sigma}{sigma} and phenotypic variance \eqn{P}{P} around the mean, the growth of the
#' population is that of a population with the same variance centred on the environment times
#' \deqn{1 - L = \exp\left(-\frac{d^2}{2 (\sigma^2 + P)}\right),}{1 - L = exp(-d^2 / (2
#' (sigma^2 + P))),}
#' so the load \eqn{L}{L} goes from 0 (the mean phenotype matches the site) towards 1. It does
#' not depend on `r_max` and it excludes the cost of the variance itself (the standing load).
#' With several evolving traits (one per niche axis), \eqn{d^2 / (\sigma^2 + P)}{d^2 / (sigma^2
#' + P)} becomes the squared Mahalanobis distance \eqn{d^\top (\Omega + P)^{-1} d}{d' (Omega +
#' P)^-1 d} over the axes with a finite niche breadth, which also counts the covariances among
#' traits; axes with `sigma = Inf` add a deviation but no load. For evolving metacommunities
#' ([simulate_evolving_metacommunity()]), the population of a species in a site can hold
#' several lineages: its phenotype is their mean weighted by abundance and \eqn{P}{P} adds the
#' variance among them, so the load is that of a single Gaussian population with that mean and
#' variance, an approximation of the mean load of the lineages.
#'
#' Maladaptation can come from gene flow (immigrants with the means of other sites; see also
#' `migrant_load` in [simulate_metacommunity()]), drift and founder effects in small
#' populations, a lag behind an environment that changes (`env_time`), genetic constraints
#' (little genetic variance in the direction of selection, see [genetic_matrices()]) and
#' plasticity of the wrong sign or size.
#'
#' @param x A simulation from [simulate_metacommunity()] or
#'   [simulate_evolving_metacommunity()] run with `evolution`.
#' @param step `NULL` (the default) for the last step, or one or more steps kept with
#'   `record_steps`.
#' @return A data frame with one row per population present: `step` (with several steps),
#'   `site`, `species`, `abundance`, then `env`, `phenotype` and `deviation` (\eqn{E - \bar
#'   p}{E - pbar}), or with several traits `env_<trait>`, `phenotype_<trait>` and
#'   `deviation_<trait>` for each trait, then `distance` (\eqn{|d| / \sqrt{\sigma^2 + P}}{|d|
#'   / sqrt(sigma^2 + P)}, or the Mahalanobis distance) and `load`.
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' land <- landscape_grid(6, disp_scale = 0.5)
#' traits <- data.frame(optimum = 50, sigma = 15, r_max = 4, row.names = "sp1")
#' alpha <- matrix(0.05, dimnames = list("sp1", "sp1"))
#' sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.1, n_steps = 100,
#'                               evolution = evolution_control(genetic_var = 4))
#' m <- maladaptation(sim)
#' head(m)
#' # sites far from the regional mean keep immigrants from the centre: more maladapted
#' plot(abs(m$env - 50), m$load)
#' @export
maladaptation <- function(x, step = NULL) {
  if (inherits(x, "mc_evolving")) x <- as_metacommunity(x)
  if (!inherits(x, "mc_sim") || is.null(x$evolution$genetic)) {
    stop("`x` must be a simulation run with `evolution`", call. = FALSE)
  }
  if (!is.null(step) && (!is.numeric(step) || !length(step) || anyNA(step))) {
    stop("`step` must be NULL or recorded steps", call. = FALSE)
  }
  if (length(step) > 1) {
    out <- lapply(step, function(s) {
      d <- maladaptation_at(x, s)
      cbind(step = rep(s, nrow(d)), d)
    })
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    return(out)
  }
  maladaptation_at(x, step)
}

maladaptation_at <- function(x, step) {
  recorded <- dimnames(x$evolution$history)
  if (!is.null(step) && step == x$n_steps &&
      !step_label(step) %in% recorded[[length(recorded)]]) {
    step <- NULL
  }
  if (length(dim(x$evolution$genetic)) == 3) return(maladaptation_mv(x, step))
  abundance <- sampled_abundance(x, step)
  sites <- rownames(abundance)
  species <- colnames(abundance)
  niche <- sim_niche(x)
  env <- evolved_environment(x, step) %||% niche$env
  state <- evolved_populations(x, sites, species, step)
  sigma <- niche$traits[species, "sigma"]
  width <- sweep(state$var, 2, sigma^2, "+")
  deviation <- env - state$mean
  distance <- abs(deviation) / sqrt(width)
  distance[, !is.finite(sigma)] <- 0
  present <- abundance > 0 & !is.na(state$mean)
  cell <- which(present)
  i <- (cell - 1) %% length(sites) + 1
  j <- (cell - 1) %/% length(sites) + 1
  out <- data.frame(site = sites[i], species = species[j], abundance = abundance[cell],
                    env = env[i], phenotype = state$mean[cell], deviation = deviation[cell],
                    distance = distance[cell], load = 1 - exp(-distance[cell]^2 / 2),
                    stringsAsFactors = FALSE)
  out <- out[order(match(out$site, sites), match(out$species, species)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

maladaptation_mv <- function(x, step) {
  ev <- x$evolution
  state <- mv_populations(x, step)
  axes <- state$axes
  q <- length(axes)
  abundance <- state$abundance
  sites <- rownames(abundance)
  species <- colnames(abundance)
  sigma <- vapply(x$traits, function(t) as_traits(t)[species, "sigma"], numeric(length(species)))
  sigma <- matrix(sigma, length(species), q)
  present <- abundance > 0 & !is.na(state$mean[, , 1])
  cell <- which(present)
  i <- (cell - 1) %% length(sites) + 1
  j <- (cell - 1) %/% length(sites) + 1
  phen <- matrix(state$mean, ncol = q)[cell, , drop = FALSE]
  dev <- state$env[i, , drop = FALSE] - phen
  distance <- numeric(length(cell))
  for (k in unique(j)) {
    rows <- which(j == k)
    s <- which(is.finite(sigma[k, ]))
    if (!length(s)) next
    g <- match(species[k], names(ev$genetic_var))
    p <- ev$genetic_var[[g]] + ev$env_var[[g]]
    d <- dev[rows, s, drop = FALSE]
    if (is.null(state$within)) {
      A <- solve(diag(sigma[k, s]^2, length(s)) + p[s, s, drop = FALSE])
      distance[rows] <- sqrt(rowSums((d %*% A) * d))
    } else {
      # the variance among lineages in each site adds to the diagonal of P
      w <- matrix(state$within, ncol = q)[cell[rows], , drop = FALSE]
      for (r in seq_along(rows)) {
        A <- solve(diag(sigma[k, s]^2 + w[r, s], length(s)) + p[s, s, drop = FALSE])
        distance[rows[r]] <- sqrt(sum((d[r, ] %*% A) * d[r, ]))
      }
    }
  }
  out <- data.frame(site = sites[i], species = species[j], abundance = abundance[cell],
                    stringsAsFactors = FALSE)
  for (a in seq_len(q)) {
    out[[paste0("env_", axes[a])]] <- state$env[i, a]
    out[[paste0("phenotype_", axes[a])]] <- phen[, a]
    out[[paste0("deviation_", axes[a])]] <- dev[, a]
  }
  out$distance <- distance
  out$load <- 1 - exp(-distance^2 / 2)
  ord <- order(match(out$site, sites), match(out$species, species))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}
