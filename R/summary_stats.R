#' Summary statistics of a metacommunity
#'
#' The same statistics for simulated and observed data, so that they can be compared or used
#' for inference. Each statistic that lacks its inputs is `NA`, so the result always has the
#' same names and results can be bound by rows.
#'
#' `alpha`, `gamma` and `beta` are the mean number of species per site, the number of species
#' in all sites and `gamma / alpha`. `occupancy` and `occupancy_sd` are the mean and the
#' standard deviation of the share of sites each species occupies, over the species present.
#' `abundance` is the mean number of individuals per site, `evenness` the mean Pielou
#' evenness of sites with two or more species, `evenness_regional` the same for all sites
#' together and `sad_sd` the standard deviation of the logarithm of regional abundances, which
#' describes the species abundance distribution. With 0/1 data, these last four describe
#' frequencies of occurrence.
#'
#' `ses_mpd` and `ses_mntd` need `tree`. They are the standardised effect sizes of the mean
#' pairwise distance and of the mean nearest taxon distance, from presences, averaged over sites
#' with two or more species. The null model shuffles species across the tips, as
#' `null.model = "taxa.labels"` of `picante::ses.mpd()`, so the pool is all the species shared
#' by `x` and `tree`, including those without individuals. Compare results only between data
#' with the same pool. This is the only part that uses random numbers.
#'
#' `decay` needs `coords`. It is the correlation between Bray-Curtis dissimilarity and distance
#' among occupied sites. `trait_env` needs `trait` and `env`: the correlation between the mean
#' trait of each site, weighted by abundance, and the environment. `trait_sd` is the standard
#' deviation of the trait among all individuals.
#'
#' The temporal statistics need at least two steps, from `history` or an array. `cv_population`,
#' `cv_site` and `cv_regional` are mean coefficients of variation over time: of each species
#' in each site, of the individuals of each site, and of each species in all sites, leaving out
#' those that are always absent. `turnover` is the mean Bray-Curtis dissimilarity of a site
#' between consecutive steps. Guzman et al. (2022) found that temporal statistics added
#' explanatory power, and that dispersal depended on them.
#'
#' @param x An `mc_sim` object, a matrix with one row per site and one column per species, or
#'   an array of sites x species x steps. In an array, the last step gives the statistics that
#'   are not temporal. From an `mc_sim`, these come from `abundance` and the temporal ones
#'   from `history`, and `env`, `coords` and `trait` default to those of the simulation, the
#'   trait being the environmental optimum.
#' @param tree Optional `phylo` object, matched to species by name.
#' @param trait Optional numeric vector with one value per species, named by species.
#' @param env Optional numeric vector with one value per site.
#' @param coords Optional matrix or data frame with one row per site and the coordinates in
#'   columns.
#' @param n_null Number of shuffles of species for `ses_mpd` and `ses_mntd`.
#' @return Named numeric vector.
#' @references Guzman, L. M., Thompson, P. L., Viana, D. S., Vanschoenwinkel, B., Horváth, Z.,
#'   Ptacnik, R., Jeliazkov, A., Gascón, S., Lemmens, P., Anton-Pardo, M., Langenheder, S., De
#'   Meester, L., & Chase, J. M. (2022). Accounting for temporal change in multiple
#'   biodiversity patterns improves the inference of metacommunity processes. *Ecology*, 103,
#'   e3683. \doi{10.1002/ecy.3683}
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(15)
#' sim <- simulate_metacommunity(landscape_grid(side = 5), species_traits(tree),
#'                               competition_matrix(tree, "phylo"), n_steps = 60,
#'                               seed_steps = 1:5, record_steps = seq(40, 60, by = 5))
#' round(summary_stats(sim, tree, n_null = 19), 3)
#' @export
summary_stats <- function(x, tree = NULL, trait = NULL, env = NULL, coords = NULL, n_null = 99) {
  stopifnot(n_null >= 2)
  series <- NULL
  if (inherits(x, "mc_sim")) {
    comm <- x$abundance
    series <- x$history
    env <- env %||% x$landscape$env %||% x$env
    coords <- coords %||% x$landscape$coords
    if (is.null(trait) && !is.null(x$traits)) {
      trait <- stats::setNames(x$traits$optimum, rownames(x$traits))
    }
  } else if (length(dim(x)) == 3) {
    series <- x
    comm <- x[, , dim(x)[3]]
  } else {
    comm <- x
  }
  comm <- check_comm(comm)
  if (!is.null(env) && (!is.numeric(env) || length(env) != nrow(comm))) {
    stop("`env` needs one number per site", call. = FALSE)
  }
  if (!is.null(coords) && nrow(coords) != nrow(comm)) {
    stop("`coords` needs one row per site", call. = FALSE)
  }
  c(community_stats(comm), phylo_stats(comm, tree, n_null), decay = decay_stat(comm, coords),
    trait_stats(comm, trait, env), temporal_stats(series))
}

community_stats <- function(comm) {
  present <- comm > 0
  richness <- rowSums(present)
  regional <- colSums(comm)
  regional <- regional[regional > 0]
  occupancy <- colMeans(present)[colSums(present) > 0]
  alpha <- mean(richness)
  gamma <- length(regional)
  c(alpha = alpha, gamma = gamma, beta = if (alpha > 0) gamma / alpha else NA_real_,
    occupancy = if (gamma > 0) mean(occupancy) else NA_real_,
    occupancy_sd = if (gamma > 1) stats::sd(occupancy) else NA_real_,
    abundance = mean(rowSums(comm)),
    evenness = mean_or_na(apply(comm, 1, pielou)), evenness_regional = pielou(regional),
    sad_sd = if (gamma > 1) stats::sd(log(regional)) else NA_real_)
}

pielou <- function(abundance) {
  abundance <- abundance[abundance > 0]
  if (length(abundance) < 2) return(NA_real_)
  p <- abundance / sum(abundance)
  -sum(p * log(p)) / log(length(p))
}

mean_or_na <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) mean(values) else NA_real_
}

phylo_stats <- function(comm, tree, n_null) {
  if (is.null(tree)) return(c(ses_mpd = NA_real_, ses_mntd = NA_real_))
  matched <- check_phylo_comm(comm, tree, allow_empty = TRUE)
  presence <- (matched$comm > 0) * 1
  distances <- ape::cophenetic.phylo(matched$tree)[colnames(presence), colnames(presence)]
  null <- lapply(seq_len(n_null), function(r) {
    moved <- sample.int(ncol(presence))
    list(mpd = mean_pairwise(presence, distances[moved, moved]),
         mntd = mean_nearest(presence, distances[moved, moved]))
  })
  effect <- function(observed, statistic) {
    values <- vapply(null, `[[`, numeric(nrow(presence)), statistic)
    values <- matrix(values, nrow(presence))
    mean_or_na((observed - rowMeans(values)) / apply(values, 1, stats::sd))
  }
  c(ses_mpd = effect(mean_pairwise(presence, distances), "mpd"),
    ses_mntd = effect(mean_nearest(presence, distances), "mntd"))
}

# per site, from presences; NA with fewer than two species
mean_pairwise <- function(presence, distances) {
  richness <- rowSums(presence)
  total <- rowSums((presence %*% distances) * presence)
  ifelse(richness > 1, total / (richness * (richness - 1)), NA_real_)
}

mean_nearest <- function(presence, distances) {
  diag(distances) <- Inf
  apply(presence, 1, function(site) {
    members <- which(site > 0)
    if (length(members) < 2) return(NA_real_)
    mean(apply(distances[members, members, drop = FALSE], 1, min))
  })
}

decay_stat <- function(comm, coords) {
  if (is.null(coords)) return(NA_real_)
  occupied <- rowSums(comm) > 0
  if (sum(occupied) < 3) return(NA_real_)
  dissimilarity <- bray_curtis(comm[occupied, , drop = FALSE])
  safe_cor(dissimilarity, as.vector(stats::dist(coords[occupied, , drop = FALSE])))
}

# lower triangle, as in stats::dist()
bray_curtis <- function(comm) {
  sums <- rowSums(comm)
  full <- as.matrix(stats::dist(comm, method = "manhattan")) / outer(sums, sums, "+")
  full[lower.tri(full)]
}

safe_cor <- function(a, b) {
  usable <- is.finite(a) & is.finite(b)
  if (sum(usable) < 3 || stats::sd(a[usable]) == 0 || stats::sd(b[usable]) == 0) return(NA_real_)
  stats::cor(a[usable], b[usable])
}

trait_stats <- function(comm, trait, env) {
  if (is.null(trait)) return(c(trait_env = NA_real_, trait_sd = NA_real_))
  if (!is.numeric(trait) || is.null(names(trait))) {
    stop("`trait` needs one number per species, named by species", call. = FALSE)
  }
  values <- trait[colnames(comm)]
  if (anyNA(values)) stop("`trait` has no value for some species of `x`", call. = FALSE)
  sums <- rowSums(comm)
  regional <- colSums(comm) / sum(comm)
  centre <- sum(regional * values)
  c(trait_env = if (is.null(env)) NA_real_ else {
    safe_cor(drop(comm %*% values)[sums > 0] / sums[sums > 0], env[sums > 0])
  }, trait_sd = if (sum(comm) > 0) sqrt(sum(regional * (values - centre)^2)) else NA_real_)
}

temporal_stats <- function(series) {
  if (is.null(series) || dim(series)[3] < 2) {
    return(c(cv_population = NA_real_, cv_site = NA_real_, cv_regional = NA_real_,
             turnover = NA_real_))
  }
  steps <- dim(series)[3]
  cv <- function(values) {
    means <- rowMeans(values)
    mean_or_na(apply(values, 1, stats::sd)[means > 0] / means[means > 0])
  }
  turnover <- vapply(seq_len(steps - 1), function(t) {
    now <- series[, , t]
    after <- series[, , t + 1]
    rowSums(abs(now - after)) / (rowSums(now) + rowSums(after))
  }, numeric(dim(series)[1]))
  c(cv_population = cv(matrix(series, ncol = steps)), cv_site = cv(apply(series, c(1, 3), sum)),
    cv_regional = cv(apply(series, c(2, 3), sum)), turnover = mean_or_na(turnover))
}
