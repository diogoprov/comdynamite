#' Spatial surrogates of predictors
#'
#' Draws predictors with the same spatial autocorrelation as the observed ones but no relation
#' to anything else, as null models for tests in landscapes with edges, where
#' [torus_permutations()] do not apply or reject too often. They are the `null` options of
#' [pcps_test()], [node_test()] and [phylo_beta_test()], which also take the list returned
#' here (for example drawn with `lonlat` or `strata`), and also serve as null environments for
#' simulations.
#'
#' `"variogram"` follows Viladomat et al. (2014): each predictor is permuted, smoothed with an
#' exponential kernel over a proportion of its nearest sites, with the proportion chosen so that
#' the variogram of the smoothed values best matches the observed one up to the 25% quantile of
#' the distances between sites, and then rescaled, with independent noise for the nugget, so
#' that the scale and intercept of the variogram match too. We use binned variograms rather than
#' their kernel-smoothed one. With several predictors, all of them share the permutation of each
#' draw, but each one gets its own smoothing, so their correlation is kept only approximately.
#'
#' `"parametric"` fits an exponential range by maximum likelihood to the normal scores of the
#' predictor, draws a Gaussian field with that range and gives the observed values to the sites
#' in the order of the field, so each draw is a permutation of sites (in `attr(, "rows")`).
#' With several predictors, the order is that of the first principal component of the scaled
#' predictors, and whole rows move together, which keeps their correlation.
#'
#' In `data-raw/validation_continuous_nulls.R`, for the site test of [pcps_test()] with one
#' predictor and no relation, on a grid of 20 x 20 sites, an L-shaped region of 280 random
#' points and a strip of 50 x 8 sites, with exponential ranges of 2, 5 and 10 sites, the rate of
#' rejection at the 5% level was 0.025-0.09 (mean 0.062) for `"variogram"` and 0.03-0.10 (mean
#' 0.066) for `"parametric"`, against 0.24-0.79 for free permutations and 0.055-0.14 for the
#' torus on the grid. Both became liberal (0.08-0.10) when the range approached the size of the
#' landscape, when few independent patches of environment are left, so with broad fields p
#' values near 0.05 deserve caution. These rates are for one predictor; with several, and in
#' [node_test()] and [phylo_beta_test()], they were not checked.
#'
#' For simulations, `env_surrogates()` gives null environments: a landscape, or an `env_time`
#' matrix, with the same spatial structure but no relation to the environment that shaped the
#' observed communities. To ask whether a simulated pattern depends on where the environment
#' is, rather than only on how much of it is autocorrelated, run the same simulation on the
#' observed landscape and on its surrogates and compare a statistic. With `env_time`,
#' `"parametric"` moves the whole series of each site to another site, so the dynamics of each
#' site are kept and only their place changes; `"variogram"` smooths each step on its own,
#' which does not keep the temporal autocorrelation. Landscapes keep their own distances, such
#' as those of a torus, though the exponential covariance may then not be positive definite.
#'
#' @param env Predictor, or matrix or data frame of predictors, with one row per site, such
#'   as an `env_time` matrix with one column per step. Or an `mc_landscape`, whose
#'   environment (or `env_axes`) is replaced, with the distances of the landscape.
#' @param coords Matrix or data frame of coordinates with one row per site. Not needed for a
#'   landscape.
#' @param n Number of surrogates, one whole number, 1 or more.
#' @param method `"variogram"` or `"parametric"`.
#' @param lonlat Are `coords` longitude and latitude in decimal degrees (`TRUE` or `FALSE`)?
#'   Distances are then great-circle distances in kilometres.
#' @param strata Optional group of each site, such as the city or the basin where it lies.
#'   Surrogates are drawn within each group on its own, so values never move between groups,
#'   which have their own autocorrelation. Each group needs at least 5 sites; a group where
#'   the predictors do not vary keeps them.
#' @return List of `n` matrices with the dimensions and column names of `env`, or of `n`
#'   landscapes when `env` is one. For
#'   `"parametric"`, `attr(, "rows")` is a matrix with one row per draw, giving for each site
#'   the row of `env` it received.
#' @references
#' Viladomat, J., Mazumder, R., McInturff, A., McCauley, D. J., & Hastie, T. (2014). Assessing
#' the significance of global and local correlations under spatial autocorrelation: a
#' nonparametric approach. *Biometrics*, 70, 409--418. \doi{10.1111/biom.12139}
#' @seealso [torus_permutations()], [rigid_permutations()].
#' @family restricted null models
#' @examples
#' set.seed(1)
#' coords <- as.matrix(expand.grid(1:12, 1:12))
#' field <- t(chol(exp(-as.matrix(dist(coords)) / 4))) %*% rnorm(144)
#' draws <- env_surrogates(field, coords, n = 3)
#' length(draws)
#' cor(field, draws[[1]])
#' # null landscapes for the simulator
#' grid <- landscape_grid(10, env_range = c(0, 100))
#' nulls <- env_surrogates(grid, n = 2, method = "parametric")
#' cor(grid$env, nulls[[1]]$env)
#' @export
env_surrogates <- function(env, coords, n = 99, method = c("variogram", "parametric"),
                           lonlat = FALSE, strata = NULL) {
  method <- match.arg(method)
  check_count(n, "n")
  if (inherits(env, "mc_landscape")) return(landscape_surrogates(env, n, method))
  check_flag(lonlat, "lonlat")
  env <- surrogate_env(env)
  coords <- as.matrix(check_coords(coords))
  if (nrow(coords) != nrow(env)) {
    stop("`coords` has ", nrow(coords), " rows and `env` has ", nrow(env), " sites",
         call. = FALSE)
  }
  distances <- function(rows) {
    if (lonlat) great_circle(coords[rows, , drop = FALSE]) else
      as.matrix(stats::dist(coords[rows, , drop = FALSE]))
  }
  if (is.null(strata)) return(surrogates_from(env, distances(seq_len(nrow(env))), n, method))
  if (length(strata) != nrow(env) || anyNA(strata)) {
    stop("`strata` needs one value per site, without NA", call. = FALSE)
  }
  # each stratum on its own, so values never move between strata
  groups <- split(seq_len(nrow(env)), strata, drop = TRUE)
  small <- vapply(groups, length, 1L) < 5
  if (any(small)) {
    stop("strata need at least 5 sites: ", toString(names(groups)[small]),
         call. = FALSE)
  }
  out <- replicate(n, env, simplify = FALSE)
  rows <- matrix(seq_len(nrow(env)), n, nrow(env), byrow = TRUE)
  for (g in groups) {
    local <- env[g, , drop = FALSE]
    if (!all(no_variance(local))) {
      draws <- surrogates_from(local, distances(g), n, method)
      for (r in seq_len(n)) out[[r]][g, ] <- draws[[r]]
      if (method == "parametric") rows[, g] <- matrix(g[attr(draws, "rows")], n)
    }
  }
  if (method == "parametric") attr(out, "rows") <- rows
  out
}

surrogate_env <- function(env) {
  if (is.data.frame(env)) {
    check_numeric_columns(env, "`env`", "Code factors as dummy variables.")
  }
  env <- as.matrix(env)
  if (!is.numeric(env) && !is.logical(env)) stop("`env` is not numeric", call. = FALSE)
  env <- env * 1
  if (anyNA(env)) stop("`env` has missing values", call. = FALSE)
  if (any(is.infinite(env))) stop("`env` has infinite values", call. = FALSE)
  if (all(no_variance(env))) stop("`env` has no variation among sites", call. = FALSE)
  env
}

surrogates_from <- function(env, d, n, method) {
  stopifnot("internal: one row of distances per site of `env`" = nrow(d) == nrow(env))
  if (nrow(env) < 5) stop("surrogates need at least 5 sites in `env`", call. = FALSE)
  if (any(d[upper.tri(d)] == 0)) stop("two sites share the same coordinates", call. = FALSE)
  variable <- !no_variance(env)
  if (method == "parametric") {
    reference <- if (sum(variable) == 1) env[, variable] else {
      stats::prcomp(env[, variable, drop = FALSE], scale. = TRUE)$x[, 1]
    }
    # the fitted range is one whose covariance is positive definite
    root <- t(chol(exp(-d / fit_exponential_range(reference, d)) + diag(1e-8, nrow(d))))
    by_rank <- order(reference)
    rows <- t(replicate(n, by_rank[rank(drop(root %*% stats::rnorm(nrow(d))),
                                        ties.method = "random")]))
    if (n == 1) rows <- matrix(rows, 1)
    draws <- lapply(seq_len(n), function(r) env[rows[r, ], , drop = FALSE])
    attr(draws, "rows") <- rows
    return(draws)
  }
  setup <- variogram_setup(d)
  lapply(seq_len(n), function(r) {
    moved <- sample.int(nrow(env))
    out <- vapply(seq_len(ncol(env)), function(j) {
      if (!variable[j]) return(env[, j])
      variogram_surrogate(env[, j], moved, setup)
    }, numeric(nrow(env)))
    matrix(out, nrow(env), dimnames = dimnames(env))
  })
}

#' landscapes with the environment (and its axes) replaced, with the distances of the landscape,
#' so a torus or a network keeps its own geometry
#' @noRd
landscape_surrogates <- function(landscape, n, method) {
  d <- landscape$dist
  if (is.null(d)) {
    if (is.null(landscape$coords)) {
      stop("the landscape has neither distances nor coordinates", call. = FALSE)
    }
    d <- stats::dist(landscape$coords)
  }
  d <- as.matrix(d)
  axes <- !is.null(landscape$env_axes)
  env <- surrogate_env(if (axes) landscape$env_axes else landscape$env)
  draws <- surrogates_from(env, d, n, method)
  out <- lapply(draws, function(e) {
    moved <- landscape
    if (axes) {
      moved$env_axes <- e
      moved$env <- unname(e[, 1])
    } else {
      moved$env <- unname(e[, 1])
    }
    moved
  })
  attr(out, "rows") <- attr(draws, "rows")
  out
}

#' exponential range by maximum likelihood on normal scores, with the variance profiled out, over
#' a grid of ranges from half the smallest to half the largest distance
#' @noRd
fit_exponential_range <- function(x, d) {
  y <- stats::qnorm(rank(x) / (length(x) + 1))
  y <- y - mean(y)
  positive <- d[upper.tri(d)]
  candidates <- exp(seq(log(min(positive) / 2), log(max(positive) / 2), length.out = 25))
  loglik <- vapply(candidates, function(range) {
    # on a torus or a network, broad ranges may not give a valid covariance
    u <- tryCatch(chol(exp(-d / range) + diag(1e-8, nrow(d))), error = function(e) NULL)
    if (is.null(u)) return(-Inf)
    q <- backsolve(u, y, transpose = TRUE)
    -length(y) / 2 * log(sum(q^2) / length(y)) - sum(log(diag(u)))
  }, numeric(1))
  if (all(loglik == -Inf)) {
    stop("the exponential covariance is not positive definite with these distances", # nocov
         call. = FALSE) # nocov
  }
  candidates[which.max(loglik)]
}

#' what depends only on the coordinates: pairs of sites up to the 25% quantile of the distances,
#' in bins, and the neighbours of each site in order of distance
#' @noRd
variogram_setup <- function(d, proportions = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25,
                                               0.35, 0.5), n_bins = 15) {
  m <- nrow(d)
  pairs <- which(upper.tri(d), arr.ind = TRUE)
  distance <- d[pairs]
  close <- distance <= stats::quantile(distance, 0.25)
  neighbours <- t(apply(d, 1, order))
  sorted <- t(apply(d, 1, sort))
  list(pairs = pairs[close, , drop = FALSE], bin = cut(distance[close], n_bins, labels = FALSE),
       neighbours = neighbours, sorted = sorted,
       k = unique(pmax(2L, as.integer(ceiling(proportions * m)))))
}

binned_variogram <- function(x, setup) {
  g <- 0.5 * (x[setup$pairs[, 1]] - x[setup$pairs[, 2]])^2
  as.numeric(tapply(g, setup$bin, mean))
}

variogram_surrogate <- function(x, moved, setup) {
  target <- binned_variogram(x, setup)
  permuted <- x[moved]
  m <- length(x)
  best <- NULL
  for (k in setup$k) {
    near <- setup$sorted[, seq_len(k), drop = FALSE]
    weights <- exp(-near / pmax(near[, k], 1e-12))
    values <- matrix(permuted[setup$neighbours[, seq_len(k), drop = FALSE]], m)
    smoothed <- rowSums(weights * values) / rowSums(weights)
    g <- binned_variogram(smoothed, setup)
    ok <- is.finite(g) & is.finite(target)
    fit <- stats::lm.fit(cbind(1, g[ok]), target[ok])
    rss <- sum(fit$residuals^2)
    if (is.null(best) || rss < best$rss) {
      best <- list(rss = rss, smoothed = smoothed, coef = fit$coefficients)
    }
  }
  slope <- max(best$coef[2], 0)
  nugget <- max(best$coef[1], 0)
  sqrt(slope) * (best$smoothed - mean(best$smoothed)) + sqrt(nugget) * stats::rnorm(m) + mean(x)
}

#' null draws of the predictors for the occupied sites: permutations of sites (free or given),
#' surrogates drawn here, or surrogates given as a list. `env_at(r)` gives draw r
#' @noRd
site_null <- function(null, draws, landscape_env, env, coords, n_perm) {
  if (is.list(null)) {
    m <- nrow(null[[1]])
    if (m == length(draws$occupied)) {
      null <- lapply(null, function(e) e[draws$occupied, , drop = FALSE])
    } else if (m != sum(draws$occupied)) {
      stop("the surrogates in `null` need one row per site (", length(draws$occupied), ")",
           call. = FALSE)
    }
    if (any(vapply(null, ncol, 1L) != ncol(env))) {
      stop("the surrogates in `null` need one column per predictor (", ncol(env), ")",
           call. = FALSE)
    }
    return(list(n = length(null), env_at = function(r) null[[r]]))
  }
  if (null == "permutation") {
    return(list(n = nrow(draws$rows),
                env_at = function(r) landscape_env[draws$rows[r, ], , drop = FALSE]))
  }
  if (is.null(coords)) stop("`null = \"", null, "\"` needs `coords`", call. = FALSE)
  coords <- as.matrix(check_coords(coords))
  if (nrow(coords) != length(draws$occupied)) {
    stop("`coords` needs one row per site (", length(draws$occupied), ")", call. = FALSE)
  }
  surrogates <- env_surrogates(env, coords[draws$occupied, , drop = FALSE], n_perm, null)
  list(n = n_perm, env_at = function(r) surrogates[[r]])
}

#' a name of a null model, or a list of surrogates from env_surrogates()
#' @noRd
check_null <- function(null, site_permutations) {
  if (is.list(null)) {
    if (length(null) == 0 || !all(vapply(null, function(e) is.matrix(e) && is.numeric(e),
                                         logical(1)))) {
      stop("`null` as a list needs numeric matrices, as from `env_surrogates()`", call. = FALSE)
    }
    if (!is.null(site_permutations)) {
      stop("give either `site_permutations` or surrogates in `null`", call. = FALSE)
    }
    return(null)
  }
  null <- match.arg(null, c("permutation", "variogram", "parametric"))
  if (null != "permutation" && !is.null(site_permutations)) {
    stop("give either `site_permutations` or `null = \"", null, "\"`", call. = FALSE)
  }
  null
}
