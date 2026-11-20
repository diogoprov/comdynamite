#' Spatial surrogates of predictors
#'
#' Draws predictors with the same spatial autocorrelation as the observed ones but no relation
#' to anything else, as null models for tests in landscapes with edges, where
#' [torus_permutations()] do not apply or reject too often. They are the `null` options of
#' [pcps_test()], [node_test()] and [phylo_beta_test()], and also serve as null environments for
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
#' @param env Predictor, or matrix or data frame of predictors, with one row per site.
#' @param coords Matrix or data frame of coordinates with one row per site.
#' @param n Number of surrogates.
#' @param method `"variogram"` or `"parametric"`.
#' @return List of `n` matrices with the dimensions and column names of `env`. For
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
#' @export
env_surrogates <- function(env, coords, n = 99, method = c("variogram", "parametric")) {
  method <- match.arg(method)
  env <- as.matrix(check_env(env, matrix(0, NROW(env), 1), collinear = TRUE))
  coords <- as.matrix(check_coords(coords))
  stopifnot(nrow(coords) == nrow(env), n >= 1, n == round(n), nrow(env) >= 5)
  d <- as.matrix(stats::dist(coords))
  if (any(d[upper.tri(d)] == 0)) stop("two sites share the same coordinates", call. = FALSE)
  if (method == "parametric") {
    reference <- if (ncol(env) == 1) env[, 1] else stats::prcomp(env, scale. = TRUE)$x[, 1]
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
    out <- vapply(seq_len(ncol(env)), function(j) variogram_surrogate(env[, j], moved, setup),
                  numeric(nrow(env)))
    out <- matrix(out, nrow(env), dimnames = dimnames(env))
    out
  })
}

# exponential range by maximum likelihood on normal scores, with the variance profiled out, over
# a grid of ranges from half the smallest to half the largest distance
fit_exponential_range <- function(x, d) {
  y <- stats::qnorm(rank(x) / (length(x) + 1))
  y <- y - mean(y)
  positive <- d[upper.tri(d)]
  candidates <- exp(seq(log(min(positive) / 2), log(max(positive) / 2), length.out = 25))
  loglik <- vapply(candidates, function(range) {
    u <- chol(exp(-d / range) + diag(1e-8, nrow(d)))
    q <- backsolve(u, y, transpose = TRUE)
    -length(y) / 2 * log(sum(q^2) / length(y)) - sum(log(diag(u)))
  }, numeric(1))
  candidates[which.max(loglik)]
}

# what depends only on the coordinates: pairs of sites up to the 25% quantile of the distances,
# in bins, and the neighbours of each site in order of distance
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

# null draws of the predictors for the occupied sites: permutations of sites (free or given) or
# surrogates. `env_at(r)` gives draw r
site_null <- function(null, draws, landscape_env, env, coords, n_perm) {
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

check_null <- function(null, site_permutations) {
  null <- match.arg(null, c("permutation", "variogram", "parametric"))
  if (null != "permutation" && !is.null(site_permutations)) {
    stop("give either `site_permutations` or `null = \"", null, "\"`", call. = FALSE)
  }
  null
}
