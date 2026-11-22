#' Spatially implicit landscape
#'
#' Sites without coordinates, with environmental values drawn from a uniform
#' distribution. Dispersal among sites is uniform: emigrants can reach any
#' site, so sites are exchangeable and permuting them is a valid null model.
#'
#' @param n_sites Number of sites: one whole number, 2 or more.
#' @param env_range Range of the uniform distribution of the environment: two finite numbers,
#'   the first smaller.
#' @return An `mc_landscape` object: a list with `env` (one value per site),
#'   and `coords`, `dist` and `dispersal` set to `NULL`.
#' @seealso [landscape_grid()] for a spatially explicit landscape.
#' @family landscapes
#' @examples
#' set.seed(1)
#' landscape_implicit(20)
#' @export
landscape_implicit <- function(n_sites, env_range = c(0, 100)) {
  check_count(n_sites, "n_sites", min = 2)
  check_env_range(env_range)
  new_landscape(env = stats::runif(n_sites, env_range[1], env_range[2]),
                coords = NULL, dist = NULL, dispersal = NULL)
}

#' Spatially explicit landscape on a square grid
#'
#' Sites on a `side` x `side` grid, optionally wrapped as a torus so that no
#' site lies on an edge. The environment is spatially autocorrelated:
#' independent normal deviates are smoothed by a Gaussian kernel with a
#' bandwidth of `env_bandwidth` cells, and then given a uniform marginal
#' distribution on `env_range` by ranks. Emigrants reach the other sites with
#' probabilities from [dispersal_kernel()].
#'
#' @param side Number of cells per side: one whole number, 2 or more.
#' @param env_bandwidth Bandwidth (in cells) of the Gaussian kernel that
#'   smooths the environment: one number greater than 0. Larger values give smoother gradients.
#' @param disp_scale Scale (in cells) of the dispersal kernel: one number greater than 0.
#' @param torus Wrap the grid as a torus? `TRUE` or `FALSE`.
#' @param env_range Range of the uniform marginal distribution of the
#'   environment: two finite numbers, the first smaller.
#' @param disp_shape Shape of the dispersal kernel, `"exponential"` (the default) or
#'   `"gaussian"`; see [dispersal_kernel()].
#' @param n_axes Number of environmental axes: one whole number, 1 or more. With more than one, each
#'   axis is a field
#'   smoothed in the same way, the axes are correlated by `axis_cor` before the conversion to
#'   ranks, and the landscape also holds `env_axes`, a matrix with one column per axis, which
#'   [simulate_metacommunity()] uses with a list of traits when `env` is not given. The first
#'   axis is `env`, the same as with `n_axes = 1` and the same seed.
#' @param axis_cor Correlation between the axes: one number for all pairs, or a correlation
#'   matrix (`n_axes` x `n_axes`). It applies to the smoothed fields, so the correlation of the
#'   uniform axes (by ranks) is slightly lower.
#' @return An `mc_landscape` object: a list with `env`, `coords` (data frame
#'   with columns `x` and `y`), `dist` (distances among sites, toroidal if
#'   `torus = TRUE`), `dispersal` (row-stochastic matrix with a zero
#'   diagonal) and, with `n_axes > 1`, `env_axes`.
#' @seealso [torus_permutations()] for site permutations that keep the spatial
#'   structure of a torus.
#' @family landscapes
#' @examples
#' set.seed(1)
#' ls <- landscape_grid(side = 8)
#' image(matrix(ls$env, 8, 8), asp = 1)
#' @export
landscape_grid <- function(side = 10, env_bandwidth = 2, disp_scale = 2,
                           torus = TRUE, env_range = c(0, 100),
                           disp_shape = c("exponential", "gaussian"), n_axes = 1,
                           axis_cor = 0) {
  disp_shape <- match.arg(disp_shape)
  check_count(side, "side", min = 2)
  check_number(env_bandwidth, "env_bandwidth", lower = 0, open = TRUE)
  check_number(disp_scale, "disp_scale", lower = 0, open = TRUE)
  check_flag(torus, "torus")
  check_env_range(env_range)
  check_count(n_axes, "n_axes")
  n_axes <- as.integer(n_axes)
  if (n_axes > 1) axis_cor <- axis_correlation(axis_cor, n_axes)
  xy <- expand.grid(x = seq_len(side) - 1, y = seq_len(side) - 1)
  dx <- abs(outer(xy$x, xy$x, "-"))
  dy <- abs(outer(xy$y, xy$y, "-"))
  if (torus) {
    dx <- pmin(dx, side - dx)
    dy <- pmin(dy, side - dy)
  }
  d <- sqrt(dx^2 + dy^2)
  smooth <- exp(-d^2 / (2 * env_bandwidth^2))
  z <- drop(smooth %*% stats::rnorm(side^2))
  env <- rank_to_marginal(z, stats::runif(side^2, env_range[1], env_range[2]))
  out <- new_landscape(env = unname(env), coords = xy, dist = d,
                       dispersal = dispersal_kernel(d, disp_scale, disp_shape), torus = torus)
  if (n_axes > 1) {
    extra <- smooth %*% matrix(stats::rnorm(side^2 * (n_axes - 1)), side^2)
    # the first column of the mixture is z itself, since chol() is upper triangular
    mixed <- cbind(z, extra) %*% chol(axis_cor)
    axes <- cbind(env, vapply(seq_len(n_axes)[-1], function(a) {
      rank_to_marginal(mixed[, a], stats::runif(side^2, env_range[1], env_range[2]))
    }, numeric(side^2)))
    dimnames(axes) <- list(NULL, paste0("axis", seq_len(n_axes)))
    out$env_axes <- axes
  }
  out
}

#' stops unless `env_range` is two finite numbers in increasing order
#' @noRd
check_env_range <- function(env_range) {
  if (!is.numeric(env_range) || length(env_range) != 2 || !all(is.finite(env_range)) ||
        env_range[1] >= env_range[2]) {
    stop("`env_range` must be two finite numbers, the first smaller than the second",
         call. = FALSE)
  }
  invisible(env_range)
}

#' stops unless `env` is a numeric vector of finite values, with `n` values when `n` is given
#' @noRd
check_env_vector <- function(env, n = NULL) {
  if (!is.numeric(env) || !all(is.finite(env)) || (!is.null(n) && length(env) != n)) {
    stop("`env` must be a vector of finite numbers",
         if (!is.null(n)) paste0(", one per site (", n, ")"), call. = FALSE)
  }
  invisible(env)
}

#' correlation matrix among n axes from one number or a matrix
#' @noRd
axis_correlation <- function(axis_cor, n) {
  if (length(axis_cor) == 1) {
    if (!is.numeric(axis_cor) || is.na(axis_cor) || axis_cor <= -1 / (n - 1) || axis_cor >= 1) {
      stop("one `axis_cor` must lie between -1 / (n_axes - 1) and 1", call. = FALSE)
    }
    m <- matrix(axis_cor, n, n)
    diag(m) <- 1
    return(m)
  }
  if (!is.matrix(axis_cor) || !identical(dim(axis_cor), c(n, n)) ||
      max(abs(axis_cor - t(axis_cor))) > 1e-10 || max(abs(diag(axis_cor) - 1)) > 1e-10 ||
      min(eigen(axis_cor, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("`axis_cor` must be a positive definite correlation matrix of n_axes x n_axes",
         call. = FALSE)
  }
  unname(axis_cor)
}

#' Landscape from site coordinates
#'
#' Builds a spatially explicit landscape from the coordinates and environment of real sites,
#' such as surveyed lakes, so that simulations share their spatial layout. Emigrants reach the
#' other sites with probabilities from [dispersal_kernel()].
#'
#' @param coords Matrix or data frame with one row per site and two columns of coordinates.
#'   With `lonlat = TRUE`, longitude and latitude in decimal degrees.
#' @param env Numeric vector with the environment of each site.
#' @param disp_scale Scale of the dispersal kernel, in the units of the distances
#'   (kilometres with `lonlat = TRUE`): one number greater than 0.
#' @param disp_shape Shape of the dispersal kernel, `"exponential"` (the default) or
#'   `"gaussian"`; see [dispersal_kernel()].
#' @param lonlat Are the coordinates longitude and latitude? `TRUE` or `FALSE`. Distances are then
#'   great-circle
#'   distances in kilometres, with a spherical Earth of radius 6371 km; otherwise Euclidean.
#' @return An `mc_landscape` object, as from [landscape_grid()], without a torus.
#' @family landscapes
#' @examples
#' xy <- cbind(lon = c(-73, -72.5, -71), lat = c(43, 43.2, 44))
#' ls <- landscape_points(xy, env = c(10, 20, 80), disp_scale = 50, lonlat = TRUE)
#' round(ls$dist)
#' @export
landscape_points <- function(coords, env, disp_scale, lonlat = FALSE,
                             disp_shape = c("exponential", "gaussian")) {
  disp_shape <- match.arg(disp_shape)
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || ncol(coords) != 2 || nrow(coords) < 2 || !all(is.finite(coords))) {
    stop("`coords` must have two columns of finite numbers and at least two rows (sites)",
         call. = FALSE)
  }
  check_env_vector(env, nrow(coords))
  check_number(disp_scale, "disp_scale", lower = 0, open = TRUE)
  check_flag(lonlat, "lonlat")
  d <- if (lonlat) great_circle(coords) else as.matrix(stats::dist(coords))
  dimnames(d) <- NULL
  coords <- as.data.frame(coords)
  rownames(coords) <- NULL
  out <- new_landscape(env = unname(env), coords = coords, dist = d,
                       dispersal = dispersal_kernel(d, disp_scale, disp_shape), torus = FALSE)
  out$points <- TRUE
  out
}

#' Landscape from a dispersal matrix
#'
#' Builds a spatially explicit landscape from a given matrix of movement among sites, for
#' connectivity that distances alone do not describe: a kernel of resistance (least-cost or
#' circuit) distances, estimated movement rates, or a network of ponds linked by corridors.
#' Rows are rescaled to sum to 1 and the diagonal is ignored, so the matrix may hold relative
#' weights.
#'
#' @param dispersal Square, non-negative matrix: row s gives the relative probabilities that an
#'   emigrant from site s reaches each other site. Every site needs at least one destination.
#'   For resistance distances `r`, `dispersal_kernel(r, scale)` gives an exponential kernel.
#' @param env Numeric vector with the environment of each site.
#' @param coords Optional matrix or data frame with two columns of site coordinates, kept for
#'   plotting and for functions that need them.
#' @param dist Optional matrix of distances among sites (for instance, the resistance distances
#'   behind `dispersal`), kept with the landscape.
#' @return An `mc_landscape` object, as from [landscape_points()].
#' @family landscapes
#' @examples
#' # four ponds in a row, where a road cuts the movement between ponds 2 and 3
#' resistance <- as.matrix(dist(1:4))
#' resistance[1:2, 3:4] <- resistance[1:2, 3:4] + 10
#' resistance[3:4, 1:2] <- resistance[3:4, 1:2] + 10
#' ls <- landscape_matrix(dispersal_kernel(resistance, scale = 1), env = c(10, 20, 60, 70),
#'                        dist = resistance)
#' round(ls$dispersal, 3)
#' @export
landscape_matrix <- function(dispersal, env, coords = NULL, dist = NULL) {
  dispersal <- as.matrix(dispersal)
  m <- nrow(dispersal)
  if (!is.numeric(dispersal) || m < 2 || ncol(dispersal) != m || !all(is.finite(dispersal)) ||
        any(dispersal < 0)) {
    stop("`dispersal` must be a square matrix of finite, non-negative numbers with at least ",
         "two sites", call. = FALSE)
  }
  check_env_vector(env, m)
  diag(dispersal) <- 0
  if (any(rowSums(dispersal) == 0)) {
    stop("every site needs at least one destination in `dispersal` (rows ",
         paste(which(rowSums(dispersal) == 0), collapse = ", "), " are empty)", call. = FALSE)
  }
  dispersal <- dispersal / rowSums(dispersal)
  dimnames(dispersal) <- NULL
  if (!is.null(coords)) {
    coords <- as.data.frame(coords)
    if (nrow(coords) != m || ncol(coords) != 2) {
      stop("`coords` must have two columns and one row per site (", m, ")", call. = FALSE)
    }
    rownames(coords) <- NULL
  }
  if (!is.null(dist)) {
    dist <- as.matrix(dist)
    if (!is.numeric(dist) || nrow(dist) != m || ncol(dist) != m) {
      stop("`dist` must be a numeric matrix of ", m, " x ", m, " sites", call. = FALSE)
    }
    dimnames(dist) <- NULL
  }
  out <- new_landscape(env = unname(as.numeric(env)), coords = coords, dist = dist,
                       dispersal = dispersal, torus = FALSE)
  out$matrix <- TRUE
  out
}

#' Dispersal kernel
#'
#' The probability that an emigrant from site s reaches site t is proportional to
#' exp(-d / `scale`) for the exponential kernel or exp(-d^2 / (2 `scale`^2)) for the Gaussian
#' one, where d is the distance between the sites. Emigrants do not stay at their own site:
#' the share that stays is set by the emigration probability (`dispersal` in
#' [simulate_metacommunity()]). For the same `scale`, the Gaussian kernel sends emigrants to
#' near sites more often and to far sites much less often (a thin tail), whereas the
#' exponential kernel keeps more long-distance dispersal. Any other kernel can be passed as a
#' matrix to [landscape_matrix()].
#'
#' @param dist Square matrix of distances among sites.
#' @param scale Scale of the kernel, in the units of `dist`: the mean distance of the
#'   exponential kernel and the standard deviation of the Gaussian one. One number greater
#'   than 0.
#' @param shape `"exponential"` (the default) or `"gaussian"`.
#' @return Row-stochastic matrix with a zero diagonal: row s gives the
#'   probabilities that an emigrant from site s reaches each site.
#' @family landscapes
#' @examples
#' d <- as.matrix(dist(cbind(1:5, 0)))
#' round(dispersal_kernel(d, scale = 1), 2)
#' round(dispersal_kernel(d, scale = 1, shape = "gaussian"), 2)
#' @export
dispersal_kernel <- function(dist, scale, shape = c("exponential", "gaussian")) {
  shape <- match.arg(shape)
  if (!is.matrix(dist) || !is.numeric(dist) || nrow(dist) != ncol(dist) || nrow(dist) < 2 ||
        anyNA(dist)) {
    stop("`dist` must be a square numeric matrix without NA, with at least two sites",
         call. = FALSE)
  }
  check_number(scale, "scale", lower = 0, open = TRUE)
  # log weights, shifted by the row maximum so that far sites do not underflow to zero
  w <- if (shape == "exponential") -dist / scale else -dist^2 / (2 * scale^2)
  diag(w) <- -Inf
  k <- exp(w - apply(w, 1, max))
  k / rowSums(k)
}

#' Site permutations that keep the structure of a torus
#'
#' The translations of a torus combined with the symmetries of the grid keep all distances
#' among sites, and therefore the spatial autocorrelation of the environment. They give a
#' restricted null model for spatially explicit metacommunities (see the vignette), as the
#' torus-translation test of Harms et al. (2001). A square grid has eight symmetries and a
#' rectangular one four. Sites are numbered with the first coordinate running fastest, as in
#' `expand.grid(x, y)` and [landscape_grid()]. A real landscape has edges, and this null
#' treats opposite edges as adjacent: see [rigid_permutations()] and [patch_permutations()].
#'
#' @param side Number of cells per side, as in [landscape_grid()], or two numbers for a
#'   rectangular grid, the first for the coordinate that runs fastest: one or two whole
#'   numbers, 2 or more.
#' @param symmetries Combine the translations with the reflections and rotations of the grid?
#'   `TRUE` or `FALSE`.
#' @param n_perm Number of permutations to keep, drawn without replacement: one whole number,
#'   1 or more. The default, `NULL`, keeps all of them.
#' @return Integer matrix with one permutation of the site indices per row, excluding the
#'   identity (\eqn{8\, \mathrm{side}^2 - 1}{8 * side^2 - 1} rows on a square grid with
#'   `side > 2`).
#' @references Harms, K. E., Condit, R., Hubbell, S. P., & Foster, R. B. (2001). Habitat
#'   associations of trees and shrubs in a 50-ha neotropical forest plot. *Journal of
#'   Ecology*, 89, 947--959.
#' @family restricted null models
#' @examples
#' p <- torus_permutations(4)
#' dim(p)
#' dim(torus_permutations(c(6, 4), symmetries = FALSE))
#' @export
torus_permutations <- function(side, symmetries = TRUE, n_perm = NULL) {
  if (!is.numeric(side) || !length(side) %in% 1:2 || !all(is.finite(side)) || any(side < 2) ||
        any(side != round(side))) {
    stop("`side` must be one or two whole numbers, 2 or more", call. = FALSE)
  }
  check_flag(symmetries, "symmetries")
  check_count(n_perm, "n_perm", null_ok = TRUE)
  sides <- rep_len(side, 2)
  xy <- expand.grid(x = seq_len(sides[1]) - 1, y = seq_len(sides[2]) - 1)
  maps <- list(function(x, y) cbind(x, y))
  if (symmetries) {
    maps <- c(maps, function(x, y) cbind(-x, y), function(x, y) cbind(x, -y),
              function(x, y) cbind(-x, -y))
    if (sides[1] == sides[2]) {
      maps <- c(maps, function(x, y) cbind(y, x), function(x, y) cbind(-y, x),
                function(x, y) cbind(y, -x), function(x, y) cbind(-y, -x))
    }
  }
  perms <- list()
  for (f in maps) for (sx in seq_len(sides[1]) - 1) for (sy in seq_len(sides[2]) - 1) {
    p <- f(xy$x, xy$y)
    p <- cbind((p[, 1] + sx) %% sides[1], (p[, 2] + sy) %% sides[2])
    perms[[length(perms) + 1]] <- as.integer(p[, 1] + sides[1] * p[, 2] + 1)
  }
  perms <- unique(do.call(rbind, perms))
  perms <- perms[!apply(perms, 1, function(p) all(p == seq_len(prod(sides)))), , drop = FALSE]
  if (is.null(n_perm)) return(perms)
  perms[sample.int(nrow(perms), min(n_perm, nrow(perms))), , drop = FALSE]
}

#' Environment that changes over time
#'
#' Adds to the environment of each site a deviation that follows an
#' independent first-order autoregressive process with lag-one
#' autocorrelation `rho` and stationary standard deviation `sd`, starting from
#' its stationary distribution. This is a simple generator of temporal
#' autocorrelation, not the one used by Thompson et al. (2020).
#'
#' @param env Environment of each site, e.g. the `env` of an `mc_landscape`.
#' @param n_steps Number of time steps: one whole number, 1 or more.
#' @param rho Lag-one autocorrelation of the deviations: one number in \[0, 1).
#' @param sd Stationary standard deviation of the deviations: one number, 0 or more.
#' @return Matrix with one row per site and one column per time step, for the
#'   `env_time` argument of [simulate_metacommunity()].
#' @family landscapes
#' @examples
#' set.seed(1)
#' e <- env_over_time(c(20, 50, 80), n_steps = 100, rho = 0.95, sd = 5)
#' matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
#' @export
env_over_time <- function(env, n_steps, rho = 0.9, sd = 5) {
  check_env_vector(env)
  check_count(n_steps, "n_steps")
  check_number(rho, "rho", 0, 1, open = c(FALSE, TRUE))
  check_number(sd, "sd", lower = 0)
  m <- length(env)
  dev <- matrix(0, m, n_steps)
  dev[, 1] <- stats::rnorm(m, 0, sd)
  innov <- sd * sqrt(1 - rho^2)
  for (t in seq_len(n_steps)[-1]) dev[, t] <- rho * dev[, t - 1] + stats::rnorm(m, 0, innov)
  env + dev
}

#' Environment autocorrelated in space and time
#'
#' Simulates the changing environment of Thompson et al. (2020): a Gaussian random field over
#' space and time with exponential covariance, held constant during a burn-in, rescaled to
#' \[0, 1\] and converted to ranks. Unlike [env_over_time()], which gives each site its own
#' AR(1) series, sites close together change together.
#'
#' The covariance between two site-steps is \eqn{\exp(-\sqrt{(d/a_s)^2 + (\Delta t/a_t)^2})},
#' with \eqn{d} the distance between the sites, \eqn{\Delta t} the number of steps between
#' them, and \eqn{a_s} and \eqn{a_t} the spatial and temporal scales. The authors' code
#' (Meta_com_framework v1.0, Zenodo 10.5281/zenodo.3708350) draws this field with
#' `RandomFields::RFsimulate(RMexp(scale = 500))` on coordinates multiplied by 10 and steps
#' in units of one, which gives the defaults \eqn{a_s = 50} and \eqn{a_t = 500} for patches
#' with coordinates between 1 and 100. RandomFields is no longer on CRAN, so the field is
#' drawn here by the spectral method (random Fourier features): the sum of `n_features`
#' cosine waves whose frequencies are drawn from the spectral density of the exponential
#' covariance, a multivariate Student t with one degree of freedom. The covariance is exact
#' in expectation, and the field is close to Gaussian for large `n_features`.
#'
#' As in their code, the field is drawn for the steps after the burn-in, the environment of
#' the first step after it is repeated over the burn-in, the whole field is rescaled to
#' \[0, 1\] and redrawn until the environment of that first step spans more than `min_range`,
#' and the values of all sites and steps, burn-in included, are converted to ranks, so that
#' the environment is uniform on \[0, 1\]. The field is not on a torus, even when the
#' landscape is.
#'
#' @param landscape An `mc_landscape` with coordinates, from [landscape_grid()] or
#'   [landscape_points()], or a two-column matrix or data frame of site coordinates.
#' @param n_steps Number of steps, burn-in included: one whole number, larger than `burn_in`.
#' @param burn_in Number of initial steps with a constant environment: one whole number, 0 or
#'   more.
#' @param spatial_scale,temporal_scale Scales \eqn{a_s} (in units of the coordinates) and
#'   \eqn{a_t} (in steps) of the exponential covariance: each one number greater than 0.
#' @param ranks `TRUE` or `FALSE`. If `TRUE`, convert the environment to ranks over all sites and
#'   steps, as in their code; if `FALSE`, keep the rescaled field.
#' @param min_range Minimum range of the rescaled environment at the first step after the
#'   burn-in; the field is redrawn until it is exceeded: one number in \[0, 1). Use 0 to accept
#'   the first draw.
#' @param n_features Number of cosine waves of the spectral method: one whole number, 1 or
#'   more.
#' @return Matrix with one row per site and one column per step, with values in \[0, 1\], for
#'   the `env_time` argument of [simulate_metacommunity()]. Multiply it to put it on the scale
#'   of the traits.
#' @references Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik, R.,
#'   Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A process-based metacommunity
#'   framework linking local and regional scale community ecology. *Ecology Letters*, 23,
#'   1314--1329. \doi{10.1111/ele.13568}
#' @family landscapes
#' @examples
#' set.seed(1)
#' xy <- cbind(x = runif(30, 1, 100), y = runif(30, 1, 100))
#' e <- env_space_time(xy, n_steps = 600, burn_in = 100)
#' matplot(t(e[1:5, ]), type = "l", lty = 1, xlab = "step", ylab = "environment")
#' @export
env_space_time <- function(landscape, n_steps, burn_in = 200, spatial_scale = 50,
                           temporal_scale = 500, ranks = TRUE, min_range = 0.6,
                           n_features = 2000) {
  xy <- if (inherits(landscape, "mc_landscape")) landscape$coords else landscape
  if (is.null(xy)) stop("`landscape` has no coordinates", call. = FALSE)
  xy <- as.matrix(xy)
  if (ncol(xy) != 2 || !is.numeric(xy) || !all(is.finite(xy))) {
    stop("coordinates in `landscape` must be a numeric matrix with two columns of finite ",
         "values", call. = FALSE)
  }
  check_count(burn_in, "burn_in", min = 0)
  check_count(n_steps, "n_steps", min = burn_in + 1)
  check_number(spatial_scale, "spatial_scale", lower = 0, open = TRUE)
  check_number(temporal_scale, "temporal_scale", lower = 0, open = TRUE)
  check_flag(ranks, "ranks")
  check_number(min_range, "min_range", 0, 1, open = c(FALSE, TRUE))
  check_count(n_features, "n_features")
  m <- nrow(xy)
  steps <- (burn_in + 1):n_steps
  # site-steps on the scale of the covariance, so that it is exp(-|h|)
  p <- cbind(xy[rep(seq_len(m), length(steps)), , drop = FALSE] / spatial_scale,
             rep(steps, each = m) / temporal_scale)
  tries <- 0
  repeat {
    tries <- tries + 1
    if (tries > 100) stop("no field met `min_range` in 100 draws; lower it", call. = FALSE)
    z <- exp_field(p, n_features)
    z <- (z - min(z)) / (max(z) - min(z))
    first <- z[seq_len(m)]
    if (max(first) - min(first) > min_range) break
  }
  env <- cbind(matrix(first, m, burn_in), matrix(z, m))
  if (ranks) env[] <- stats::ecdf(env)(env)
  unname(env)
}

#' the environment of env_seasonal() and env_pulse() as a matrix of sites x steps, from a
#' vector repeated over `n_steps` steps or from a matrix
#' @noRd
env_steps <- function(env, n_steps) {
  check_count(n_steps, "n_steps", null_ok = TRUE)
  if (is.matrix(env)) {
    if (!is.null(n_steps) && n_steps != ncol(env)) {
      stop("`n_steps` differs from the number of columns of `env`", call. = FALSE)
    }
    return(env)
  }
  if (!is.numeric(env) || anyNA(env)) {
    stop("`env` must be a numeric vector without NA, or a matrix (sites x steps)",
         call. = FALSE)
  }
  if (is.null(n_steps)) stop("`n_steps` is needed when `env` is a vector", call. = FALSE)
  matrix(env, length(env), n_steps)
}

#' Gaussian field with covariance exp(-|h|) at the rows of `p` (points in any dimension), by
#' random Fourier features: frequencies from the spectral density of the exponential
#' covariance, a multivariate Student t with one degree of freedom. Unit variance.
#' @noRd
exp_field <- function(p, n_features) {
  d <- ncol(p)
  w <- matrix(stats::rnorm(d * n_features), n_features) / abs(stats::rnorm(n_features))
  phase <- stats::runif(n_features, 0, 2 * pi)
  z <- numeric(nrow(p))
  for (s in seq(1, nrow(p), by = 20000)) {
    i <- s:min(nrow(p), s + 19999)
    arg <- p[i, , drop = FALSE] %*% t(w) + rep(phase, each = length(i))
    z[i] <- rowSums(cos(arg))
  }
  z * sqrt(2 / n_features)
}

#' Seasonal cycle of the environment
#'
#' Adds a seasonal cycle to the environment of each site,
#' \eqn{a_s \sin(2 \pi (t - \phi) / P)}{a_s sin(2 pi (t - phi) / P)}, with an amplitude per site
#' (for instance larger in shallow ponds or in headwaters) and, optionally, an anomaly per year
#' shared by all sites (wet and dry years). The environment can already change in time (an
#' `env_time` matrix, e.g. from [env_space_time()] with `ranks = FALSE` or [env_river()] with
#' `n_steps`), which then gives variation among years that differs among sites; floods can be
#' added with [env_pulse()] or, as losses of populations, with [river_floods()].
#'
#' @param env Environment of each site (a vector), or a matrix (sites x steps).
#' @param n_steps Number of steps: one whole number, 1 or more. Taken from `env` when it is a
#'   matrix.
#' @param period Steps per cycle (a year): one number greater than 0.
#' @param amplitude Amplitude of the cycle: one finite value, or one per site.
#' @param phase Step at which the cycle crosses its mean on the way up: one finite number.
#' @param year_sd Standard deviation of an anomaly drawn for each year (each block of `period`
#'   steps) and added to all sites: one number, 0 or more; 0 (the default) for none.
#' @return Matrix with one row per site and one column per step, for `env_time` in
#'   [simulate_metacommunity()].
#' @family landscapes
#' @examples
#' set.seed(1)
#' e <- env_seasonal(c(10, 20, 30), n_steps = 100, period = 25, amplitude = c(1, 3, 6),
#'                   year_sd = 1)
#' matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
#' @export
env_seasonal <- function(env, n_steps = NULL, period, amplitude, phase = 0, year_sd = 0) {
  out <- env_steps(env, n_steps)
  m <- nrow(out)
  steps <- seq_len(ncol(out))
  check_number(period, "period", lower = 0, open = TRUE)
  if (!is.numeric(amplitude) || !length(amplitude) %in% c(1, m) || !all(is.finite(amplitude))) {
    stop("`amplitude` must be one finite number or one per site (", m, ")", call. = FALSE)
  }
  check_number(phase, "phase")
  check_number(year_sd, "year_sd", lower = 0)
  cycle <- outer(rep_len(amplitude, m), sin(2 * pi * (steps - phase) / period))
  out <- out + cycle
  if (year_sd > 0) {
    year <- ceiling(steps / period)
    anomaly <- stats::rnorm(max(year), 0, year_sd)
    out <- out + rep(anomaly[year], each = m)
  }
  out
}

#' Pulse disturbance of the environment
#'
#' Sets the environment of some sites to a new value from step `start` for `duration` steps,
#' after which it returns to its previous course, as the pulse disturbances of
#' Hajian-Forooshani & Chase (2026). Disturbance then acts on growth rates through the niche
#' of each species, not by removing individuals as `disturbance` in [simulate_metacommunity()]
#' does. Their disturbed environment (0.1) lies below the optima of all species (0.5 to 1), so
#' that no species grows well in disturbed sites; a value inside the range of optima instead
#' shifts which species are favoured, and `duration = Inf` makes the change permanent.
#'
#' @param env Environment of each site (a vector, e.g. the `env` of an `mc_landscape`), or a
#'   matrix (sites x steps) with the environment at each step, e.g. from [env_over_time()].
#' @param n_steps Number of time steps: one whole number, 1 or more. Taken from `env` when it
#'   is a matrix.
#' @param sites Disturbed sites: indices or a logical vector with one value per site. With no
#'   site (e.g. `integer(0)`), `env` is returned unchanged, as an undisturbed control.
#' @param start First disturbed step: one whole number, 1 or more.
#' @param duration Number of disturbed steps of each pulse: one whole number, 1 or more, or
#'   `Inf` for a permanent change.
#' @param every Optional period of recurrent pulses (one whole number, 1 or more): a new pulse
#'   starts every `every` steps
#'   after `start` until the last step, so that the frequency of disturbance is `1 / every`.
#'   `NULL` (the default) for a single pulse. Needs a finite `duration` shorter than `every`.
#' @param value Environment of the disturbed sites: one value, or one per disturbed site.
#' @return Matrix with one row per site and one column per time step, for the `env_time`
#'   argument of [simulate_metacommunity()].
#' @references Hajian-Forooshani, Z., & Chase, J. M. (2026). Disturbance interacts with
#'   dispersal and niche breadth to shape scale-dependent diversity change in metacommunities.
#'   *The American Naturalist*, 207, 838--850. \doi{10.1086/739955}
#' @family landscapes
#' @examples
#' set.seed(1)
#' env <- runif(6, 0.5, 1)
#' e <- env_pulse(env, n_steps = 30, sites = 1:3, start = 10, duration = 10, value = 0.1)
#' matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
#' # recurrent pulses of 3 steps every 10 steps
#' e <- env_pulse(env, n_steps = 60, sites = 1:3, start = 5, duration = 3, value = 0.1,
#'                every = 10)
#' which(e[1, ] == 0.1)
#' @export
env_pulse <- function(env, n_steps = NULL, sites, start, duration = Inf, value,
                      every = NULL) {
  out <- env_steps(env, n_steps)
  m <- nrow(out)
  if (is.logical(sites)) {
    if (length(sites) != m) stop("a logical `sites` needs one value per site", call. = FALSE)
    sites <- which(sites)
  }
  if (!all(sites %in% seq_len(m)) || anyDuplicated(sites)) {
    stop("`sites` must be distinct indices of sites", call. = FALSE)
  }
  check_count(start, "start")
  if (start > ncol(out)) {
    stop("`start` is after the last step (", ncol(out), ")", call. = FALSE)
  }
  # Inf makes the change permanent
  if (!identical(duration, Inf)) check_count(duration, "duration")
  if (!is.null(every)) {
    check_count(every, "every")
    if (duration >= every) stop("`every` must be longer than `duration`", call. = FALSE)
  }
  n_values <- if (length(sites) > 0) c(1, length(sites)) else 1
  if (!is.numeric(value) || anyNA(value) || !length(value) %in% n_values) {
    stop("`value` must be one number, or one per disturbed site, without NA", call. = FALSE)
  }
  if (length(sites) == 0) return(out)
  starts <- if (is.null(every)) start else seq(start, ncol(out), by = every)
  steps <- unlist(lapply(starts, function(s) seq(s, min(ncol(out), s + duration - 1))))
  out[sites, steps] <- value
  out
}

new_landscape <- function(env, coords, dist, dispersal, torus = NA) {
  structure(list(env = env, coords = coords, dist = dist, dispersal = dispersal,
                 torus = torus), class = "mc_landscape")
}

#' @export
print.mc_landscape <- function(x, ...) {
  type <- if (is.null(x$dispersal)) "spatially implicit, uniform dispersal" else
    if (!is.null(x$network)) paste0("river network (", sum(x$network$headwater),
                                    " headwaters)") else
    if (isTRUE(x$matrix)) "dispersal matrix" else
    paste0(if (isTRUE(x$points)) "points" else "grid", if (isTRUE(x$torus)) " (torus)" else "",
           ", exponential dispersal kernel")
  cat("<mc_landscape>", length(x$env), "sites;", type, "\n")
  cat("environment: ", paste(format(range(x$env), digits = 3), collapse = " to "), "\n")
  invisible(x)
}

#' great-circle distances in kilometres between points given as longitude and latitude
#' @noRd
great_circle <- function(coords) {
  radians <- as.matrix(coords) * pi / 180
  lon <- radians[, 1]
  lat <- radians[, 2]
  h <- sin(outer(lat, lat, "-") / 2)^2 +
    outer(cos(lat), cos(lat)) * sin(outer(lon, lon, "-") / 2)^2
  matrix(2 * 6371 * asin(sqrt(pmin(h, 1))), nrow(radians))
}
