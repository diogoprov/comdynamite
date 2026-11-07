#' Spatially implicit landscape
#'
#' Sites without coordinates, with environmental values drawn from a uniform
#' distribution. Dispersal among sites is uniform: emigrants can reach any
#' site, so sites are exchangeable and permuting them is a valid null model.
#'
#' @param n_sites Number of sites.
#' @param env_range Range of the uniform distribution of the environment.
#' @return An `mc_landscape` object: a list with `env` (one value per site),
#'   and `coords`, `dist` and `dispersal` set to `NULL`.
#' @seealso [landscape_grid()] for a spatially explicit landscape.
#' @examples
#' set.seed(1)
#' landscape_implicit(20)
#' @export
landscape_implicit <- function(n_sites, env_range = c(0, 100)) {
  stopifnot(n_sites >= 2, length(env_range) == 2, env_range[1] < env_range[2])
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
#' @param side Number of cells per side.
#' @param env_bandwidth Bandwidth (in cells) of the Gaussian kernel that
#'   smooths the environment. Larger values give smoother gradients.
#' @param disp_scale Scale (in cells) of the exponential dispersal kernel.
#' @param torus Wrap the grid as a torus?
#' @param env_range Range of the uniform marginal distribution of the
#'   environment.
#' @return An `mc_landscape` object: a list with `env`, `coords` (data frame
#'   with columns `x` and `y`), `dist` (distances among sites, toroidal if
#'   `torus = TRUE`) and `dispersal` (row-stochastic matrix with a zero
#'   diagonal).
#' @seealso [torus_permutations()] for site permutations that keep the spatial
#'   structure of a torus.
#' @examples
#' set.seed(1)
#' ls <- landscape_grid(side = 8)
#' image(matrix(ls$env, 8, 8), asp = 1)
#' @export
landscape_grid <- function(side = 10, env_bandwidth = 2, disp_scale = 2,
                           torus = TRUE, env_range = c(0, 100)) {
  stopifnot(side >= 2, env_bandwidth > 0, disp_scale > 0)
  xy <- expand.grid(x = seq_len(side) - 1, y = seq_len(side) - 1)
  dx <- abs(outer(xy$x, xy$x, "-"))
  dy <- abs(outer(xy$y, xy$y, "-"))
  if (torus) {
    dx <- pmin(dx, side - dx)
    dy <- pmin(dy, side - dy)
  }
  d <- sqrt(dx^2 + dy^2)
  z <- drop(exp(-d^2 / (2 * env_bandwidth^2)) %*% stats::rnorm(side^2))
  env <- rank_to_marginal(z, stats::runif(side^2, env_range[1], env_range[2]))
  new_landscape(env = unname(env), coords = xy, dist = d,
                dispersal = dispersal_kernel(d, disp_scale), torus = torus)
}

#' Exponential dispersal kernel
#'
#' The probability that an emigrant from site s reaches site t is proportional
#' to exp(-d / `scale`), where d is the distance between the sites. Emigrants
#' do not stay at their own site.
#'
#' @param dist Square matrix of distances among sites.
#' @param scale Scale of the kernel, in the units of `dist`.
#' @return Row-stochastic matrix with a zero diagonal: row s gives the
#'   probabilities that an emigrant from site s reaches each site.
#' @examples
#' d <- as.matrix(dist(cbind(1:5, 0)))
#' round(dispersal_kernel(d, scale = 1), 2)
#' @export
dispersal_kernel <- function(dist, scale) {
  stopifnot(is.matrix(dist), nrow(dist) == ncol(dist), scale > 0)
  k <- exp(-dist / scale)
  diag(k) <- 0
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
#'   rectangular grid, the first for the coordinate that runs fastest.
#' @param symmetries Combine the translations with the reflections and rotations of the grid?
#' @param n_perm Number of permutations to keep, drawn without replacement. The default keeps
#'   all of them.
#' @return Integer matrix with one permutation of the site indices per row, excluding the
#'   identity (`8 * side^2 - 1` rows on a square grid with `side > 2`).
#' @references Harms, K. E., Condit, R., Hubbell, S. P., & Foster, R. B. (2001). Habitat
#'   associations of trees and shrubs in a 50-ha neotropical forest plot. *Journal of
#'   Ecology*, 89, 947--959.
#' @examples
#' p <- torus_permutations(4)
#' dim(p)
#' dim(torus_permutations(c(6, 4), symmetries = FALSE))
#' @export
torus_permutations <- function(side, symmetries = TRUE, n_perm = NULL) {
  stopifnot(is.numeric(side), length(side) %in% 1:2, side >= 2, side == round(side))
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
  stopifnot(n_perm >= 1)
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
#' @param n_steps Number of time steps.
#' @param rho Lag-one autocorrelation of the deviations, in \[0, 1).
#' @param sd Stationary standard deviation of the deviations.
#' @return Matrix with one row per site and one column per time step, for the
#'   `env_time` argument of [simulate_metacommunity()].
#' @examples
#' set.seed(1)
#' e <- env_over_time(c(20, 50, 80), n_steps = 100, rho = 0.95, sd = 5)
#' matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
#' @export
env_over_time <- function(env, n_steps, rho = 0.9, sd = 5) {
  stopifnot(n_steps >= 1, rho >= 0, rho < 1, sd >= 0)
  m <- length(env)
  dev <- matrix(0, m, n_steps)
  dev[, 1] <- stats::rnorm(m, 0, sd)
  innov <- sd * sqrt(1 - rho^2)
  for (t in seq_len(n_steps)[-1]) dev[, t] <- rho * dev[, t - 1] + stats::rnorm(m, 0, innov)
  env + dev
}

new_landscape <- function(env, coords, dist, dispersal, torus = NA) {
  structure(list(env = env, coords = coords, dist = dist, dispersal = dispersal,
                 torus = torus), class = "mc_landscape")
}

#' @export
print.mc_landscape <- function(x, ...) {
  type <- if (is.null(x$dispersal)) "spatially implicit, uniform dispersal" else
    paste0("grid", if (isTRUE(x$torus)) " (torus)" else "", ", exponential dispersal kernel")
  cat("<mc_landscape>", length(x$env), "sites;", type, "\n")
  cat("environment: ", paste(format(range(x$env), digits = 3), collapse = " to "), "\n")
  invisible(x)
}
