#' Species traits, with or without phylogenetic signal
#'
#' Environmental optima evolve by Brownian motion with Pagel's lambda and are
#' rescaled to `optima_range` (or, with `uniform_optima = TRUE`, replaced by
#' ranks with evenly spaced values on that range, so that only their
#' arrangement on the tree depends on the tree). Niche breadths are drawn from
#' a normal distribution truncated at zero, and maximum growth rates from a
#' uniform distribution. Either can carry phylogenetic signal: the values are
#' then assigned to species by the ranks of a Brownian-motion trait, which
#' keeps their marginal distribution.
#'
#' Any of the three traits can come from elsewhere instead, such as optima simulated by
#' mvMORPH or measured in the field, through `optima`, `sigma` and `r_max`. The other traits are
#' still drawn, and with the same seed they are the same as without the external ones. External
#' optima are rescaled to `optima_range`, or ranked with `uniform_optima`, unless
#' `rescale_optima = FALSE`.
#'
#' @param species A `phylo` tree, or the number of species (one whole number, 1 or more) for
#'   traits without phylogenetic structure.
#' @param lambda Pagel's lambda of the optima (0 = no signal, 1 = Brownian
#'   motion): one number between 0 and 1. Ignored when `species` is a number.
#' @param optima_range Range of the optima: two finite numbers, the first not larger than the
#'   second.
#' @param uniform_optima Give the optima evenly spaced values by rank? `TRUE` or `FALSE`.
#' @param sigma_mean,sigma_sd Mean and standard deviation of the normal
#'   distribution (truncated at zero) of niche breadths: one finite number each, `sigma_sd` 0
#'   or more.
#' @param rmax_range Range of the uniform distribution of maximum growth
#'   rates: two finite numbers, 0 or more, the first not larger than the second.
#' @param sigma_signal,rmax_signal Should niche breadth or maximum growth rate
#'   carry phylogenetic signal (Brownian motion)? `TRUE` or `FALSE` each.
#' @param optima,sigma,r_max Optional values of each trait, named by species or in the order
#'   of the species, as a vector or a one-column matrix or data frame. They replace the values
#'   drawn for that trait.
#' @param rescale_optima Rescale external optima to `optima_range`? `TRUE` or `FALSE`.
#' @param nodes Also give the optima at the internal nodes of the tree, as [trait_bm()] does?
#'   `TRUE` or `FALSE`. They are drawn after all other values, so the traits are the same as with
#'   `nodes = FALSE` and the same seed, and rescaled as the optima of the tips, so they may fall
#'   outside `optima_range`. Not available with external optima or `uniform_optima`.
#' @return Data frame with columns `optimum`, `sigma` and `r_max`, one row per
#'   species, with species names as row names. With a tree, the attribute `tree` holds it, and
#'   with `nodes = TRUE` the attribute `nodes` holds the optima at its internal nodes.
#' @family species traits
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(10)
#' species_traits(tr, lambda = 1, sigma_signal = TRUE)
#' species_traits(5)
#' # optima from another model, here an Ornstein-Uhlenbeck process
#' ou <- trait_ou(tr, alpha = 1, theta = 0)
#' species_traits(tr, optima = ou[, 1])
#' @export
species_traits <- function(species, lambda = 1, optima_range = c(0, 100),
                           uniform_optima = FALSE, sigma_mean = 10, sigma_sd = 10,
                           rmax_range = c(2, 8), sigma_signal = FALSE,
                           rmax_signal = FALSE, optima = NULL, sigma = NULL, r_max = NULL,
                           rescale_optima = TRUE, nodes = FALSE) {
  check_number(lambda, "lambda", 0, 1)
  check_range_pair(optima_range, "optima_range")
  check_range_pair(rmax_range, "rmax_range", lower = 0)
  check_number(sigma_mean, "sigma_mean")
  check_number(sigma_sd, "sigma_sd", lower = 0)
  if (sigma_sd == 0 && sigma_mean <= 0) {
    stop("with `sigma_sd = 0`, `sigma_mean` must be greater than 0", call. = FALSE)
  }
  for (flag in c("uniform_optima", "sigma_signal", "rmax_signal", "rescale_optima", "nodes")) {
    check_flag(get(flag), flag)
  }
  if (inherits(species, "phylo")) {
    k <- ape::Ntip(species)
    sp <- species$tip.label
  } else {
    if (!is.numeric(species) || length(species) != 1 || !is.finite(species) || species < 1 ||
          species != round(species)) {
      stop("`species` must be a `phylo` tree or the number of species (one whole number, 1 or ",
           "more)", call. = FALSE)
    }
    if (sigma_signal || rmax_signal)
      stop("phylogenetic signal needs a tree: pass a 'phylo' object as `species`", call. = FALSE)
    k <- as.integer(species)
    sp <- paste0("sp", seq_len(k))
  }
  if (nodes && (!inherits(species, "phylo") || !is.null(optima) || uniform_optima)) {
    stop("`nodes` needs a tree and optima drawn by `species_traits()` without `uniform_optima`",
         call. = FALSE)
  }
  given <- list(optima = external_trait(optima, sp, "optima"),
                sigma = external_trait(sigma, sp, "sigma", finite = FALSE),
                r_max = external_trait(r_max, sp, "r_max"))
  if (!is.null(given$sigma) && any(given$sigma <= 0)) {
    stop("`sigma` must be positive", call. = FALSE)
  }
  if (!is.null(given$r_max) && any(given$r_max < 0)) {
    stop("`r_max` cannot be negative", call. = FALSE)
  }
  # values are drawn even when given, so the other traits do not depend on what was given
  u <- if (inherits(species, "phylo")) trait_bm(species, lambda) else stats::rnorm(k)
  drawn <- u
  if (!is.null(given$optima)) u <- given$optima
  if (is.null(given$optima) || rescale_optima) {
    u <- if (uniform_optima) {
      rank_to_marginal(u, seq(optima_range[1], optima_range[2], length.out = k))
    } else {
      rescale(u, optima_range)
    }
  }
  breadth <- numeric(0)
  while (length(breadth) < k) {
    s <- stats::rnorm(2 * k, sigma_mean, sigma_sd)
    breadth <- c(breadth, s[s > 0])
  }
  breadth <- breadth[seq_len(k)]
  growth <- stats::runif(k, rmax_range[1], rmax_range[2])
  if (sigma_signal) breadth <- rank_to_marginal(trait_bm(species, 1), breadth)
  if (rmax_signal) growth <- rank_to_marginal(trait_bm(species, 1), growth)
  sigma <- given$sigma %||% breadth
  r_max <- given$r_max %||% growth
  out <- data.frame(optimum = unname(u), sigma = unname(sigma), r_max = unname(r_max),
                    row.names = sp)
  if (inherits(species, "phylo")) attr(out, "tree") <- species
  if (nodes) {
    at_nodes <- bm_nodes(species, lambda, drawn)
    attr(out, "nodes") <- optima_range[1] + (at_nodes - min(drawn)) / diff(range(drawn)) *
      diff(optima_range)
  }
  out
}

#' Trait evolution by Brownian motion with Pagel's lambda
#'
#' Draws from a multivariate normal distribution with covariance
#' \eqn{\lambda C + (1 - \lambda)\, \mathrm{diag}(C)}{lambda * C + (1 - lambda) * diag(C)},
#' where \eqn{C} is the phylogenetic covariance matrix of `tree`.
#'
#' @param tree A `phylo` tree with branch lengths.
#' @param lambda Pagel's lambda: one number between 0 and 1.
#' @param nodes Also draw the values at the internal nodes? `TRUE` or `FALSE`. They are drawn after
#'   the tips, given them, so the tips are the same as with `nodes = FALSE` and the same seed. Under
#'   lambda, the nodes are those of the tree whose internal branches are multiplied by lambda and
#'   whose terminal branches are lengthened to keep the tips where they were.
#' @return Named numeric vector with the simulated trait value of each tip. With `nodes = TRUE`,
#'   the attribute `nodes` holds the values at the internal nodes, named by their numbers in
#'   `tree$edge`, the root first.
#' @family species traits
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' trait_bm(treeshapesim::sim_tree_shape(5), lambda = 0.5)
#' @export
trait_bm <- function(tree, lambda = 1, nodes = FALSE) {
  check_tree(tree)
  check_number(lambda, "lambda", 0, 1)
  check_flag(nodes, "nodes")
  cv <- ape::vcv.phylo(tree)
  cv <- lambda * cv + (1 - lambda) * diag(diag(cv))
  x <- drop(crossprod(chol(cv), stats::rnorm(ape::Ntip(tree))))
  x <- stats::setNames(x, tree$tip.label)
  if (nodes) attr(x, "nodes") <- bm_nodes(tree, lambda, x)
  x
}

#' Growth rates of each species at each site
#'
#' Density-independent growth rate of species \eqn{j} at site \eqn{i},
#' \eqn{r_j \exp(-(E_i - \theta_j)^2 / (2 \sigma_j^2))}{r_max[j] * exp(-(E[i] - optimum[j])^2 /
#' (2 * sigma[j]^2))}, with \eqn{r_j} = `r_max`, \eqn{\theta_j}{theta_j} = `optimum` and
#' \eqn{\sigma_j}{sigma_j} = `sigma`: a Gaussian response whose standard deviation is the niche
#' breadth `sigma`. This `sigma` equals \eqn{\sqrt{2}}{sqrt(2)} times the niche breadth
#' parameter of Thompson et al. (2020).
#'
#' With several environmental axes, each species has a response to each axis and the growth
#' rate is the product of the responses, as in [simulate_static()]: a Gaussian niche in
#' several dimensions whose axes are independent (no correlation between the responses to
#' different axes), with the optimum and breadth of each axis in its own set of traits. The
#' maxima multiply too, so give the traits of all but one axis `r_max` of 1, for example with
#' `species_traits(rmax_range = c(1, 1))`.
#'
#' A list of traits passed through [niche_correlation()] gives the responses to the axes a
#' correlation \eqn{R}{R}: the growth rate is then \eqn{\prod_a r_a \exp(-d^\top R^{-1} d /
#' 2)}{prod(r_max) exp(-t(d) solve(R) d / 2)}, with \eqn{d_a = (E_a - \theta_a)/\sigma_a}{d_a =
#' (E_a - theta_a) / sigma_a}, a multivariate Gaussian niche whose covariance has the breadths
#' \eqn{\sigma_a}{sigma_a} on the diagonal. A positive correlation favours sites where both
#' axes depart from the optima in the same direction. With \eqn{R}{R} the identity, it is the
#' product above.
#'
#' @param env Environment of each site; for several axes, a matrix, data frame or list with
#'   one axis per column or element, in the order of `traits`.
#' @param traits Data frame from [species_traits()], or any data frame with
#'   columns `optimum`, `sigma` and `r_max`; for several axes, a list of them with the same
#'   species.
#' @return Numeric matrix with one row per site and one column per species:
#'   the density-independent growth rate of each species at each site.
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' g <- growth_rates(seq(0, 100, 10), species_traits(3))
#' round(g, 2)
#' @export
growth_rates <- function(env, traits) {
  if (is.list(traits) && !is.data.frame(traits)) {
    traits_in <- traits
    traits <- lapply(traits, as_traits)
    axes <- if (is.matrix(env)) as.data.frame(env) else as.list(env)
    if (length(axes) != length(traits)) {
      stop("`env` has ", length(axes), " axes and `traits` ", length(traits), " sets of traits",
           call. = FALSE)
    }
    species <- rownames(traits[[1]])
    if (!all(vapply(traits, function(x) setequal(rownames(x), species), logical(1)))) {
      stop("growth_rates(): sets of traits have different species", call. = FALSE)
    }
    sites <- if (is.matrix(env) || is.data.frame(env)) rownames(env)
    if (is.data.frame(env) && .row_names_info(env) < 0) sites <- NULL
    traits <- lapply(traits, function(x) x[species, , drop = FALSE])
    cor <- attr(traits_in, "niche_cor")
    r <- if (is.null(cor)) {
      Reduce(`*`, Map(function(e, x) growth_rates(unname(e), x), unname(as.list(axes)), traits))
    } else {
      correlated_growth(axes, traits, cor, species)
    }
    rownames(r) <- sites
    return(r)
  }
  if (is.data.frame(env) || is.matrix(env)) {
    if (ncol(env) != 1) stop("several axes in `env` need a list of traits", call. = FALSE)
    env <- stats::setNames(if (is.data.frame(env)) env[[1]] else env[, 1], rownames(env))
  }
  traits <- as_traits(traits)
  r <- exp(-sweep(outer(env, traits$optimum, "-")^2, 2, 2 * traits$sigma^2, "/"))
  r <- sweep(r, 2, traits$r_max, "*")
  dimnames(r) <- list(names(env), rownames(traits))
  r
}

#' multivariate Gaussian response: standardised departures d (sites x axes) for each species,
#' exp(-t(d) solve(R) d / 2), times the product of the maxima
#' @noRd
correlated_growth <- function(axes, traits, cor, species) {
  q <- length(traits)
  env <- vapply(unname(as.list(axes)), function(e) as.numeric(unname(e)), numeric(length(axes[[1]])))
  env <- matrix(env, ncol = q)
  cor <- niche_cor_species(cor, species, q)
  r <- matrix(0, nrow(env), length(species))
  for (j in seq_along(species)) {
    theta <- vapply(traits, function(x) x$optimum[j], numeric(1))
    sigma <- vapply(traits, function(x) x$sigma[j], numeric(1))
    rmax <- prod(vapply(traits, function(x) x$r_max[j], numeric(1)))
    d <- sweep(sweep(env, 2, theta, "-"), 2, sigma, "/")
    d[, !is.finite(sigma)] <- 0
    m <- rowSums((d %*% solve(cor[[j]])) * d)
    r[, j] <- rmax * exp(-m / 2)
  }
  colnames(r) <- species
  r
}

#' the correlation matrix of each species (a list in the order of `species`) from the attribute
#' of niche_correlation(); NULL without it
#' @noRd
niche_cor_species <- function(cor, species, q) {
  if (is.null(cor)) return(NULL)
  if (is.list(cor)) {
    if (!all(species %in% names(cor))) {
      stop("the correlation of the niche lacks matrices for ",
           species_list(setdiff(species, names(cor))), call. = FALSE)
    }
    cor <- unname(cor[species])
  } else {
    cor <- rep(list(cor), length(species))
  }
  ok <- vapply(cor, function(m) is.matrix(m) && is.numeric(m) && identical(dim(m), c(q, q)) &&
                 !anyNA(m), logical(1))
  if (!all(ok)) {
    stop("the correlation of the niche needs ", q, " x ", q, " matrices, from ",
         "`niche_correlation()`", call. = FALSE)
  }
  lapply(cor, unname)
}

#' covariance Omega of the Gaussian niche of one species on its axes with a finite breadth
#' `s`: D R D with D = diag(sigma) when all are finite; otherwise the inverse of that block of
#' the precision, the limit of an infinite breadth, as in correlated_growth()
#' @noRd
niche_omega <- function(sigma, cor, s = which(is.finite(sigma))) {
  d <- diag(sigma[s], length(s))
  r <- if (length(s) == length(sigma)) cor else solve(solve(cor)[s, s, drop = FALSE])
  om <- d %*% r %*% d
  (om + t(om)) / 2
}

#' precision of the niche over all the axes, with zeros on those without a finite breadth
#' @noRd
niche_precision <- function(sigma, cor) {
  s <- which(is.finite(sigma))
  out <- matrix(0, length(sigma), length(sigma))
  if (length(s)) out[s, s] <- solve(niche_omega(sigma, cor, s))
  out
}

#' the matrix M of the overlap exp(-d' M d / 2) of two Gaussian niches with precisions a and b
#' whose optima differ by d: (Omega_a + Omega_b)^-1 when both covariances exist
#' @noRd
niche_overlap_matrix <- function(a, b) {
  u <- which(diag(a) > 0 | diag(b) > 0)
  m <- matrix(0, nrow(a), ncol(a))
  if (length(u)) {
    au <- a[u, u, drop = FALSE]
    bu <- b[u, u, drop = FALSE]
    mu <- au %*% solve(au + bu, bu)
    m[u, u] <- (mu + t(mu)) / 2
  }
  m
}

#' Correlation between the axes of a niche
#'
#' Attaches to a list of traits (one set per environmental axis, as for [growth_rates()] and
#' [simulate_metacommunity()]) a correlation between the responses of each species to the
#' axes, so that the niche is a multivariate Gaussian with breadths `sigma` on each axis and
#' this correlation. A species with a positive correlation between two axes grows best where
#' both depart from its optima in the same direction, as a species that tolerates warm water
#' only if it is also well oxygenated.
#'
#' With `evolution` in [simulate_metacommunity()] and in [simulate_evolving_metacommunity()],
#' the correlation enters the matrix of Gaussian selection, \eqn{\Omega = D R D}{Omega = D R
#' D} with \eqn{D}{D} the diagonal matrix of the breadths (see the section Several traits of
#' [evolution_control()]): growth, the response to selection, the selection of immigrants and
#' the overlap of niches (`niche_overlap`, and `competition_matrix(type = "trait")`) all use
#' it. New species take the correlation of their mother. An axis with an infinite breadth
#' exerts no selection, whatever its correlation; a breadth of 0 cannot be combined with a
#' correlation.
#'
#' @param traits List of data frames from [species_traits()], one per axis, with the same
#'   species.
#' @param cor A correlation matrix (axes x axes), shared by all species, or a list of them
#'   named by species. Each must be symmetric, with ones on the diagonal and positive
#'   definite.
#' @return `traits`, with the attribute `niche_cor`.
#' @family species traits
#' @examples
#' set.seed(1)
#' tr <- list(temperature = species_traits(2, optima = c(40, 60), sigma = c(10, 10),
#'                                         r_max = c(5, 5), rescale_optima = FALSE),
#'            oxygen = species_traits(2, optima = c(40, 60), sigma = c(10, 10),
#'                                    r_max = c(1, 1), rescale_optima = FALSE))
#' env <- cbind(temperature = c(30, 50), oxygen = c(30, 50))
#' growth_rates(env, tr)
#' growth_rates(env, niche_correlation(tr, matrix(c(1, 0.8, 0.8, 1), 2)))
#' @export
niche_correlation <- function(traits, cor) {
  if (!is.list(traits) || is.data.frame(traits) || length(traits) < 2) {
    stop("`traits` must be a list of at least two sets of traits, one per axis", call. = FALSE)
  }
  q <- length(traits)
  species <- rownames(as_traits(traits[[1]]))
  check <- function(m, what) {
    if (!is.matrix(m) || !identical(dim(m), c(q, q)) || !is.numeric(m) || anyNA(m)) {
      stop(what, " must be a ", q, " x ", q, " numeric matrix", call. = FALSE)
    }
    if (max(abs(m - t(m))) > 1e-10 || max(abs(diag(m) - 1)) > 1e-10) {
      stop(what, " must be symmetric with ones on the diagonal", call. = FALSE)
    }
    if (min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop(what, " must be positive definite", call. = FALSE)
    }
    unname(m)
  }
  if (is.list(cor)) {
    if (is.null(names(cor)) || !setequal(names(cor), species)) {
      stop("a list of `cor` needs one matrix per species, named by species", call. = FALSE)
    }
    cor <- lapply(stats::setNames(names(cor), names(cor)),
                  function(sp) check(cor[[sp]], paste0("`cor` of ", sp)))
  } else {
    cor <- check(cor, "`cor`")
  }
  attr(traits, "niche_cor") <- cor
  traits
}

#' traits of species_traits() as a data frame with names of species, also from a matrix or a
#' tibble, which has no row names
#' a trait given from outside, in the order of `species`; NULL when not given
#' @noRd
external_trait <- function(x, species, name, finite = TRUE) {
  if (is.null(x)) return(NULL)
  if (!is.null(dim(x))) {
    if (ncol(x) != 1) stop("`", name, "` needs one column", call. = FALSE)
    labels <- rownames(x)
    if (is.data.frame(x) && .row_names_info(x) < 0) labels <- NULL
    x <- stats::setNames(if (is.data.frame(x)) x[[1]] else x[, 1], labels)
  }
  if (!is.numeric(x) || anyNA(x)) stop("`", name, "` must be numbers without NA", call. = FALSE)
  if (finite && any(is.infinite(x))) stop("`", name, "` must be finite numbers", call. = FALSE)
  if (is.null(names(x))) {
    if (length(x) != length(species)) {
      stop("`", name, "` has ", length(x), " values and there are ", length(species), " species",
           call. = FALSE)
    }
    return(unname(x))
  }
  missing_species <- setdiff(species, names(x))
  if (length(missing_species) > 0) {
    stop("`", name, "` lacks values for ", species_list(missing_species), call. = FALSE)
  }
  unname(x[species])
}

as_traits <- function(traits) {
  traits <- as.data.frame(traits)
  missing_traits <- setdiff(c("optimum", "sigma", "r_max"), names(traits))
  if (length(missing_traits) > 0) {
    stop("`traits` needs the columns of `species_traits()`; missing: ",
         toString(missing_traits), call. = FALSE)
  }
  if (.row_names_info(traits) < 0) rownames(traits) <- paste0("sp", seq_len(nrow(traits)))
  traits
}

#' stops unless `x` is two finite numbers, the first not larger than the second, both `lower` or
#' more
#' @noRd
check_range_pair <- function(x, name, lower = -Inf) {
  if (!is.numeric(x) || length(x) != 2 || !all(is.finite(x)) || x[1] > x[2] || any(x < lower)) {
    stop("`", name, "` must be two finite numbers, the first not larger than the second",
         if (is.finite(lower)) paste0(", both ", lower, " or more"), call. = FALSE)
  }
  invisible(x)
}

rank_to_marginal <- function(x, target) {
  stats::setNames(sort(target)[rank(x, ties.method = "random")], names(x))
}

rescale <- function(x, range) {
  if (length(x) == 1 || diff(base::range(x)) == 0) return(rep(mean(range), length(x)))
  range[1] + (x - min(x)) / diff(base::range(x)) * diff(range)
}
