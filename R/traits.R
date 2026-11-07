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
#' @param species A `phylo` tree, or the number of species for traits without
#'   phylogenetic structure.
#' @param lambda Pagel's lambda of the optima (0 = no signal, 1 = Brownian
#'   motion). Ignored when `species` is a number.
#' @param optima_range Range of the optima.
#' @param uniform_optima Give the optima evenly spaced values by rank?
#' @param sigma_mean,sigma_sd Mean and standard deviation of the normal
#'   distribution (truncated at zero) of niche breadths.
#' @param rmax_range Range of the uniform distribution of maximum growth
#'   rates.
#' @param sigma_signal,rmax_signal Should niche breadth or maximum growth rate
#'   carry phylogenetic signal (Brownian motion)?
#' @param optima,sigma,r_max Optional values of each trait, named by species or in the order
#'   of the species, as a vector or a one-column matrix or data frame. They replace the values
#'   drawn for that trait.
#' @param rescale_optima Rescale external optima to `optima_range`?
#' @return Data frame with columns `optimum`, `sigma` and `r_max`, one row per
#'   species, with species names as row names.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(10)
#' species_traits(tr, lambda = 1, sigma_signal = TRUE)
#' species_traits(5)
#' # optima from another model, here an Ornstein-Uhlenbeck process
#' ou <- traits_ou(tr, alpha = 1, theta = 0)
#' species_traits(tr, optima = ou[, 1])
#' @export
species_traits <- function(species, lambda = 1, optima_range = c(0, 100),
                           uniform_optima = FALSE, sigma_mean = 10, sigma_sd = 10,
                           rmax_range = c(2, 8), sigma_signal = FALSE,
                           rmax_signal = FALSE, optima = NULL, sigma = NULL, r_max = NULL,
                           rescale_optima = TRUE) {
  if (inherits(species, "phylo")) {
    k <- ape::Ntip(species)
    sp <- species$tip.label
  } else {
    stopifnot(length(species) == 1, species >= 1)
    if (sigma_signal || rmax_signal)
      stop("phylogenetic signal needs a tree: pass a 'phylo' object as `species`")
    k <- as.integer(species)
    sp <- paste0("sp", seq_len(k))
  }
  given <- list(optima = external_trait(optima, sp, "optima"),
                sigma = external_trait(sigma, sp, "sigma"),
                r_max = external_trait(r_max, sp, "r_max"))
  if (!is.null(given$sigma) && any(given$sigma <= 0)) {
    stop("`sigma` must be positive", call. = FALSE)
  }
  if (!is.null(given$r_max) && any(given$r_max < 0)) {
    stop("`r_max` cannot be negative", call. = FALSE)
  }
  # values are drawn even when given, so the other traits do not depend on what was given
  u <- if (inherits(species, "phylo")) trait_bm(species, lambda) else stats::rnorm(k)
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
  out
}

#' Trait evolution by Brownian motion with Pagel's lambda
#'
#' Draws from a multivariate normal distribution with covariance
#' `lambda * C + (1 - lambda) * diag(C)`, where `C` is the phylogenetic
#' covariance matrix of `tree`.
#'
#' @param tree A `phylo` tree with branch lengths.
#' @param lambda Pagel's lambda, in \[0, 1\].
#' @return Named numeric vector with the simulated trait value of each tip.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' trait_bm(treeshapesim::sim_tree_shape(5), lambda = 0.5)
#' @export
trait_bm <- function(tree, lambda = 1) {
  stopifnot(lambda >= 0, lambda <= 1)
  cv <- ape::vcv.phylo(tree)
  cv <- lambda * cv + (1 - lambda) * diag(diag(cv))
  x <- drop(crossprod(chol(cv), stats::rnorm(ape::Ntip(tree))))
  stats::setNames(x, tree$tip.label)
}

#' Growth rates of each species at each site
#'
#' Density-independent growth rate of species j at site i,
#' `r_max[j] * exp(-(E[i] - optimum[j])^2 / (2 * sigma[j]^2))`: a Gaussian
#' response whose standard deviation is the niche breadth `sigma`. This
#' `sigma` equals `sqrt(2)` times the niche breadth parameter of Thompson et
#' al. (2020).
#'
#' @param env Environment of each site.
#' @param traits Data frame from [species_traits()], or any data frame with
#'   columns `optimum`, `sigma` and `r_max`.
#' @return Numeric matrix with one row per site and one column per species:
#'   the density-independent growth rate of each species at each site.
#' @examples
#' set.seed(1)
#' g <- growth_rates(seq(0, 100, 10), species_traits(3))
#' round(g, 2)
#' @export
growth_rates <- function(env, traits) {
  traits <- as_traits(traits)
  r <- exp(-sweep(outer(env, traits$optimum, "-")^2, 2, 2 * traits$sigma^2, "/"))
  r <- sweep(r, 2, traits$r_max, "*")
  dimnames(r) <- list(names(env), rownames(traits))
  r
}

# traits of species_traits() as a data frame with names of species, also from a matrix or a
# tibble, which has no row names
# a trait given from outside, in the order of `species`; NULL when not given
external_trait <- function(x, species, name) {
  if (is.null(x)) return(NULL)
  if (!is.null(dim(x))) {
    if (ncol(x) != 1) stop("`", name, "` needs one column", call. = FALSE)
    labels <- rownames(x)
    if (is.data.frame(x) && .row_names_info(x) < 0) labels <- NULL
    x <- stats::setNames(if (is.data.frame(x)) x[[1]] else x[, 1], labels)
  }
  if (!is.numeric(x) || anyNA(x)) stop("`", name, "` must be numbers without NA", call. = FALSE)
  if (is.null(names(x))) {
    if (length(x) != length(species)) {
      stop("`", name, "` has ", length(x), " values and there are ", length(species), " species",
           call. = FALSE)
    }
    return(unname(x))
  }
  missing_species <- setdiff(species, names(x))
  if (length(missing_species) > 0) {
    stop("`", name, "` has no value for ", species_list(missing_species), call. = FALSE)
  }
  unname(x[species])
}

as_traits <- function(traits) {
  traits <- as.data.frame(traits)
  missing_traits <- setdiff(c("optimum", "sigma", "r_max"), names(traits))
  if (length(missing_traits) > 0) {
    stop("`traits` needs the columns of `species_traits()`; missing: ",
         paste(missing_traits, collapse = ", "), call. = FALSE)
  }
  if (.row_names_info(traits) < 0) rownames(traits) <- paste0("sp", seq_len(nrow(traits)))
  traits
}

rank_to_marginal <- function(x, target) {
  stats::setNames(sort(target)[rank(x, ties.method = "random")], names(x))
}

rescale <- function(x, range) {
  if (length(x) == 1 || diff(base::range(x)) == 0) return(rep(mean(range), length(x)))
  range[1] + (x - min(x)) / diff(base::range(x)) * diff(range)
}
