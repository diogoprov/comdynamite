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
#' @return Data frame with columns `optimum`, `sigma` and `r_max`, one row per
#'   species, with species names as row names.
#' @examples
#' set.seed(1)
#' tr <- ape::rcoal(10)
#' species_traits(tr, lambda = 1, sigma_signal = TRUE)
#' species_traits(5)
#' @export
species_traits <- function(species, lambda = 1, optima_range = c(0, 100),
                           uniform_optima = FALSE, sigma_mean = 10, sigma_sd = 10,
                           rmax_range = c(2, 8), sigma_signal = FALSE,
                           rmax_signal = FALSE) {
  if (inherits(species, "phylo")) {
    k <- ape::Ntip(species)
    sp <- species$tip.label
    u <- trait_bm(species, lambda)
  } else {
    stopifnot(length(species) == 1, species >= 1)
    if (sigma_signal || rmax_signal)
      stop("phylogenetic signal needs a tree: pass a 'phylo' object as `species`")
    k <- as.integer(species)
    sp <- paste0("sp", seq_len(k))
    u <- stats::rnorm(k)
  }
  u <- if (uniform_optima) rank_to_marginal(u, seq(optima_range[1], optima_range[2], length.out = k))
       else rescale(u, optima_range)
  sigma <- numeric(0)
  while (length(sigma) < k) {
    s <- stats::rnorm(2 * k, sigma_mean, sigma_sd)
    sigma <- c(sigma, s[s > 0])
  }
  sigma <- sigma[seq_len(k)]
  r_max <- stats::runif(k, rmax_range[1], rmax_range[2])
  if (sigma_signal) sigma <- rank_to_marginal(trait_bm(species, 1), sigma)
  if (rmax_signal) r_max <- rank_to_marginal(trait_bm(species, 1), r_max)
  data.frame(optimum = unname(u), sigma = unname(sigma), r_max = unname(r_max),
             row.names = sp)
}

#' Trait evolution by Brownian motion with Pagel's lambda
#'
#' Draws from a multivariate normal distribution with covariance
#' `lambda * C + (1 - lambda) * diag(C)`, where `C` is the phylogenetic
#' covariance matrix of `tree`.
#'
#' @param tree A `phylo` tree with branch lengths.
#' @param lambda Pagel's lambda, in \[0, 1\].
#' @return Named numeric vector, one value per tip.
#' @examples
#' set.seed(1)
#' trait_bm(ape::rcoal(5), lambda = 0.5)
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
#' @return Matrix with one row per site and one column per species.
#' @examples
#' set.seed(1)
#' g <- growth_rates(seq(0, 100, 10), species_traits(3))
#' round(g, 2)
#' @export
growth_rates <- function(env, traits) {
  stopifnot(all(c("optimum", "sigma", "r_max") %in% names(traits)))
  r <- exp(-sweep(outer(env, traits$optimum, "-")^2, 2, 2 * traits$sigma^2, "/"))
  r <- sweep(r, 2, traits$r_max, "*")
  dimnames(r) <- list(names(env), rownames(traits))
  r
}

rank_to_marginal <- function(x, target) {
  stats::setNames(sort(target)[rank(x, ties.method = "random")], names(x))
}

rescale <- function(x, range) {
  if (length(x) == 1 || diff(base::range(x)) == 0) return(rep(mean(range), length(x)))
  range[1] + (x - min(x)) / diff(base::range(x)) * diff(range)
}
