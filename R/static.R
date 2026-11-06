#' Simulate a static metacommunity
#'
#' Draws species abundances directly from unimodal responses to an
#' environmental gradient, without population dynamics, competition or
#' dispersal. The expected abundance of species j at site i is
#' `r_max[j] * exp(-(env[i] - optimum[j])^2 / (2 * sigma[j]^2))`, the Gaussian
#' response used in the simulations of Peres-Neto et al. (2012), and observed
#' abundances are Poisson counts with that mean. Because nothing but the
#' gradient structures the data, this is a baseline for
#' [simulate_metacommunity()] and a fast generator for evaluating methods.
#'
#' @param traits Data frame from [species_traits()], or any data frame with
#'   columns `optimum`, `sigma` and `r_max`. Here `r_max` is the maximum
#'   expected abundance, not a growth rate, so set `rmax_range` accordingly.
#' @param env Environmental value of each site.
#' @param incidence Return presences and absences instead of abundances?
#' @return List of class `mc_sim` with `abundance`, a matrix with one row per
#'   site and one column per species, plus `env` and `traits`. Sites and species
#'   with no individuals are kept, so remove them if a method requires it.
#' @references
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects
#' of spatial contingency and environmental filtering on metacommunity
#' phylogenetics. *Ecology*, 93, S14--S30. \doi{10.1890/11-0494.1}
#' @examples
#' set.seed(1)
#' tr <- ape::rcoal(20)
#' traits <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
#'                          rmax_range = c(0, 30))
#' sim <- simulate_static(traits, env = runif(100, 0, 100))
#' sim
#' @export
simulate_static <- function(traits, env, incidence = FALSE) {
  stopifnot(is.numeric(env), length(env) >= 1, all(traits$r_max >= 0))
  expected <- growth_rates(env, traits)
  abundance <- matrix(stats::rpois(length(expected), expected), nrow(expected),
                      ncol(expected), dimnames = list(
                        rownames(expected) %||% paste0("site", seq_along(env)),
                        rownames(traits)))
  if (incidence) abundance[] <- (abundance > 0) * 1L
  structure(list(abundance = abundance, history = NULL, env = env, traits = traits),
            class = "mc_sim")
}
