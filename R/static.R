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
#' With several gradients, for example an environmental and a spatial one, each
#' species has its own response to each gradient and the expected abundance is
#' the product of the responses, as in Peres-Neto et al. (2012).
#'
#' @param traits Data frame from [species_traits()], or any data frame with
#'   columns `optimum`, `sigma` and `r_max`, or a list of them, one per
#'   gradient, with the same species. Here `r_max` is the maximum expected
#'   abundance, not a growth rate, so set `rmax_range` accordingly. With several
#'   gradients the maxima multiply, so a range such as `c(1, 1)` in all but one
#'   of them keeps abundances on the scale of the first.
#' @param env Value of the gradient at each site, or a list or data frame with
#'   one gradient per element, in the order of `traits`.
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
#'
#' # an environmental and a spatial gradient
#' space <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
#'                         rmax_range = c(1, 1))
#' gradients <- data.frame(env = runif(100, 0, 100), space = 1:100)
#' simulate_static(list(traits, space), gradients)
#' @export
simulate_static <- function(traits, env, incidence = FALSE) {
  trait_list <- if (is.data.frame(traits)) list(traits) else traits
  env_list <- if (is.list(env)) as.list(env) else list(env)
  stopifnot(length(trait_list) == length(env_list),
            all(vapply(env_list, is.numeric, logical(1))),
            length(unique(lengths(env_list))) == 1, length(env_list[[1]]) >= 1,
            all(vapply(trait_list, function(x) all(x$r_max >= 0), logical(1))),
            all(vapply(trait_list, function(x) identical(rownames(x), rownames(trait_list[[1]])),
                       logical(1))))
  expected <- Reduce(`*`, Map(growth_rates, unname(env_list), trait_list))
  abundance <- matrix(stats::rpois(length(expected), expected), nrow(expected),
                      ncol(expected), dimnames = list(
                        rownames(expected) %||% paste0("site", seq_along(env_list[[1]])),
                        rownames(trait_list[[1]])))
  if (incidence) abundance[] <- (abundance > 0) * 1L
  structure(list(abundance = abundance, history = NULL, env = env, traits = traits),
            class = "mc_sim")
}
