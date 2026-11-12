#' Evolution of local populations
#'
#' Settings for the evolution of the niche optimum within each population (a species at a site)
#' during [simulate_metacommunity()] or [mc_dynamics()]. Without them, every population of a
#' species keeps the optimum of the species. With them, each population has its own genetic
#' mean \eqn{\bar z}{zbar}, which responds to selection by the local environment, drifts, mixes
#' with immigrants and starts from a sample of its founders.
#'
#' At each step, in the order of the model:
#'
#' 1. **Seeds** carry the optimum of the species in `traits`, the mean of a regional pool that
#'    does not evolve. They mix with residents in proportion to their numbers; a population
#'    founded by seeds alone starts as below.
#' 2. **Phenotype and growth.** The phenotype is \eqn{p = \bar z + b (E - E_0)}{p = zbar + b
#'    (E - E0)}, a linear reaction norm with slope `plasticity` (b) around the reference
#'    environment \eqn{E_0}{E0}, the mean environment of the sites at the first step. Growth
#'    rates follow [growth_rates()] with \eqn{p} in place of the optimum: \eqn{r = r_{max}
#'    \exp(-(E - p)^2 / (2\sigma^2))}{r = r_max exp(-(E - p)^2 / (2 sigma^2))}.
#' 3. **Selection and drift**, in each population present after reproduction (\eqn{N > 0}):
#'    \deqn{\Delta \bar z = \frac{1}{g} \frac{G}{\sigma^2 + P} (E - p) + \epsilon, \quad
#'    \epsilon \sim N\left(0, \frac{G}{g N_e}\right), \quad N_e = \nu N,}{Delta zbar = (1/g)
#'    G / (sigma^2 + P) (E - p) + e, e ~ N(0, G / (g Ne)), Ne = nu N,}
#'    the response to Gaussian stabilising selection of Lande (1976), with additive genetic
#'    variance \eqn{G}{G} (`genetic_var`), phenotypic variance \eqn{P = G + V_E}{P = G + VE}
#'    (`env_var`), the niche breadth \eqn{\sigma}{sigma} of the species as the width of the
#'    fitness function, \eqn{g}{g} steps per generation (`steps_per_generation`) and the ratio
#'    \eqn{\nu}{nu} of effective to census size (`ne_ratio`). Growth uses the mean phenotype and
#'    ignores the load of the variance around it, so that with \eqn{G = 0}{G = 0} the model is
#'    the one without evolution.
#' 4. **Gene flow.** After dispersal, the genetic mean of a population is the mean of residents
#'    and immigrants weighted by their numbers. Immigrants carry the mean of the populations
#'    they come from, weighted by the emigrants each sends through the dispersal kernel (or, on
#'    spatially implicit landscapes, by the emigrants of all other sites). A population founded
#'    at this step (absent at the end of the last one) starts from its founders: with `founder
#'    = TRUE`, the mean of the \eqn{k}{k} founders is drawn from \eqn{N(\bar z_{source},
#'    G/k)}{N(zbar_source, G / k)}, the sampling of a few individuals from the source.
#'    Migrants keep the phenotype of their site of origin until they reproduce, so the
#'    phenotype recorded after dispersal mixes plastic responses to different sites; their
#'    offspring develop in the new site.
#' 5. Populations lost to extirpation, competition or the `pool` lose their mean.
#'
#' Selection and drift draw no random numbers when \eqn{G = 0}{G = 0} or `ne_ratio = Inf`, so
#' with `genetic_var = 0` and `plasticity = 0` the abundances are the same as without
#' evolution for the same seed.
#'
#' Two cases have known answers, checked in the tests of the package. A population alone,
#' without selection (very wide niche) and of constant size \eqn{N}{N}, accumulates a variance
#' of \eqn{G t / (\nu N)}{G t / (nu N)} after \eqn{t}{t} generations (Lande 1976). Two
#' populations of the same size in environments \eqn{E_1}{E1} and \eqn{E_2}{E2}, exchanging a
#' share \eqn{m}{m} of individuals each generation, with selection \eqn{a = G / (\sigma^2 +
#' P)}{a = G / (sigma^2 + P)} and no drift, reach the divergence \eqn{d^* = (1 - 2m) a (E_1 -
#' E_2) / (1 - (1 - 2m)(1 - a))}{d* = (1 - 2m) a (E1 - E2) / (1 - (1 - 2m)(1 - a))}, which
#' follows from the recursion of the steps above.
#'
#' @param genetic_var Additive genetic variance of the optimum \eqn{G}{G}, in squared units of
#'   the environment, one value or one per species.
#' @param ne_ratio Ratio of effective to census population size \eqn{\nu}{nu}; `Inf` turns drift
#'   off.
#' @param plasticity Slope of the reaction norm b, one value or one per species: 0 for a fixed
#'   phenotype, 1 for a phenotype that tracks the environment fully around \eqn{E_0}{E0}.
#' @param env_var Environmental (non-heritable) variance of the phenotype \eqn{V_E}{VE}, one
#'   value or one per species.
#' @param steps_per_generation Time steps per generation \eqn{g}{g}; selection and drift per
#'   step are divided by it.
#' @param founder Should populations founded at a step start from the sampled mean of their
#'   founders (`TRUE`) or from the mean of their source (`FALSE`)?
#' @param env,traits The environment of each site (a vector, or a function of the step that
#'   returns one) and the traits of the species (from [species_traits()]). Filled in by
#'   [simulate_metacommunity()]; needed only with [mc_dynamics()].
#' @return An `mc_evolution` object, to pass as `evolution`.
#' @references Lande, R. (1976). Natural selection and random genetic drift in phenotypic
#'   evolution. *Evolution*, 30, 314--334. \doi{10.1111/j.1558-5646.1976.tb00911.x}
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' ls <- landscape_implicit(20)
#' tr <- treeshapesim::sim_tree_shape(10)
#' traits <- species_traits(tr)
#' sim <- simulate_metacommunity(ls, traits, competition_matrix(tr, "phylo"),
#'                               dispersal = 0.001, n_steps = 300,
#'                               evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5))
#' # local genetic means follow the environment of each site
#' head(sim$evolution$genetic[, 1:4])
#' @export
evolution_control <- function(genetic_var, ne_ratio = 1, plasticity = 0, env_var = 0,
                              steps_per_generation = 1, founder = TRUE, env = NULL,
                              traits = NULL) {
  if (missing(genetic_var)) stop("`genetic_var` is needed", call. = FALSE)
  stopifnot(is.numeric(genetic_var), all(genetic_var >= 0),
            is.numeric(ne_ratio), length(ne_ratio) == 1, ne_ratio > 0,
            is.numeric(plasticity), is.numeric(env_var), all(env_var >= 0),
            is.numeric(steps_per_generation), length(steps_per_generation) == 1,
            steps_per_generation > 0, is.logical(founder), length(founder) == 1)
  structure(list(genetic_var = genetic_var, ne_ratio = ne_ratio, plasticity = plasticity,
                 env_var = env_var, steps_per_generation = steps_per_generation,
                 founder = founder, env = env, traits = traits), class = "mc_evolution")
}

#' @export
print.mc_evolution <- function(x, ...) {
  cat("<mc_evolution> G =", format(x$genetic_var), "| Ne/N =", format(x$ne_ratio),
      "| plasticity =", format(x$plasticity), "| VE =", format(x$env_var),
      "| steps per generation =", format(x$steps_per_generation),
      "| founder sampling:", x$founder, "\n")
  invisible(x)
}

# the state of evolution inside mc_dynamics(): checks the settings against the sites and
# species and returns the functions that the loop calls
evolution_state <- function(evolution, m, k_all, species) {
  if (!inherits(evolution, "mc_evolution")) {
    stop("`evolution` must come from `evolution_control()`", call. = FALSE)
  }
  if (is.null(evolution$env) || is.null(evolution$traits)) {
    stop("`evolution` needs `env` and `traits` with `mc_dynamics()`", call. = FALSE)
  }
  traits <- as_traits(evolution$traits)
  if (nrow(traits) != k_all) {
    stop("`traits` of `evolution` need one row per species (", k_all, ")", call. = FALSE)
  }
  if (!is.null(species) && all(species %in% rownames(traits))) traits <- traits[species, ]
  env_fun <- if (is.function(evolution$env)) evolution$env else {
    function(step) evolution$env
  }
  e1 <- env_fun(1)
  if (length(e1) != m) stop("`env` of `evolution` needs one value per site (", m, ")",
                            call. = FALSE)
  each <- function(x) {
    if (!length(x) %in% c(1, k_all)) stop("settings of `evolution` take one value or one per ",
                                          "species (", k_all, ")", call. = FALSE)
    rep_len(x, k_all)
  }
  list(env = function(step) unname(as.numeric(env_fun(step))), ref = mean(e1),
       optimum = traits$optimum, sigma = traits$sigma, r_max = traits$r_max,
       g = each(evolution$genetic_var), b = each(evolution$plasticity),
       ve = each(evolution$env_var), nu = evolution$ne_ratio,
       per = evolution$steps_per_generation, founder = evolution$founder)
}

# growth rates of the phenotypes p (sites x species) in environment e, as in growth_rates()
evolution_growth <- function(p, e, sigma, r_max) {
  r <- exp(-sweep((e - p)^2, 2, 2 * sigma^2, "/"))
  sweep(r, 2, r_max, "*")
}

# mean trait of the immigrants into each site, from the emigrants em and their trait x
# (sites x species), through the dispersal matrix or, without it, from all other sites
immigrant_mean <- function(em, x, disp_matrix, fallback) {
  x[is.na(x)] <- 0
  if (!is.null(disp_matrix)) {
    num <- crossprod(disp_matrix, em * x)
    den <- crossprod(disp_matrix, em)
  } else {
    num <- sweep(-(em * x), 2, colSums(em * x), "+")
    den <- sweep(-em, 2, colSums(em), "+")
  }
  out <- num / den
  out[!is.finite(out)] <- fallback[!is.finite(out)]
  out
}
