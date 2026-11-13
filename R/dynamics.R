#' Simulate a metacommunity
#'
#' Runs the process-based model of Thompson et al. (2020) on a landscape. At
#' each time step, abundances are drawn as
#' \eqn{N' \sim \mathrm{Poisson}(r N / (1 + N \alpha'))}{N' ~ Poisson(r N / (1 + N alpha'))},
#' where \eqn{r} holds the growth rates of each species at each site ([growth_rates()]) and
#' \eqn{\alpha}{alpha} the competition coefficients ([competition_matrix()]). Each individual
#' then emigrates with
#' probability `dispersal`, and emigrants are redistributed among sites:
#'
#' * spatially implicit landscapes ([landscape_implicit()]): the emigrants of
#'   each species are distributed multinomially, each site receiving them with
#'   probability proportional to the number of emigrants from the other sites
#'   (a uniform kernel; an individual may return to its own site);
#' * spatially explicit landscapes ([landscape_grid()]): the number of
#'   immigrants of each species at each site is drawn from a Poisson
#'   distribution whose mean is the expected flow of emigrants through the
#'   dispersal kernel. This approximation does not conserve the number of
#'   dispersing individuals exactly, but has the right expectation. The code of
#'   Thompson et al. (2020) does the same. Few emigrants may then produce no
#'   immigrant, so with `extirpation` and low dispersal, recolonisation is rarer
#'   than in mcomsimr, which conserves dispersers. With `immigration = "multinomial"`,
#'   the emigrants of each species at each site are instead distributed among the other
#'   sites multinomially, with the probabilities of the kernel, which conserves them. It
#'   costs about as much as the Poisson draw at low dispersal and up to 20 times more at high
#'   dispersal, where the two give nearly the same results.
#'
#' After dispersal, each population is lost with probability `extirpation`,
#' independently of its size. Thompson et al. (2020) used 0.002 in their
#' competition-colonisation scenario and no extirpation in the others. With `disturbance`,
#' each site then loses all its populations at once with that probability, as the patch
#' extinctions of Leibold et al. (2019); the site is recolonised by dispersal.
#'
#' During colonisation, Poisson(`seed_mean`) individuals of every species are
#' added to every site at each of the `seed_steps`. After the last of these,
#' species absent from all sites can never return and are no longer computed.
#'
#' Distributions from elsewhere, such as ranges, biogeographic reconstructions or the output of
#' gen3sis, enter as `pool` and `initial`. With `pool`, a species can only occur where its value
#' is `TRUE` (or 1): seeds go only there, and individuals that disperse to other sites are lost,
#' as in a hostile matrix. A species without any site in `pool` dies out. `initial` replaces the
#' first draw of Poisson(`seed_mean`) individuals with given abundances, or with one individual
#' where a matrix of presences has 1; seeding still follows `seed_steps`, so use
#' `seed_steps = numeric(0)` to start from `initial` alone. Both take one row per site and one
#' column per species, named by species. Rows are matched by name when they have names (`site1`,
#' `site2` and so on in [simulate_metacommunity()]) and by position otherwise. Without them,
#' results with the same seed are the same as before they existed. A tree from elsewhere, such as
#' one simulated with gen3sis or a model of state-dependent diversification, enters through
#' [species_traits()]; drop extinct lineages first, since they have no distribution.
#'
#' @section Colonisation history:
#' The simulation records, for every population (a species at a site) present at the last step,
#' when and how it was founded, for arrival-dependent trait deviations in
#' [sample_individuals()] (founder effects and monopolization):
#'
#' * `arrival`, the step at which the population was last founded: the first step at whose end
#'   it is present after being absent at the end of the previous one. Populations present at
#'   the end of step 1, from `initial` or the first draw of seeds, have arrival 1; a population
#'   lost to extirpation or competition and founded again gets the step of the new arrival.
#' * `propagule`, the number of founders: the seeds (or `initial` individuals) and immigrants
#'   that reached the site at the step of arrival.
#' * `immigrants`, the number of immigrants and seeds that reached the population after its
#'   arrival.
#' * `migrant_load`, the sum over the steps since arrival of the share of the population made
#'   up of newcomers: immigrants over the abundance after dispersal, plus seeds over the
#'   abundance before reproduction. Each newcomer carries the species mean, so under complete
#'   mixing a deviation of the population shrinks by a factor \eqn{1 - m_t \approx
#'   e^{-m_t}}{1 - m_t ~ exp(-m_t)} at a step where newcomers make up a share \eqn{m_t}{m_t},
#'   and by about \eqn{e^{-L}}{exp(-L)} over the steps with migrant load \eqn{L}.
#'
#' Each is a sites x species matrix, `NA` for populations absent at the last step. Recording
#' them takes no random draws, so abundances are the same as without them for the same seed.
#' [colonisation_history()] returns them as a table.
#'
#' @section Article and archived code:
#' The article of Thompson et al. (2020) and the code archived with it
#' (\doi{10.5281/zenodo.3708350}) differ in points that change results. Beta is
#' gamma / alpha in the code and in their figures, but gamma - alpha in their
#' Table 1. Their environment is converted to ranks, so it is uniform between 0
#' and 1, and the lowest dispersal is 0.00001, not 0.0001. Emigrants are Poisson in
#' the code and binomial in the article, which comdynamite follows. In their
#' competition-colonisation scenario, the code also halves the immigrants. See
#' the article "Validating the model of Thompson et al. (2020)" on the package
#' website for the complete list.
#'
#' @param landscape An `mc_landscape` object.
#' @param traits Species traits from [species_traits()], or, for a niche with several
#'   environmental axes, a list of them, one per column of `env`, with the same species; see
#'   [growth_rates()].
#' @param env Optional environmental axes of each site, a matrix or data frame with one row per
#'   site and one column per element of `traits` (such as the `env` of several landscapes on
#'   the same grid). Needed when `traits` is a list; replaces the `env` of `landscape`. With
#'   `env_time`, it gives the axes (names and number), and the environment of each step comes
#'   from `env_time`. Columns without names take the names of `traits`, or `axis1`, `axis2`,
#'   ...
#' @param alpha Competition matrix from [competition_matrix()]. Species are matched to the
#'   rows of `traits` by name, or by position if `alpha` has no names.
#' @param dispersal Emigration probability, one value or one per species.
#' @param extirpation Probability that a population is lost at each time step,
#'   one value or one per species.
#' @param disturbance Probability that a site loses all its populations at each time step
#'   (one value; 0, the default, for none).
#' @param n_steps Number of time steps.
#' @param env_time Optional matrix (sites x `n_steps`) with the environment at
#'   each step, e.g. from [env_over_time()]. If `NULL`, the environment of the
#'   landscape is constant. With `env` (several axes), a list with one such matrix per column
#'   of `env` (matched by name when named), or an array of sites x steps x axes.
#' @param seed_steps Steps at which individuals are added during colonisation.
#' @param seed_mean Mean number of individuals of each species added to each
#'   site at each seeding step (and at the start).
#' @param record_steps Optional steps at which abundances are recorded.
#' @param immigration How emigrants reach other sites on spatially explicit landscapes:
#'   `"poisson"` (the default, as in the code of Thompson et al. 2020) or `"multinomial"`,
#'   which conserves them. Spatially implicit landscapes are always multinomial.
#' @param pool Optional matrix (sites x species), logical or 0/1, with the sites where each
#'   species can occur.
#' @param pool_from Optional matrix (sites x species) with the first step at which each species
#'   can occur in each site (1 for the start), for example to open a new site to two species at
#'   different times, as in the colonisation-time experiments of Urban & De Meester (2009).
#'   Combined with `pool`; until then, immigrants and seeds of that species do not establish
#'   there.
#' @param initial Optional matrix (sites x species) of abundances or 0/1 presences at the start,
#'   only where `pool` allows.
#' @param evolution Optional settings from [evolution_control()] for the evolution of the niche
#'   optimum within each population: selection by the local environment, drift, gene flow,
#'   founder sampling and plasticity. `NULL` (the default) keeps the optimum of each species in
#'   all its populations, as before. With a niche of several axes, one trait evolves per axis,
#'   with covariance matrices \eqn{G}{G} and \eqn{E}{E} (see the section Several traits of
#'   [evolution_control()] and [genetic_matrices()]).
#' @return An `mc_sim` object: a list with `abundance` (sites x species matrix
#'   at the last step), `history` (sites x species x recorded steps array, or
#'   `NULL`), and the `landscape`, `traits`, `alpha`, `dispersal`, `extirpation` and `pool`
#'   used, `env` when given, and `tree` when the traits came from [species_traits()] with a
#'   tree, the colonisation history (`arrival`, `propagule`, `immigrants`, `migrant_load`;
#'   see the section Colonisation history) and `n_steps`. With `evolution`, also `evolution`, a
#'   list with the genetic mean of the optimum in each population (`genetic`), the plastic part
#'   of its phenotype at the last step (`plastic`), the phenotype after dispersal (`phenotype`,
#'   which mixes residents and immigrants), the genetic means and the environment at
#'   `record_steps` (`history`, `environment_history`), the reference environment and the
#'   settings, and, with `mutation_var`, the genetic variance of each population at the last
#'   step (`genetic_variance`); `NA` where a population is absent. With several evolving traits, `genetic`,
#'   `plastic` and `phenotype` are arrays of sites x species x traits, `history` has the
#'   recorded steps as a fourth dimension, and `genetic_var`, `env_var` (one matrix per species)
#'   and `plasticity` (species x traits) are kept. With a niche of several axes, [as_mc_data()] keeps every axis (`env` and one
#'   `optimum_` column of traits per axis), while [summary_stats()] uses the first axis unless
#'   given `trait` and `env`, and [sample_individuals()] needs a single axis unless the traits
#'   evolved (`source = "evolution"`, one column per trait).
#' @references Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z.,
#'   Ptacnik, R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A
#'   process-based metacommunity framework linking local and regional scale
#'   community ecology. *Ecology Letters*, 23, 1314--1329.
#'   \doi{10.1111/ele.13568}
#'
#'   Urban, M. C., & De Meester, L. (2009). Community monopolization: local adaptation enhances
#'   priority effects in an evolving metacommunity. *Proceedings of the Royal Society B*, 276,
#'   4129--4138. \doi{10.1098/rspb.2009.1382}
#'
#'   Leibold, M. A., Urban, M. C., De Meester, L., Klausmeier, C. A., & Vanoverbeke, J.
#'   (2019). Regional neutrality evolves through local adaptive niche evolution.
#'   *Proceedings of the National Academy of Sciences*, 116, 2612--2617.
#'   \doi{10.1073/pnas.1808615116}
#' @family simulate metacommunities
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' ls <- landscape_implicit(30)
#' tr <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(ls, species_traits(tr),
#'                               competition_matrix(tr, "phylo"),
#'                               dispersal = 0.01, n_steps = 200)
#' sim
#' # a niche with two environmental axes
#' ls2 <- landscape_implicit(30)
#' niche2 <- species_traits(tr, rmax_range = c(1, 1))
#' sim2 <- simulate_metacommunity(ls, list(species_traits(tr), niche2),
#'                                competition_matrix(tr, "phylo"), n_steps = 200,
#'                                env = cbind(ls$env, ls2$env))
#' @export
simulate_metacommunity <- function(landscape, traits, alpha, dispersal = 0.01,
                                   n_steps = 1000, env_time = NULL,
                                   seed_steps = seq(10, 100, 10), seed_mean = 0.5,
                                   record_steps = NULL, extirpation = 0, pool = NULL,
                                   initial = NULL, immigration = c("poisson", "multinomial"),
                                   env = NULL, evolution = NULL, disturbance = 0,
                                   pool_from = NULL) {
  stopifnot(inherits(landscape, "mc_landscape"))
  immigration <- match.arg(immigration)
  several <- is.list(traits) && !is.data.frame(traits)
  if (several || !is.null(env)) {
    trait_list <- if (several) lapply(traits, as_traits) else list(as_traits(traits))
    if (is.null(env)) stop("a list of traits needs `env`, one column per element", call. = FALSE)
    unnamed <- is.null(colnames(env)) && !is.data.frame(env)
    env <- as.data.frame(env)
    if (unnamed) {
      names(env) <- (if (several) names(traits)) %||% paste0("axis", seq_len(ncol(env)))
    }
    if (nrow(env) != length(landscape$env) || ncol(env) != length(trait_list)) {
      stop("`env` needs one row per site (", length(landscape$env), ") and one column per ",
           "set of traits (", length(trait_list), ")", call. = FALSE)
    }
    traits <- trait_list[[1]]
    if (!is.null(env_time)) env_time <- env_time_axes(env_time, env, n_steps)
  }
  if (!identical(dim(alpha), rep(nrow(traits), 2L))) {
    stop("`alpha` needs one row and one column per species (", nrow(traits), ")", call. = FALSE)
  }
  species <- rownames(traits)
  if (!is.null(rownames(alpha)) && !is.null(colnames(alpha))) {
    if (!setequal(rownames(alpha), species) || !setequal(colnames(alpha), species)) {
      stop("names of species in `alpha` and `traits` differ", call. = FALSE)
    }
    alpha <- alpha[species, species]
  }
  sites <- paste0("site", seq_along(landscape$env))
  r <- if (!is.null(env) && !is.null(env_time)) {
    rownames(env) <- sites
    function(step) {
      now <- as.data.frame(lapply(env_time, function(x) x[, step]), row.names = sites)
      growth_rates(now, if (several) trait_list else traits)
    }
  } else if (!is.null(env)) {
    rownames(env) <- sites
    growth_rates(env, if (several) trait_list else traits)
  } else if (is.null(env_time)) {
    growth_rates(stats::setNames(landscape$env, sites), traits)
  } else {
    stopifnot(nrow(env_time) == length(landscape$env), ncol(env_time) >= n_steps)
    function(step) growth_rates(stats::setNames(env_time[, step], sites), traits)
  }
  if (!is.null(evolution)) {
    if (!inherits(evolution, "mc_evolution")) {
      stop("`evolution` must come from `evolution_control()`", call. = FALSE)
    }
    mv <- is.matrix(evolution$genetic_var) || is.list(evolution$genetic_var) ||
      is.matrix(evolution$env_var) || is.list(evolution$env_var)
    if (several || mv) {
      if (!several) {
        stop("`evolution` with covariance matrices needs a niche of several axes (a list of ",
             "`traits` and `env`), one trait per axis", call. = FALSE)
      }
      evolution$traits <- stats::setNames(trait_list, names(env))
      evolution$env <- if (is.null(env_time)) as.matrix(env) else {
        function(step) vapply(env_time, function(x) x[, step], numeric(nrow(env)))
      }
    } else {
    evolution$traits <- traits
    evolution$env <- if (!is.null(env) && !is.null(env_time)) {
      function(step) env_time[[1]][, step]
    } else if (!is.null(env)) env[[1]] else if (is.null(env_time)) landscape$env else
      function(step) env_time[, step]
    }
  }
  out <- mc_dynamics(r, alpha, dispersal, n_steps = n_steps, seed_steps = seed_steps,
                     seed_mean = seed_mean, disp_matrix = landscape$dispersal,
                     record_steps = record_steps, sites = sites, species = rownames(traits),
                     extirpation = extirpation, pool = pool, initial = initial,
                     immigration = immigration, evolution = evolution,
                     disturbance = disturbance, pool_from = pool_from)
  tree <- attr(traits, "tree")
  if (several) traits <- trait_list
  structure(c(out, list(landscape = landscape, env = env, traits = traits, alpha = alpha,
                        dispersal = dispersal, extirpation = extirpation,
                        disturbance = disturbance, pool = out$pool, tree = tree)),
            class = "mc_sim")
}

# the environment through time on several axes: a list of matrices (sites x steps), one per
# column of `env` (matched by name when named), or an array of sites x steps x axes
env_time_axes <- function(env_time, env, n_steps) {
  if (is.array(env_time) && length(dim(env_time)) == 3) {
    env_time <- lapply(seq_len(dim(env_time)[3]), function(a) env_time[, , a])
  } else if (is.matrix(env_time) && ncol(env) == 1) {
    env_time <- list(env_time)
  }
  if (!is.list(env_time) || length(env_time) != ncol(env)) {
    stop("with `env`, `env_time` needs one matrix (sites x steps) per column of `env` (",
         ncol(env), "), as a list or an array of sites x steps x axes", call. = FALSE)
  }
  if (!is.null(names(env_time)) && all(names(env) %in% names(env_time))) {
    env_time <- env_time[names(env)]
  }
  ok <- vapply(env_time, function(x) {
    is.numeric(x) && is.matrix(x) && nrow(x) == nrow(env) && ncol(x) >= n_steps && !anyNA(x)
  }, logical(1))
  if (!all(ok)) {
    stop("each matrix of `env_time` needs one row per site (", nrow(env), ") and at least ",
         "`n_steps` (", n_steps, ") columns, without NA", call. = FALSE)
  }
  stats::setNames(lapply(env_time, unname), names(env))
}

#' Metacommunity dynamics from growth rates
#'
#' The engine of [simulate_metacommunity()], for users who want to supply
#' their own growth rates (e.g. responses to several environmental variables).
#'
#' @inheritParams simulate_metacommunity
#' @param r Sites x species matrix of density-independent growth rates, or a
#'   function of the time step that returns such a matrix.
#' @param disp_matrix Row-stochastic dispersal matrix (sites x sites), or
#'   `NULL` for uniform dispersal.
#' @param sites,species Optional names of sites and species.
#' @return List with `abundance`, `history`, the colonisation history (`arrival`,
#'   `propagule`, `immigrants`, `migrant_load`) and `n_steps`, as in
#'   [simulate_metacommunity()], and `pool` as a logical matrix, or `NULL` without it.
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' # one species alone converges to (r - 1) / alpha = (5 - 1) / 0.05 = 80
#' out <- mc_dynamics(matrix(5, 20, 1), matrix(0.05), dispersal = 0, n_steps = 200)
#' mean(out$abundance)
#' @export
mc_dynamics <- function(r, alpha, dispersal, n_steps = 1000,
                        seed_steps = seq(10, 100, 10), seed_mean = 0.5,
                        disp_matrix = NULL, record_steps = NULL,
                        sites = NULL, species = NULL, extirpation = 0, pool = NULL,
                        initial = NULL, immigration = c("poisson", "multinomial"),
                        evolution = NULL, disturbance = 0, pool_from = NULL) {
  # the state lives in this frame; the preparation, each step and the end are evaluated in it, so
  # that simulate_evolving_metacommunity() can run the same steps and change species in between
  mc_setup(environment())
  for (step in seq_len(n_steps)) mc_step(environment(), step)
  mc_finish(environment())
}

# preparation of mc_dynamics(), evaluated in the frame that holds its arguments
mc_setup <- function(st) {
  eval(quote({
    immigration <- match.arg(immigration)
    r_fun <- if (is.function(r)) r else function(step) r
    r1 <- r_fun(1)
    m <- nrow(r1)
    k_all <- ncol(r1)
    stopifnot(all(dim(alpha) == k_all), n_steps >= 1)
    if (!is.null(disp_matrix)) stopifnot(all(dim(disp_matrix) == m))
    dispersal <- rep_len(dispersal, k_all)
    extirpation <- rep_len(extirpation, k_all)
    stopifnot(all(dispersal >= 0 & dispersal <= 1), all(extirpation >= 0 & extirpation <= 1))
    if (!is.numeric(disturbance) || length(disturbance) != 1 || is.na(disturbance) ||
        disturbance < 0 || disturbance > 1) {
      stop("`disturbance` must be one probability between 0 and 1", call. = FALSE)
    }
    if (m == 1 && any(dispersal > 0)) {
      stop("dispersal needs at least two sites: use `dispersal = 0` with one site", call. = FALSE)
    }
    sites <- sites %||% rownames(r1) %||% paste0("site", seq_len(m))
    species <- species %||% colnames(r1) %||% paste0("sp", seq_len(k_all))
    if (!is.null(pool)) {
      pool <- site_species(pool, sites, species, "`pool`")
      if (!all(pool %in% c(0, 1))) stop("`pool` must be logical or 0/1", call. = FALSE)
      pool <- pool == 1
    }
    if (!is.null(initial)) {
      initial <- site_species(initial, sites, species, "`initial`")
      if (any(initial < 0) || any(initial != round(initial))) {
        stop("`initial` must hold whole numbers of individuals, zero or more", call. = FALSE)
      }
      if (!is.null(pool) && any(initial[!pool] > 0)) {
        stop("`initial` has individuals outside `pool`", call. = FALSE)
      }
      storage.mode(initial) <- "integer"
    }
    record_steps <- sort(unique(record_steps[record_steps >= 1 & record_steps <= n_steps]))
    history <- if (length(record_steps))
      array(0L, c(m, k_all, length(record_steps)),
            dimnames = list(sites, species, step_label(record_steps))) else NULL

    alive <- seq_len(k_all)
    k <- k_all
    p_em <- rep(dispersal, each = m)
    p_ex <- rep(extirpation, each = m)
    at <- t(alpha)
    if (!is.null(pool_from)) {
      pool_from <- site_species(pool_from, sites, species, "`pool_from`")
      if (any(is.na(pool_from)) || any(pool_from < 1)) {
        stop("`pool_from` must hold steps of 1 or more", call. = FALSE)
      }
      if (!is.null(initial) && any(initial[pool_from > 1] > 0)) {
        stop("`initial` has individuals where `pool_from` opens after the start", call. = FALSE)
      }
    }
    # sites open to each species at a step: `pool` and, with `pool_from`, from that step on
    inside_at <- function(step, cols) {
      x <- if (is.null(pool)) matrix(TRUE, m, length(cols)) else pool[, cols, drop = FALSE]
      if (!is.null(pool_from)) x <- x & pool_from[, cols, drop = FALSE] <= step
      x
    }
    inside <- if (is.null(pool) && is.null(pool_from)) NULL else inside_at(1, seq_len(k_all))
    n <- if (is.null(initial)) matrix(stats::rpois(m * k, seed_mean), m, k) else unname(initial)
    if (!is.null(inside)) n[!inside] <- 0L
    last_seed <- max(c(0, seed_steps))
    competition_hook <- NULL
    ev <- if (!is.null(evolution)) {
      if (!inherits(evolution, "mc_evolution")) {
        stop("`evolution` must come from `evolution_control()`", call. = FALSE)
      }
      if (evolution_is_mv(evolution)) evolution_state_mv(evolution, m, k_all, species) else
        evolution_state(evolution, m, k_all, species)
    }
    if (isTRUE(ev$imm_sel) && (is.null(disp_matrix) || immigration != "poisson")) {
      stop("`immigrant_selection` needs a spatially explicit landscape and ",
           "`immigration = \"poisson\"`", call. = FALSE)
    }
    if (isTRUE(ev$mv)) {
      # several traits: genetic means and phenotypes as sites x species x traits
      z <- ph <- array(NA_real_, c(m, k_all, ev$q))
      z_history <- if (length(record_steps))
        array(NA_real_, c(m, k_all, ev$q, length(record_steps)),
              dimnames = list(sites, species, ev$axes, step_label(record_steps))) else NULL
      fixed_z <- ev$fixed
      pool_mean <- array(rep(ev$optimum, each = m), c(m, k_all, ev$q))
      start_pop <- n > 0
      if (any(start_pop)) z <- mv_found(ev, z, pool_mean, n, start_pop, seq_len(k_all))
    } else if (!is.null(ev)) {
      # genetic mean and phenotype (after dispersal) of every population; NA where absent
      z <- matrix(NA_real_, m, k_all)
      ph <- matrix(NA_real_, m, k_all)
      z_history <- if (length(record_steps))
        array(NA_real_, c(m, k_all, length(record_steps)),
              dimnames = list(sites, species, step_label(record_steps))) else NULL
      # mean of k founders drawn from sources with mean `from` and, when they come from several
      # sources, variance `between` among their means, for the populations in `cells`
      found <- function(from, k, cells, cols, between = 0, g_cells = NULL) {
        g <- if (is.null(g_cells)) matrix(ev$g[cols], m, length(cols), byrow = TRUE)[cells] else
          g_cells[cells]
        v <- g + (between + 0 * from)[cells]
        out <- from[cells]
        if (ev$founder && any(g > 0)) {
          draw <- g > 0
          out[draw] <- out[draw] + stats::rnorm(sum(draw), 0, sqrt(v[draw] / k[cells][draw]))
        }
        out
      }
      # species without genetic variance keep the optimum of the pool exactly
      fixed_z <- ev$g == 0
      pool_mean <- matrix(ev$optimum, m, k_all, byrow = TRUE)
      start_pop <- n > 0
      if (any(start_pop)) z[start_pop] <- found(pool_mean, n, start_pop, seq_len(k_all))
      # genetic variance of every population, with `mutation_var` (NA where absent)
      pool_g <- matrix(ev$g, m, k_all, byrow = TRUE)
      vg <- NULL
      if (isTRUE(ev$dyn)) {
        vg <- matrix(NA_real_, m, k_all)
        vg[start_pop] <- pool_g[start_pop]
      }
    }
    # colonisation history of every population (sites x all species): step of the last arrival,
    # number of founders, and immigrants and migrant load since then; no random draws, so seeds
    # give the same abundances as before
    arrival <- matrix(NA_integer_, m, k_all)
    propagule <- matrix(NA_integer_, m, k_all)
    immigrants <- matrix(0, m, k_all)
    load <- matrix(0, m, k_all)
    present <- matrix(FALSE, m, k_all)
  }), st)
  invisible(st)
}

# one time step of mc_dynamics(), evaluated in the state `st`
mc_step <- function(st, step) {
  assign("step", step, envir = st)
  eval(quote({
    added <- 0L
    if (!is.null(pool_from)) inside <- inside_at(step, alive)
    if (step %in% seed_steps) {
      seeds <- stats::rpois(m * k, seed_mean)
      added <- if (is.null(inside)) seeds else seeds * inside
      if (isTRUE(ev$mv)) {
        added_m <- added + matrix(0L, m, k)
        za <- z[, alive, , drop = FALSE]
        pm <- pool_mean[, alive, , drop = FALSE]
        mix <- n > 0 & added_m > 0
        w_mix <- ifelse(mix, n / (n + added_m), 1)
        for (a in seq_len(ev$q)) za[, , a] <- w_mix * za[, , a] + (1 - w_mix) * pm[, , a]
        new_seed <- n == 0 & added_m > 0
        if (any(new_seed)) za <- mv_found(ev, za, pm, added_m, new_seed, alive)
        z[, alive, ] <- za
      } else if (!is.null(ev)) {
        # seeds carry the mean of the regional pool and mix with residents; populations
        # founded by seeds alone start from a sample of them
        added_m <- added + matrix(0L, m, k)
        za <- z[, alive, drop = FALSE]
        pm <- pool_mean[, alive, drop = FALSE]
        mix <- n > 0 & added_m > 0
        za[mix] <- (n[mix] * za[mix] + added_m[mix] * pm[mix]) / (n[mix] + added_m[mix])
        new_seed <- n == 0 & added_m > 0
        if (any(new_seed)) za[new_seed] <- found(pm, added_m, new_seed, alive)
        z[, alive] <- za
        if (isTRUE(ev$dyn)) {
          va <- vg[, alive, drop = FALSE]
          pg <- pool_g[, alive, drop = FALSE]
          w <- n[mix] / (n[mix] + added_m[mix])
          va[mix] <- w * va[mix] + (1 - w) * pg[mix]
          va[new_seed] <- pg[new_seed]
          vg[, alive] <- va
        }
      }
      n <- n + added
    }
    keep <- colSums(n) > 0
    if (step > last_seed && !all(keep)) {
      alive <- alive[keep]
      n <- n[, keep, drop = FALSE]
      at <- at[keep, keep, drop = FALSE]
      k <- length(alive)
      p_em <- rep(dispersal[alive], each = m)
      p_ex <- rep(extirpation[alive], each = m)
      if (!is.null(inside)) inside <- inside_at(step, alive)
    }
    if (k > 0) {
      if (is.null(ev)) {
        r_now <- r_fun(step)[, alive, drop = FALSE]
      } else if (isTRUE(ev$mv)) {
        e_now <- ev$env(step)
        za <- z[, alive, , drop = FALSE]
        pm <- pool_mean[, alive, , drop = FALSE]
        za[is.na(za)] <- pm[is.na(za)]
        if (any(fixed_z[alive])) za[, fixed_z[alive], ] <- pm[, fixed_z[alive], , drop = FALSE]
        pa <- mv_phenotype(ev, za, e_now, alive)
        r_now <- mv_growth(ev, pa, e_now, alive)
      } else {
        e_now <- ev$env(step)
        za <- z[, alive, drop = FALSE]
        za[is.na(za)] <- pool_mean[, alive, drop = FALSE][is.na(za)]
        if (any(fixed_z[alive])) za[, fixed_z[alive]] <- pool_mean[, alive, drop = FALSE][, fixed_z[alive]]
        pa <- za + sweep(matrix(e_now - ev$ref, m, k), 2, ev$b[alive], "*")
        if (isTRUE(ev$dyn)) {
          ga <- vg[, alive, drop = FALSE]
          ga[is.na(ga)] <- pool_g[, alive, drop = FALSE][is.na(ga)]
          r_now <- evolution_growth_cells(pa, e_now, ev$sigma[alive], ev$r_max[alive],
                                          sweep(ga, 2, ev$ve[alive], "+"))
        } else {
          r_now <- evolution_growth(pa, e_now, ev$sigma[alive], ev$r_max[alive],
                                    ev$g[alive] + ev$ve[alive])
        }
      }
      start <- n
      # crowding; simulate_evolving_metacommunity() may replace it with competition that
      # depends on the phenotypes in each site
      crowding <- if (is.null(competition_hook)) n %*% at else competition_hook(n, pa)
      lam <- r_now * n / (1 + crowding)
      n <- matrix(stats::rpois(m * k, lam), m, k)
      if (isTRUE(ev$mv)) {
        za <- mv_select(ev, za, pa, e_now, n, alive)
        pa <- mv_phenotype(ev, za, e_now, alive)
      } else if (!is.null(ev)) {
        # selection toward the local environment and drift, in populations present after
        # reproduction (Lande 1976); no draws without drift
        g <- if (isTRUE(ev$dyn)) ga else matrix(ev$g[alive], m, k, byrow = TRUE)
        width <- g + matrix(ev$sigma[alive]^2 + ev$ve[alive], m, k, byrow = TRUE)
        live <- n > 0
        za[live] <- za[live] + (g[live] / width[live]) * (e_now - pa)[live] / ev$per
        drift <- live & g > 0
        if (is.finite(ev$nu) && any(drift)) {
          za[drift] <- za[drift] +
            stats::rnorm(sum(drift), 0, sqrt(g[drift] / (ev$per * ev$nu * n[drift])))
        }
        za[!live] <- NA_real_
        if (isTRUE(ev$dyn)) {
          # mutation-drift dynamics of the genetic variance (Lynch & Hill 1986): a loss of
          # 1 / (2 Ne) and a gain of Vm per generation, equilibrium 2 Ne Vm
          vm <- matrix(ev$vm[alive], m, k, byrow = TRUE)
          keep_v <- pmax(1 - 1 / (2 * ev$nu * n), 0)^(1 / ev$per)
          ga[live] <- ga[live] * keep_v[live] + vm[live] / ev$per
          ga[!live] <- NA_real_
          vg[, alive] <- ga
        }
        # phenotype of the offspring, which develop in their site
        pa <- za + sweep(matrix(e_now - ev$ref, m, k), 2, ev$b[alive], "*")
      }
      em <- matrix(stats::rbinom(m * k, n, p_em), m, k)
      tot <- colSums(em)
      imm <- 0L
      if (any(tot > 0)) {
        if (!is.null(disp_matrix) && immigration == "multinomial") {
          imm <- place_emigrants(em, disp_matrix)
        } else if (!is.null(disp_matrix) && isTRUE(ev$imm_sel)) {
          # immigrants filtered by their fitness in the site they reach
          g_src <- if (isTRUE(ev$dyn)) ifelse(is.na(ga), 0, ga) else
            matrix(ev$g[alive], m, k, byrow = TRUE)
          sel <- immigrants_selected(em, ifelse(is.na(za), 0, za), ifelse(is.na(pa), 0, pa),
                                     e_now, g_src + matrix(ev$sigma[alive]^2 + ev$ve[alive], m,
                                                           k, byrow = TRUE),
                                     disp_matrix, g_src)
          imm <- matrix(stats::rpois(m * k, sel$lambda), m, k)
        } else if (!is.null(disp_matrix)) {
          imm <- matrix(stats::rpois(m * k, crossprod(disp_matrix, em)), m, k)
        } else {
          w <- sweep(-em, 2, tot, "+")
          imm <- matrix(0L, m, k)
          for (j in which(tot > 0)) imm[, j] <- stats::rmultinom(1, tot[j], w[, j])
        }
        if (isTRUE(ev$mv)) {
          flowed <- mv_flow(ev, za, pa, n, em, imm, disp_matrix, alive)
          za <- flowed$za
          pa <- flowed$pa
        } else if (!is.null(ev)) {
          # gene flow: residents and immigrants mixed by their numbers; immigrants carry the
          # mean (and the phenotype) of their sources; new populations start from a sample of
          # their founders
          imm_m <- imm + matrix(0L, m, k)
          res <- n - em
          if (isTRUE(ev$imm_sel)) {
            zi <- sel$mean
            between <- sel$between
            p_imm <- sel$phenotype
            vi <- sel$variance
          } else {
            zi <- immigrant_mean(em, za, disp_matrix, za)
            between <- pmax(immigrant_mean(em, za^2, disp_matrix, za^2) - zi^2, 0)
            p_imm <- immigrant_mean(em, pa, disp_matrix, pa)
            if (isTRUE(ev$dyn)) {
              g0 <- ifelse(is.na(ga), 0, ga)
              vi <- immigrant_mean(em, g0, disp_matrix, g0)
            }
          }
          pa_now <- pa
          both <- res > 0 & imm_m > 0
          if (isTRUE(ev$dyn)) {
            # variance of the merged population: the mean of the variances of residents and
            # immigrants; the difference between their means is not added, as if selection
            # removed it
            w <- res[both] / (res[both] + imm_m[both])
            ga[both] <- w * ga[both] + (1 - w) * vi[both]
          }
          za[both] <- (res[both] * za[both] + imm_m[both] * zi[both]) / (res[both] + imm_m[both])
          pa_now[both] <- (res[both] * pa[both] + imm_m[both] * p_imm[both]) /
            (res[both] + imm_m[both])
          founded <- res == 0 & imm_m > 0
          if (any(founded)) {
            if (isTRUE(ev$dyn)) {
              # founders carry the variance of their sources, less 1 / (2k) for k founders
              za[founded] <- found(zi, imm_m, founded, alive, g_cells = vi + between)
              ga[founded] <- (vi * pmax(1 - 1 / (2 * imm_m), 0))[founded]
            } else {
              za[founded] <- found(zi, imm_m, founded, alive, between)
            }
            pa_now[founded] <- p_imm[founded]
          }
          if (isTRUE(ev$dyn)) vg[, alive] <- ga
          pa <- pa_now
        }
        n <- n - em + imm
      }
      # no draws without extirpation, so that seeds give the same results as before
      if (any(p_ex > 0)) n[stats::runif(m * k) < p_ex] <- 0L
      if (disturbance > 0) n[stats::runif(m) < disturbance, ] <- 0L
      if (!is.null(inside)) n[!inside] <- 0L
      if (isTRUE(ev$mv)) {
        if (any(fixed_z[alive])) {
          za[, fixed_z[alive], ] <- pool_mean[, alive, , drop = FALSE][, fixed_z[alive], ,
                                                                         drop = FALSE]
        }
        gone <- rep(n == 0, ev$q)
        za[gone] <- NA_real_
        pa[gone] <- NA_real_
        z[, alive, ] <- za
        ph[, alive, ] <- pa
      } else if (!is.null(ev)) {
        if (any(fixed_z[alive])) za[, fixed_z[alive]] <- pool_mean[, alive, drop = FALSE][, fixed_z[alive]]
        za[n == 0] <- NA_real_
        pa[n == 0] <- NA_real_
        if (isTRUE(ev$dyn)) {
          va <- vg[, alive, drop = FALSE]
          va[n == 0] <- NA_real_
          vg[, alive] <- va
        }
        z[, alive] <- za
        ph[, alive] <- pa
      }
      was <- present[, alive, drop = FALSE]
      now <- n > 0
      new_pop <- now & !was
      stay <- now & was
      if (any(new_pop)) {
        cols <- arrival[, alive, drop = FALSE]
        cols[new_pop] <- step
        arrival[, alive] <- cols
        # a population absent at the end of the last step starts from its seeds (or `initial`)
        # and the immigrants of this step
        cols <- propagule[, alive, drop = FALSE]
        cols[new_pop] <- as.integer((start + imm)[new_pop])
        propagule[, alive] <- cols
        cols <- immigrants[, alive, drop = FALSE]
        cols[new_pop] <- 0
        immigrants[, alive] <- cols
        cols <- load[, alive, drop = FALSE]
        cols[new_pop] <- 0
        load[, alive] <- cols
      }
      if (any(stay) && (any(imm > 0) || any(added > 0))) {
        # immigrants and seeds into populations already there: their number, and the share of
        # the population they make up, summed over steps (the migrant load)
        imm <- imm + matrix(0L, m, k)
        added <- added + matrix(0L, m, k)
        arriving <- (imm + added)[stay]
        share <- imm[stay] / n[stay] + added[stay] / start[stay]
        cols <- immigrants[, alive, drop = FALSE]
        cols[stay] <- cols[stay] + arriving
        immigrants[, alive] <- cols
        cols <- load[, alive, drop = FALSE]
        cols[stay] <- cols[stay] + share
        load[, alive] <- cols
      }
      lost <- was & !now
      if (any(lost)) {
        cols <- arrival[, alive, drop = FALSE]
        cols[lost] <- NA_integer_
        arrival[, alive] <- cols
      }
      present[, alive] <- now
    }
    if (!is.null(history) && step %in% record_steps) {
      history[, alive, match(step, record_steps)] <- n
      if (isTRUE(ev$mv)) {
        z_history[, alive, , match(step, record_steps)] <- z[, alive, ]
      } else if (!is.null(ev)) {
        z_history[, alive, match(step, record_steps)] <- z[, alive]
      }
    }
  }), st)
  invisible(st)
}

# the output of mc_dynamics() from the state `st`
mc_finish <- function(st) {
  eval(quote({
    abundance <- matrix(0L, m, k_all, dimnames = list(sites, species))
    abundance[, alive] <- n
    absent <- abundance == 0
    arrival[absent] <- NA_integer_
    propagule[absent] <- NA_integer_
    immigrants[absent] <- NA
    load[absent] <- NA
    labels <- list(sites, species)
    dimnames(arrival) <- dimnames(propagule) <- dimnames(immigrants) <- dimnames(load) <- labels
    out <- list(abundance = abundance, history = history, pool = pool, arrival = arrival,
                propagule = propagule, immigrants = immigrants, migrant_load = load,
                n_steps = n_steps)
    if (isTRUE(ev$mv)) {
      gone <- rep(absent, ev$q)
      z[gone] <- NA_real_
      ph[gone] <- NA_real_
      e_end <- ev$env(n_steps)
      plastic <- z
      for (a in seq_len(ev$q)) plastic[, , a] <- outer(e_end[, a] - ev$ref[a], ev$b[, a])
      plastic[gone] <- NA_real_
      dimnames(z) <- dimnames(ph) <- dimnames(plastic) <- list(sites, species, ev$axes)
      env_history <- if (length(record_steps)) {
        array(vapply(record_steps, function(s) as.vector(ev$env(s)), numeric(m * ev$q)),
              c(m, ev$q, length(record_steps)),
              dimnames = list(sites, ev$axes, step_label(record_steps)))
      }
      out$evolution <- list(genetic = z, plastic = plastic, phenotype = ph, history = z_history,
                            environment = matrix(e_end, m, dimnames = list(sites, ev$axes)),
                            environment_history = env_history,
                            reference = stats::setNames(ev$ref, ev$axes),
                            genetic_var = stats::setNames(ev$G, species),
                            env_var = stats::setNames(ev$E, species),
                            plasticity = matrix(ev$b, k_all, dimnames = list(species, ev$axes)),
                            settings = structure(evolution[setdiff(names(evolution),
                                                                   c("env", "traits"))],
                                                 class = "mc_evolution"))
    } else if (!is.null(ev)) {
      z[absent] <- NA_real_
      ph[absent] <- NA_real_
      plastic <- sweep(matrix(ev$env(n_steps) - ev$ref, m, k_all), 2, ev$b, "*")
      plastic[absent] <- NA_real_
      dimnames(z) <- dimnames(ph) <- dimnames(plastic) <- labels
      env_history <- if (length(record_steps)) {
        matrix(vapply(record_steps, ev$env, numeric(m)), m,
               dimnames = list(sites, step_label(record_steps)))
      }
      if (!is.null(vg)) {
        vg[absent] <- NA_real_
        dimnames(vg) <- labels
      }
      out$evolution <- list(genetic = z, plastic = plastic, phenotype = ph, history = z_history,
                            environment = ev$env(n_steps), environment_history = env_history,
                            reference = ev$ref, genetic_variance = vg,
                            settings = structure(evolution[setdiff(names(evolution),
                                                                   c("env", "traits"))],
                                                 class = "mc_evolution"))
    }
    out
  }), st)
}

# immigrants when the emigrants of each species at each site go to the other sites
# multinomially, with the probabilities of the rows of `disp_matrix`
place_emigrants <- function(em, disp_matrix) {
  m <- nrow(em)
  immigrants <- matrix(0L, m, ncol(em))
  for (cell in which(em > 0)) {
    from <- (cell - 1L) %% m + 1L
    species <- (cell - 1L) %/% m + 1L
    immigrants[, species] <- immigrants[, species] +
      drop(stats::rmultinom(1, em[cell], disp_matrix[from, ]))
  }
  immigrants
}

# a sites x species matrix in the order of `sites` and `species`: rows by name when they have
# names, by position otherwise; columns by species name
site_species <- function(x, sites, species, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !(is.numeric(x) || is.logical(x))) {
    stop(name, " must be a matrix of sites x species", call. = FALSE)
  }
  if (is.null(colnames(x))) stop(name, " needs species names in its columns", call. = FALSE)
  missing <- setdiff(species, colnames(x))
  if (length(missing)) {
    stop(name, " has no column for the species ", species_list(missing), call. = FALSE)
  }
  extra <- setdiff(colnames(x), species)
  if (length(extra)) {
    warning(name, " has species that are not simulated, dropped: ", species_list(extra),
            call. = FALSE)
  }
  rows <- if (is.null(rownames(x))) {
    if (nrow(x) != length(sites)) {
      stop(name, " has ", nrow(x), " rows and the landscape ", length(sites), " sites",
           call. = FALSE)
    }
    seq_along(sites)
  } else {
    match(sites, rownames(x))
  }
  if (anyNA(rows)) {
    stop(name, " has no row for the sites ", species_list(sites[is.na(rows)]), call. = FALSE)
  }
  if (anyNA(x)) stop(name, " has missing values", call. = FALSE)
  out <- x[rows, species, drop = FALSE]
  dimnames(out) <- list(sites, species)
  out
}

#' Colonisation history of a simulated metacommunity
#'
#' The arrival step, number of founders, immigrants since arrival and migrant load of every
#' population present at the last step of a simulation, recorded by
#' [simulate_metacommunity()] (see its section Colonisation history).
#'
#' @param x An `mc_sim` object from [simulate_metacommunity()].
#' @return Data frame with one row per population present at the last step and the columns
#'   `site`, `species`, `abundance`, `arrival`, `age` (`n_steps - arrival`), `propagule`,
#'   `immigrants` and `migrant_load`, ordered by site and species.
#' @family simulate metacommunities
#' @seealso [sample_individuals()], whose `founder_sd` and `monopolization` use this history.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(10)
#' sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree),
#'                               competition_matrix(tree, "phylo"), n_steps = 100,
#'                               extirpation = 0.01)
#' head(colonisation_history(sim))
#' @export
colonisation_history <- function(x) {
  check_history(x)
  a <- x$abundance
  present <- which(a > 0)
  m <- nrow(a)
  site_of <- (present - 1) %% m + 1
  species_of <- (present - 1) %/% m + 1
  by_site <- order(site_of, species_of)
  present <- present[by_site]
  data.frame(site = rownames(a)[site_of[by_site]], species = colnames(a)[species_of[by_site]],
             abundance = a[present], arrival = x$arrival[present],
             age = x$n_steps - x$arrival[present], propagule = x$propagule[present],
             immigrants = x$immigrants[present], migrant_load = x$migrant_load[present],
             stringsAsFactors = FALSE)
}

check_history <- function(x) {
  if (!inherits(x, "mc_sim")) {
    stop("`x` must be a simulation from simulate_metacommunity()", call. = FALSE)
  }
  if (is.null(x$arrival)) {
    stop("`x` has no colonisation history: it comes from simulate_metacommunity() only, ",
         "and simulations made before comdynamite 0.6.3.9000 must be run again", call. = FALSE)
  }
  invisible(x)
}

#' @export
print.mc_sim <- function(x, ...) {
  a <- x$abundance
  cat("<mc_sim>", nrow(a), "sites x", ncol(a), "species\n")
  cat("species persisting:", sum(colSums(a) > 0), "of", ncol(a), "\n")
  cat("mean local richness:", format(mean(rowSums(a > 0)), digits = 3), "\n")
  if (!is.null(x$history)) cat("recorded steps:", dim(x$history)[3], "\n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
