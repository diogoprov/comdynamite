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
#'   dispersal, where the two give nearly the same results. Use it for low dispersal (below
#'   about 0.01), especially with `extirpation` or an environment that changes in time
#'   (`env_time`), and for studies of long-term persistence: against mcomsimr with an
#'   environment from [env_space_time()], the Poisson draw gave 16% to 26% lower local
#'   richness at dispersal 0.001 on a grid, and the multinomial draw 6% or less, of either
#'   sign (see the validation article on the package website).
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
#' @param dispersal Emigration probability, between 0 and 1: one value or one per species.
#' @param extirpation Probability that a population is lost at each time step, between 0 and
#'   1: one value or one per species.
#' @param disturbance Probability that a site is disturbed at each time step, losing all its
#'   populations unless `disturbance_sensitivity` is below 1:
#'   one value (0, the default, for none), or a matrix of sites x steps with the probability of
#'   each site at each step, for disturbances that vary in space and time, such as the
#'   seasonal floods of [river_floods()].
#' @param n_steps Number of time steps, one whole number, 1 or more.
#' @param env_time Optional matrix (sites x `n_steps`) with the environment at
#'   each step, e.g. from [env_over_time()]. If `NULL`, the environment of the
#'   landscape is constant. With `env` (several axes), a list with one such matrix per column
#'   of `env` (matched by name when named), or an array of sites x steps x axes.
#' @param seed_steps Steps at which individuals are added during colonisation, a numeric
#'   vector without `NA` (`numeric(0)` for none).
#' @param seed_mean Mean number of individuals of each species added to each
#'   site at each seeding step (and at the start), one finite number, 0 or more.
#' @param record_steps Optional steps at which abundances are recorded, a numeric vector
#'   without `NA`; steps outside 1 to `n_steps` are dropped.
#' @param immigration How emigrants reach other sites on spatially explicit landscapes:
#'   `"poisson"` (the default, as in the code of Thompson et al. 2020) or `"multinomial"`,
#'   which conserves them and is recommended for low dispersal (see Details). Spatially
#'   implicit landscapes are always multinomial.
#' @param pool Optional matrix (sites x species), logical or 0/1, with the sites where each
#'   species can occur.
#' @param pool_from Optional matrix (sites x species) with the first step at which each species
#'   can occur in each site (1 for the start), for example to open a new site to two species at
#'   different times, as in the colonisation-time experiments of Urban & De Meester (2009).
#'   Combined with `pool`. Until then, seeds of that species do not establish there and the site
#'   is not a destination of its emigrants: those that would go there stay in their site, so
#'   that a closed site is not a sink (sites outside `pool` are).
#' @param recruitment How populations grow at each step: `"beverton_holt"` (the default), the
#'   Beverton-Holt competition of Thompson et al. (2020), or `"lottery"`, recruitment into a
#'   fixed number of microsites per site, as in the individual-based models of Urban & De
#'   Meester (2009) and Leibold et al. (2019). With a lottery, each adult survives the step with
#'   probability `1 - mortality`; adults produce Poisson offspring with mean equal to their
#'   growth rate (fecundity times the fitness of their phenotype, so that selection acts on
#'   recruits); offspring emigrate with probability `dispersal`; and the microsites left vacant
#'   by the survivors in each site are filled by a draw without replacement among the offspring
#'   present there (local and immigrant), so that emigration costs the site of origin little
#'   while offspring outnumber the vacancies. Competition is through the shared microsites, so
#'   `alpha` is not used. One trait at most with `evolution`; not available in
#'   [simulate_evolving_metacommunity()].
#' @param microsites With `recruitment = "lottery"`, the number of microsites of each site (one
#'   value or one per site), whole numbers, 1 or more. `initial` must fit in them; seeds take vacant
#'   microsites by lottery with the residents. In the colonisation history, `immigrants` and the
#'   migrant load count immigrant recruits, and `propagule` the individuals of a new population
#'   after recruitment.
#' @param dispersal_time Optional change of emigration through time: a vector with one
#'   multiplier of `dispersal` per step (probabilities above 1 are set to 1), or a matrix with
#'   one row per step and one column, or one per species, of emigration probabilities. For
#'   example, a dispersal pulse during floods or a breeding migration in some steps. `NULL`
#'   (the default) keeps `dispersal` at every step.
#' @param kernels,kernel_step Optional change of the routes of dispersal through time:
#'   `kernels` is a list of dispersal matrices (sites x sites, e.g. from [dispersal_kernel()]
#'   or the `dispersal` of other landscapes on the same sites; rows are rescaled to sum to 1
#'   and the diagonal is ignored), and `kernel_step` gives, for each step, the index or name of
#'   the matrix used. For example, a matrix of connected ponds for flood steps and one of
#'   isolated ponds for dry steps. Only for spatially explicit landscapes; `NULL` (the
#'   default) keeps the dispersal matrix of `landscape`.
#' @param kernel_group Optional routes of dispersal by group of species: one group per species,
#'   given as a column name or index of `kernel_step`, which is then a matrix with one row per
#'   step and one column per group, holding the index or name of the matrix of `kernels` used
#'   by that group at that step. For example, long-distance migrants and sedentary fish of the
#'   same river, each with a matrix for high and one for low water. A run in which every group
#'   uses the same matrix gives the same results as a vector `kernel_step`. With `evolution`
#'   (one or several traits, with or without immigrant selection), gene flow follows the
#'   routes of each species. `NULL` (the default) for one route for all species.
#' @param dispersal_site Optional change of emigration among sites: a matrix with one row per
#'   site and one column per step of multipliers of the emigration probability (probabilities
#'   above 1 are set to 1), applied on top of `dispersal` and `dispersal_time`. For example, a
#'   dispersal pulse only in the floodplain at high water. `NULL` (the default) for the same
#'   emigration at every site.
#' @param disturbance_sensitivity Probability that a population is lost when its site is disturbed,
#'   between 0 and 1: one value or one per species. 1 (the default) removes all the populations of a
#'   disturbed site, as without this argument; lower values let species that resist the disturbance
#'   (for instance drought-tolerant or air-breathing fish when a lake dries out) keep their
#'   populations.
#' @param mortality With `recruitment = "lottery"`, the probability that an adult dies at each step,
#'   whatever its phenotype, one number greater than 0 and up to 1: 1 (the default) for discrete
#'   generations, lower values for overlapping generations, with a generation time of about
#'   `1 / mortality` steps. With `evolution`, selection, mutation and the loss of variance to
#'   selection act on recruits, so that their effect per step scales with the share of recruits
#'   in the population; drift comes from the sampling of both survivors and recruits, so that it
#'   approaches that of a Moran model at low mortality; founders are sampled once, when they
#'   recruit; `steps_per_generation` must stay 1.
#' @param initial Optional matrix (sites x species) of abundances or 0/1 presences at the start,
#'   only where `pool` allows.
#' @param evolution Optional settings from [evolution_control()] for the evolution of the niche
#'   optimum within each population: selection by the local environment, drift, gene flow,
#'   founder sampling and plasticity. `NULL` (the default) keeps the optimum of each species in
#'   all its populations, as before. With a niche of several axes, one trait evolves per axis,
#'   with covariance matrices \eqn{G}{G} and \eqn{E}{E} (see the section Several traits of
#'   [evolution_control()] and [genetic_matrices()]); a correlation between the axes of the
#'   niche ([niche_correlation()]) enters growth and selection.
#' @return An `mc_sim` object: a list with `abundance` (sites x species matrix
#'   at the last step), `history` (sites x species x recorded steps array, or
#'   `NULL`), and the `landscape`, `traits`, `alpha`, `dispersal`, `extirpation` and `pool`
#'   used, `env` when given, and `tree` when the traits came from [species_traits()] with a
#'   tree, the colonisation history (`arrival`, `propagule`, `immigrants`, `migrant_load`;
#'   see the section Colonisation history), `n_steps`, and the settings `dispersal_time`,
#'   `kernels` (rescaled), `kernel_step`, `kernel_group`, `dispersal_site`,
#'   `disturbance_sensitivity` (one per species), `disturbance`, `recruitment`, `microsites` and
#'   `mortality`. With `evolution`, also `evolution`, a
#'   list with the genetic mean of the optimum in each population (`genetic`), the plastic part
#'   of its phenotype at the last step (`plastic`), the phenotype after dispersal (`phenotype`,
#'   which mixes residents and immigrants), the genetic means and the environment at
#'   `record_steps` (`history`, `environment_history`), the steps since the environment of each
#'   site last changed on any axis (`since_change`, a list with `end`, one value per site at
#'   the last step, and `history`, sites x recorded steps; 0 in the step of a change and the
#'   step itself where the environment never changed), the reference environment and the
#'   settings, and, with `mutation_var`, the genetic variance of each population at the last
#'   step (`genetic_variance`); `NA` where a population is absent. With `loci`,
#'   `allele_frequency` is an array of sites x species x loci with the frequency of allele 1 at
#'   each locus in each population at the last step, `NA` where a population is absent, and
#'   `genetic_variance` is the variance that follows from them. With several evolving traits,
#'   `genetic`,
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
                                   pool_from = NULL,
                                   recruitment = c("beverton_holt", "lottery"),
                                   microsites = NULL, mortality = 1,
                                   dispersal_time = NULL, kernels = NULL,
                                   kernel_step = NULL, kernel_group = NULL,
                                   dispersal_site = NULL, disturbance_sensitivity = 1) {
  cite_method("simulate_metacommunity")
  if (!inherits(landscape, "mc_landscape")) {
    stop("`landscape` must be an `mc_landscape` object, such as one from ",
         "`landscape_implicit()` or `landscape_grid()`", call. = FALSE)
  }
  check_count(n_steps, "n_steps")
  immigration <- match.arg(immigration)
  recruitment <- match.arg(recruitment)
  inp <- mc_sim_inputs(landscape, traits, alpha, env, env_time, evolution, n_steps)
  out <- mc_dynamics(inp$r, inp$alpha, dispersal, n_steps = n_steps, seed_steps = seed_steps,
                     seed_mean = seed_mean, disp_matrix = landscape$dispersal,
                     record_steps = record_steps, sites = inp$sites,
                     species = rownames(inp$traits), extirpation = extirpation, pool = pool,
                     initial = initial, immigration = immigration,
                     evolution = inp$evolution, disturbance = disturbance,
                     pool_from = pool_from, recruitment = recruitment,
                     microsites = microsites, mortality = mortality,
                     dispersal_time = dispersal_time, kernels = kernels,
                     kernel_step = kernel_step, kernel_group = kernel_group,
                     dispersal_site = dispersal_site,
                     disturbance_sensitivity = disturbance_sensitivity)
  structure(c(out, list(landscape = landscape, env = inp$env, traits = inp$niche,
                        alpha = inp$alpha, dispersal = dispersal,
                        dispersal_time = dispersal_time, kernels = kernels,
                        kernel_step = kernel_step, kernel_group = kernel_group,
                        dispersal_site = dispersal_site,
                        disturbance_sensitivity = disturbance_sensitivity,
                        extirpation = extirpation, disturbance = disturbance, pool = out$pool,
                        tree = attr(inp$traits, "tree"), recruitment = recruitment,
                        microsites = if (recruitment == "lottery") microsites,
                        mortality = if (recruitment == "lottery") mortality)),
            class = "mc_sim")
}

#' the inputs of mc_dynamics() from those of simulate_metacommunity(): the growth rates `r`
#' (a matrix, or a function of the step), `alpha` in the order of the species, the names of the
#' `sites`, the environment `env` of a niche of several axes, the `traits` of the first axis,
#' the whole `niche` (the traits, or a list of them) and the settings of `evolution`
#' @noRd
mc_sim_inputs <- function(landscape, traits, alpha, env, env_time, evolution, n_steps) {
  several <- is.list(traits) && !is.data.frame(traits)
  niche <- traits
  if (several || !is.null(env)) {
    axes <- mc_niche_axes(traits, env, landscape, env_time, n_steps, several)
    env <- axes$env
    env_time <- axes$env_time
    traits <- axes$trait_list[[1]]
    niche <- if (several) axes$trait_list else traits
  }
  alpha <- mc_alpha_species(alpha, traits)
  sites <- paste0("site", seq_along(landscape$env))
  if (!is.null(env)) rownames(env) <- sites
  r <- mc_growth(landscape, traits, niche, env, env_time, n_steps, sites)
  if (!is.null(evolution)) {
    evolution <- mc_evolution_inputs(evolution, several, traits, niche, env, env_time,
                                     landscape)
  }
  list(r = r, alpha = alpha, sites = sites, env = env, traits = traits, niche = niche,
       evolution = evolution)
}

#' the niche of several axes: the tables of traits (one per axis), the environment of the
#' sites on each axis as a data frame and, with `env_time`, its matrices through time
#' @noRd
mc_niche_axes <- function(traits, env, landscape, env_time, n_steps, several) {
  trait_list <- if (several) lapply(traits, as_traits) else list(as_traits(traits))
  attr(trait_list, "niche_cor") <- attr(traits, "niche_cor")
  if (is.null(env)) env <- landscape$env_axes
  if (is.null(env)) stop("simulate_metacommunity(): a list of traits needs `env`, one ",
                         "column per element", call. = FALSE)
  unnamed <- is.null(colnames(env)) && !is.data.frame(env)
  env <- as.data.frame(env)
  if (unnamed) {
    names(env) <- (if (several) names(traits)) %||% paste0("axis", seq_len(ncol(env)))
  }
  if (nrow(env) != length(landscape$env) || ncol(env) != length(trait_list)) {
    stop("simulate_metacommunity(): `env` needs one row per site (",
         length(landscape$env), ") and one column per ",
         "set of traits (", length(trait_list), ")", call. = FALSE)
  }
  if (!is.null(env_time)) env_time <- env_time_axes(env_time, env, n_steps)
  list(trait_list = trait_list, env = env, env_time = env_time)
}

#' `alpha` with its rows and columns in the order of the species of `traits`, when it names
#' them; stops unless it has one row and one column per species
#' @noRd
mc_alpha_species <- function(alpha, traits) {
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
  alpha
}

#' growth rates of the species in the sites: a matrix of sites x species or, when the
#' environment changes through time (`env_time`), a function of the step that returns one
#' @noRd
mc_growth <- function(landscape, traits, niche, env, env_time, n_steps, sites) {
  force(sites)
  if (!is.null(env)) {
    if (is.null(env_time)) return(growth_rates(env, niche))
    return(function(step) {
      now <- as.data.frame(lapply(env_time, function(x) x[, step]), row.names = sites)
      growth_rates(now, niche)
    })
  }
  if (is.null(env_time)) return(growth_rates(stats::setNames(landscape$env, sites), traits))
  check_env_time_matrix(env_time, length(landscape$env), n_steps)
  function(step) growth_rates(stats::setNames(env_time[, step], sites), traits)
}

#' stops unless `env_time` is a numeric matrix with one row per site and a column per step
#' @noRd
check_env_time_matrix <- function(env_time, n_sites, n_steps) {
  if (!is.matrix(env_time) || !is.numeric(env_time) ||
        nrow(env_time) != n_sites || ncol(env_time) < n_steps) {
    stop("`env_time` must be a numeric matrix with one row per site (",
         n_sites, ") and at least `n_steps` (", n_steps, ") columns",
         call. = FALSE)
  }
  invisible(env_time)
}

#' the settings of `evolution` with the traits and the environment (fixed, or a function of
#' the step) that the engine needs, for one trait or for one trait per axis of the niche
#' @noRd
mc_evolution_inputs <- function(evolution, several, traits, niche, env, env_time, landscape) {
  check_evolution(evolution)
  mv <- is_covariance(evolution$genetic_var) || is_covariance(evolution$env_var)
  if (several || mv) {
    if (!several) {
      stop("simulate_metacommunity(): `evolution` with covariance matrices needs a niche of ",
           "several axes (a list of ",
           "`traits` and `env`), one trait per axis", call. = FALSE)
    }
    evolution$traits <- stats::setNames(niche, names(env))
    evolution$env <- if (is.null(env_time)) as.matrix(env) else {
      function(step) vapply(env_time, function(x) x[, step], numeric(nrow(env)))
    }
    return(evolution)
  }
  evolution$traits <- traits
  evolution$env <- if (!is.null(env)) {
    if (is.null(env_time)) env[[1]] else function(step) env_time[[1]][, step]
  } else if (is.null(env_time)) landscape$env else function(step) env_time[, step]
  evolution
}

#' the environment through time on several axes: a list of matrices (sites x steps), one per
#' column of `env` (matched by name when named), or an array of sites x steps x axes
#' @noRd
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
  check_env_time(env_time, env, n_steps)
  stats::setNames(lapply(env_time, unname), names(env))
}

#' stops unless each matrix of the list `env_time` has one row per site and a column per step
#' @noRd
check_env_time <- function(env_time, env, n_steps) {
  ok <- vapply(env_time, function(x) {
    is.numeric(x) && is.matrix(x) && nrow(x) == nrow(env) && ncol(x) >= n_steps && !anyNA(x)
  }, logical(1))
  if (!all(ok)) {
    stop("each matrix of `env_time` needs one row per site (", nrow(env), ") and at least ",
         "`n_steps` (", n_steps, ") columns, without NA", call. = FALSE)
  }
  invisible(env_time)
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
                        evolution = NULL, disturbance = 0, pool_from = NULL,
                        recruitment = c("beverton_holt", "lottery"), microsites = NULL,
                        mortality = 1, dispersal_time = NULL, kernels = NULL,
                        kernel_step = NULL, kernel_group = NULL, dispersal_site = NULL,
                        disturbance_sensitivity = 1) {
  # the state lives in this frame; the preparation, each step and the end are evaluated in it, so
  # that simulate_evolving_metacommunity() can run the same steps and change species in between
  mc_setup(environment())
  for (step in seq_len(n_steps)) mc_step(environment(), step)
  mc_finish(environment())
}

#' preparation of mc_dynamics(), evaluated in the frame that holds its arguments: each part
#' below is evaluated in that frame, in this order, so that they share its variables
#' @noRd
mc_setup <- function(st) {
  mc_setup_arguments(st)
  mc_setup_dispersal(st)
  mc_setup_start(st)
  mc_setup_abundance(st)
  mc_setup_evolution(st)
  eval(quote({
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

#' checks the arguments of mc_dynamics() on steps, growth rates, competition, dispersal,
#' extirpation, disturbance and recruitment, in the frame `st`, and sets the numbers of sites
#' (`m`) and species (`k_all`)
#' @noRd
mc_setup_arguments <- function(st) {
  eval(quote({
    immigration <- match.arg(immigration)
    check_count(n_steps, "n_steps")
    check_number(seed_mean, "seed_mean", lower = 0)
    mc_check_steps(seed_steps, record_steps)
    r_fun <- if (is.function(r)) r else function(step) r
    r1 <- r_fun(1)
    if (!is.matrix(r1) || !is.numeric(r1)) {
      stop("`r` must be a numeric matrix of sites x species, or a function of the step that ",
           "returns one", call. = FALSE)
    }
    m <- nrow(r1)
    k_all <- ncol(r1)
    mc_check_matrices(alpha, disp_matrix, m, k_all)
    probabilities_per_species(dispersal, "dispersal", k_all)
    probabilities_per_species(extirpation, "extirpation", k_all)
    dispersal <- rep_len(dispersal, k_all)
    extirpation <- rep_len(extirpation, k_all)
    mc_check_disturbance(disturbance, m, n_steps)
    recruitment <- match.arg(recruitment, c("beverton_holt", "lottery"))
    lottery <- recruitment == "lottery"
    if (lottery) {
      microsites <- mc_lottery_microsites(microsites, m)
      mc_check_mortality(mortality)
    }
    if (m == 1 && any(dispersal > 0)) {
      stop("dispersal needs at least two sites: use `dispersal = 0` with one site", call. = FALSE)
    }
  }), st)
  invisible(st)
}

#' stops unless `seed_steps` and `record_steps` hold numbers of steps without NA
#' @noRd
mc_check_steps <- function(seed_steps, record_steps) {
  if (!is.numeric(seed_steps) || anyNA(seed_steps)) {
    stop("`seed_steps` must be numbers of steps without NA (`numeric(0)` for none)",
         call. = FALSE)
  }
  if (!is.null(record_steps) && (!is.numeric(record_steps) || anyNA(record_steps))) {
    stop("`record_steps` must be `NULL` or numbers of steps without NA", call. = FALSE)
  }
  invisible(seed_steps)
}

#' stops unless `alpha` has one row and one column per species (`k_all`) and `disp_matrix`,
#' when given, one row and one column per site (`m`)
#' @noRd
mc_check_matrices <- function(alpha, disp_matrix, m, k_all) {
  # one species may come with a single number, as `alpha[sp, sp]` gives
  one_value <- k_all == 1 && is.null(dim(alpha)) && length(alpha) == 1
  if (!is.numeric(alpha) ||
        !(one_value || identical(as.integer(dim(alpha)), c(k_all, k_all)))) {
    stop("`alpha` must be a numeric matrix with one row and one column per species (",
         k_all, ")", call. = FALSE)
  }
  if (!is.null(disp_matrix) && !identical(as.integer(dim(disp_matrix)), c(m, m))) {
    stop("`disp_matrix` must be a matrix with one row and one column per site (", m, ")",
         call. = FALSE)
  }
  invisible(alpha)
}

#' stops unless `disturbance` is one probability, or a matrix of probabilities with one row
#' per site and one column per step
#' @noRd
mc_check_disturbance <- function(disturbance, m, n_steps) {
  # one value, or one row per site and one column per step
  shaped <- function() {
    by_step <- is.matrix(disturbance) && nrow(disturbance) == m && ncol(disturbance) >= n_steps
    length(disturbance) == 1 || by_step
  }
  bad <- !is.numeric(disturbance) || anyNA(disturbance) || any(disturbance < 0)
  if (bad || any(disturbance > 1) || !shaped()) {
    stop("`disturbance` must be one probability between 0 and 1, or a matrix of ",
         "probabilities with one row per site and one column per step (e.g. from ",
         "river_floods())", call. = FALSE)
  }
  invisible(disturbance)
}

#' the microsites of each of the `m` sites for a lottery, as integers, from one value or one
#' per site; stops unless they are whole numbers of 1 or more
#' @noRd
mc_lottery_microsites <- function(microsites, m) {
  bad <- is.null(microsites) || !is.numeric(microsites) || !length(microsites) %in% c(1, m)
  bad <- bad || anyNA(microsites) || any(microsites < 1)
  bad <- bad || any(microsites != round(microsites))
  if (bad || any(microsites > .Machine$integer.max)) {
    stop("`recruitment = \"lottery\"` needs `microsites`: whole numbers of 1 or more, one ",
         "value or one per site (", m, ")", call. = FALSE)
  }
  rep_len(as.integer(microsites), m)
}

#' stops unless `mortality` is one probability above 0 and up to 1
#' @noRd
mc_check_mortality <- function(mortality) {
  bad <- !is.numeric(mortality) || length(mortality) != 1 || is.na(mortality)
  if (bad || mortality <= 0 || mortality > 1) {
    stop("`mortality` must be one probability above 0 and up to 1", call. = FALSE)
  }
  invisible(mortality)
}

#' prepares, in the frame `st`, dispersal through time and by kernel, the multipliers of
#' emigration by site, the sensitivity to disturbance and the names of sites and species
#' @noRd
mc_setup_dispersal <- function(st) {
  eval(quote({
    # emigration probabilities at each step (steps x species), from `dispersal_time`
    disp_prob_time <- dispersal_by_step(dispersal_time, dispersal, n_steps, k_all)
    if (!is.null(disp_prob_time) && m == 1 && any(disp_prob_time > 0)) {
      stop("`dispersal_time` needs at least two sites", call. = FALSE)
    }
    kernels <- kernels_by_step(kernels, kernel_step, n_steps, m, disp_matrix, kernel_group)
    kernel_group <- kernel_groups(kernel_group, kernel_step, k_all)
    if (!is.null(dispersal_site)) mc_check_dispersal_site(dispersal_site, m, n_steps)
    disturbance_sensitivity <- mc_sensitivity(disturbance_sensitivity, k_all)
    sites <- sites %||% rownames(r1) %||% paste0("site", seq_len(m))
    species <- species %||% colnames(r1) %||% paste0("sp", seq_len(k_all))
  }), st)
  invisible(st)
}

#' stops unless `dispersal_site` holds multipliers of emigration by site and step, in a
#' landscape of two or more sites
#' @noRd
mc_check_dispersal_site <- function(dispersal_site, m, n_steps) {
  bad <- !is.matrix(dispersal_site) || !is.numeric(dispersal_site)
  bad <- bad || nrow(dispersal_site) != m || ncol(dispersal_site) < n_steps
  bad <- bad || anyNA(dispersal_site) || any(dispersal_site < 0)
  if (bad || !all(is.finite(dispersal_site))) {
    stop("`dispersal_site` must be a matrix of multipliers of emigration (0 or more) with ",
         "one row per site (", m, ") and one column per step (", n_steps, ")",
         call. = FALSE)
  }
  if (m == 1) stop("`dispersal_site` needs at least two sites", call. = FALSE)
  invisible(dispersal_site)
}

#' the sensitivity to disturbance of each of the `k_all` species, from one probability or one
#' per species
#' @noRd
mc_sensitivity <- function(disturbance_sensitivity, k_all) {
  bad <- !is.numeric(disturbance_sensitivity) || anyNA(disturbance_sensitivity)
  bad <- bad || any(disturbance_sensitivity < 0) || any(disturbance_sensitivity > 1)
  if (bad || !length(disturbance_sensitivity) %in% c(1, k_all)) {
    stop("`disturbance_sensitivity` must be one probability between 0 and 1, or one per ",
         "species (", k_all, ")", call. = FALSE)
  }
  rep_len(as.numeric(disturbance_sensitivity), k_all)
}

#' prepares, in the frame `st`, the pools of the species, the records and the species alive
#' @noRd
mc_setup_start <- function(st) {
  eval(quote({
    if (!is.null(pool)) pool <- mc_pool(pool, sites, species)
    if (!is.null(initial)) initial <- mc_initial(initial, pool, sites, species)
    record_steps <- sort(unique(record_steps[record_steps >= 1 & record_steps <= n_steps]))
    history <- if (length(record_steps))
      array(0L, c(m, k_all, length(record_steps)),
            dimnames = list(sites, species, step_label(record_steps))) else NULL

    alive <- seq_len(k_all)
    k <- k_all
    p_em <- rep(dispersal, each = m)
    p_ex <- rep(extirpation, each = m)
    at <- t(alpha)
    if (!is.null(pool_from)) pool_from <- mc_pool_from(pool_from, initial, sites, species)
    # sites open to each species at a step: `pool` and, with `pool_from`, from that step on
    inside_at <- function(step, cols) {
      x <- if (is.null(pool)) matrix(TRUE, m, length(cols)) else pool[, cols, drop = FALSE]
      if (!is.null(pool_from)) x <- x & pool_from[, cols, drop = FALSE] <= step
      x
    }
    inside <- if (is.null(pool) && is.null(pool_from)) NULL else inside_at(1, seq_len(k_all))
  }), st)
  invisible(st)
}

#' draws, in the frame `st`, the abundances at the start (`n`), unless `initial` gives them
#' @noRd
mc_setup_abundance <- function(st) {
  eval(quote({
    n <- if (is.null(initial)) matrix(stats::rpois(m * k, seed_mean), m, k) else unname(initial)
    if (!is.null(inside)) n[!inside] <- 0L
    if (lottery) {
      if (!is.null(initial) && any(rowSums(n) > microsites)) {
        stop("`initial` has more individuals than `microsites` in some sites", call. = FALSE)
      }
      # seeds at the start take the microsites by lottery
      n <- lottery_draw(n, microsites)
    }
    last_seed <- max(c(0, seed_steps))
    competition_hook <- NULL
  }), st)
  invisible(st)
}

#' the pool of each species as a logical matrix of sites x species
#' @noRd
mc_pool <- function(pool, sites, species) {
  pool <- site_species(pool, sites, species, "`pool`")
  if (!all(pool %in% c(0, 1))) stop("`pool` must be logical or 0/1", call. = FALSE)
  pool == 1
}

#' the abundances at the start as an integer matrix of sites x species; stops unless they are
#' whole numbers, zero or more, inside `pool`
#' @noRd
mc_initial <- function(initial, pool, sites, species) {
  initial <- site_species(initial, sites, species, "`initial`")
  if (any(initial < 0) || any(initial != round(initial))) {
    stop("`initial` must hold whole numbers of individuals, zero or more", call. = FALSE)
  }
  if (!is.null(pool) && any(initial[!pool] > 0)) {
    stop("`initial` has individuals outside `pool`", call. = FALSE)
  }
  storage.mode(initial) <- "integer"
  initial
}

#' the step from which each site opens to each species, as a matrix of sites x species; stops
#' when `initial` has individuals where the site opens after the start
#' @noRd
mc_pool_from <- function(pool_from, initial, sites, species) {
  pool_from <- site_species(pool_from, sites, species, "`pool_from`")
  if (anyNA(pool_from) || any(pool_from < 1)) {
    stop("`pool_from` must hold steps of 1 or more", call. = FALSE)
  }
  if (!is.null(initial) && any(initial[pool_from > 1] > 0)) {
    stop("`initial` has individuals where `pool_from` opens after the start", call. = FALSE)
  }
  pool_from
}

#' prepares, in the frame `st`, the evolutionary state (`ev`) and the genetic means of the
#' populations at the start
#' @noRd
mc_setup_evolution <- function(st) {
  eval(quote({
    ev <- if (!is.null(evolution)) {
      check_evolution(evolution)
      if (evolution_is_mv(evolution)) evolution_state_mv(evolution, m, k_all, species) else
        evolution_state(evolution, m, k_all, species)
    }
    mc_check_evolution_state(ev, lottery, disp_matrix, immigration)
    if (isTRUE(ev$mv)) {
      mc_setup_traits_mv(environment())
    } else if (!is.null(ev)) {
      mc_setup_traits(environment())
    }
  }), st)
  invisible(st)
}

#' stops when the evolutionary state `ev` cannot go with the recruitment or the immigration
#' @noRd
mc_check_evolution_state <- function(ev, lottery, disp_matrix, immigration) {
  if (lottery) mc_check_lottery_evolution(ev)
  if (isTRUE(ev$imm_sel) && (is.null(disp_matrix) || immigration != "poisson")) {
    stop("`immigrant_selection` needs a spatially explicit landscape and ",
         "`immigration = \"poisson\"`", call. = FALSE)
  }
  invisible(ev)
}

#' stops when the evolutionary state `ev` cannot go with recruitment by lottery
#' @noRd
mc_check_lottery_evolution <- function(ev) {
  if (!is.null(ev) && !isTRUE(ev$mv) && ev$per != 1) {
    stop("with `recruitment = \"lottery\"`, keep `steps_per_generation = 1`: overlapping ",
         "generations come from `mortality`", call. = FALSE)
  }
  if (isTRUE(ev$mv)) {
    stop("`recruitment = \"lottery\"` is not yet available with several traits", call. = FALSE)
  }
  invisible(ev)
}

#' the genetic means and phenotypes at the start with several traits, in the frame `st`
#' @noRd
mc_setup_traits_mv <- function(st) {
  eval(quote({
    # several traits: genetic means and phenotypes as sites x species x traits
    z <- ph <- array(NA_real_, c(m, k_all, ev$q))
    z_history <- if (length(record_steps))
      array(NA_real_, c(m, k_all, ev$q, length(record_steps)),
            dimnames = list(sites, species, ev$axes, step_label(record_steps))) else NULL
    fixed_z <- ev$fixed
    pool_mean <- array(rep(ev$optimum, each = m), c(m, k_all, ev$q))
    start_pop <- n > 0
    if (any(start_pop)) z <- mv_found(ev, z, pool_mean, n, start_pop, seq_len(k_all))
  }), st)
  invisible(st)
}

#' the genetic means, phenotypes and genetic variances at the start with one trait, and the
#' function `found()` that draws the mean of a population from its founders, in the frame `st`
#' @noRd
mc_setup_traits <- function(st) {
  eval(quote({
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
    vd <- NULL
    if (isTRUE(ev$dyn)) {
      vg <- matrix(NA_real_, m, k_all)
      vg[start_pop] <- pool_g[start_pop]
    }
    # linkage disequilibrium of every population, with variance_selection = "sexual"
    if (isTRUE(ev$sexual)) {
      vd <- matrix(NA_real_, m, k_all)
      vd[start_pop] <- 0
    }
    # allele frequencies of every population, with `loci`
    if (isTRUE(ev$loci_on)) loci_setup(environment())
  }), st)
  invisible(st)
}

#' stops unless `evolution` is a set of settings from evolution_control()
#' @noRd
check_evolution <- function(evolution) {
  if (!inherits(evolution, "mc_evolution")) {
    stop("`evolution` must come from `evolution_control()`", call. = FALSE)
  }
  invisible(evolution)
}

#' stops unless `x` holds probabilities between 0 and 1, one value or one per species (`k`)
#' @noRd
probabilities_per_species <- function(x, name, k) {
  if (!is.numeric(x) || anyNA(x) || any(x < 0 | x > 1) || !length(x) %in% c(1, k)) {
    stop("`", name, "` must be one probability between 0 and 1, or one per species (", k, ")",
         call. = FALSE)
  }
  invisible(x)
}

#' one time step of mc_dynamics(), evaluated in the state `st`: each phase below is evaluated
#' in that frame, in this order, so that they share its variables
#' @noRd
mc_step <- function(st, step) {
  assign("step", step, envir = st)
  mc_step_seeds(st)
  mc_step_prune(st)
  if (st$k > 0) {
    mc_step_growth(st)
    mc_step_reproduce(st)
    mc_step_select(st)
    mc_step_emigrate(st)
    if (any(st$tot > 0)) {
      mc_step_immigrate(st)
      mc_step_gene_flow(st)
      st$n <- st$n - st$em + st$imm
    }
    if (st$lottery) mc_step_recruit(st)
    mc_step_losses(st)
    mc_step_store(st)
    mc_step_colonisation(st)
  }
  mc_step_record(st)
  invisible(st)
}

#' seeds from the regional pool, at the steps of `seed_steps`: sets `added` (0 without seeds)
#' and adds them to `n` in the state `st`
#' @noRd
mc_step_seeds <- function(st) {
  eval(quote({
    added <- 0L
    if (!is.null(pool_from)) inside <- inside_at(step, alive)
    if (step %in% seed_steps) {
      seeds <- stats::rpois(m * k, seed_mean)
      added <- if (is.null(inside)) seeds else seeds * inside
      if (isTRUE(ev$mv)) {
        mc_seed_traits_mv(environment())
      } else if (!is.null(ev)) {
        mc_seed_traits(environment())
      }
      n <- n + added
      # with a lottery, seeds compete with residents for the microsites
      if (lottery) n <- lottery_draw(n, microsites)
    }
  }), st)
  invisible(st)
}

#' genetic means after the seeds of this step, with several traits: seeds carry the mean of
#' the regional pool and mix with residents; changes `z` in the state `st`
#' @noRd
mc_seed_traits_mv <- function(st) {
  eval(quote({
    added_m <- added + matrix(0L, m, k)
    za <- z[, alive, , drop = FALSE]
    pm <- pool_mean[, alive, , drop = FALSE]
    mix <- n > 0 & added_m > 0
    w_mix <- ifelse(mix, n / (n + added_m), 1)
    for (a in seq_len(ev$q)) za[, , a] <- w_mix * za[, , a] + (1 - w_mix) * pm[, , a]
    new_seed <- n == 0 & added_m > 0
    if (any(new_seed)) za <- mv_found(ev, za, pm, added_m, new_seed, alive)
    z[, alive, ] <- za
  }), st)
  invisible(st)
}

#' genetic means (and variances) after the seeds of this step, with one trait; changes `z` and
#' `vg` in the state `st`
#' @noRd
mc_seed_traits <- function(st) {
  eval(quote({
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
      if (isTRUE(ev$sexual)) {
        # seeds carry no disequilibrium
        dv <- vd[, alive, drop = FALSE]
        dv[mix] <- w * dv[mix]
        dv[new_seed] <- 0
        vd[, alive] <- dv
      }
    }
    if (isTRUE(ev$loci_on)) loci_step_seed(environment())
  }), st)
  invisible(st)
}

#' drops from the state `st` the species with no individual left, once seeding is over
#' @noRd
mc_step_prune <- function(st) {
  eval(quote({
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
  }), st)
  invisible(st)
}

#' growth rates of this step (`r_now`) and, with evolution, the genetic means (`za`), the
#' phenotypes (`pa`) and the environment (`e_now`) of the alive species, in the state `st`
#' @noRd
mc_step_growth <- function(st) {
  eval(quote({
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
      if (any(fixed_z[alive])) {
        za[, fixed_z[alive]] <- pool_mean[, alive, drop = FALSE][, fixed_z[alive]]
      }
      if (isTRUE(ev$loci_on)) loci_step_start(environment())
      pa <- za + sweep(matrix(e_now - ev$ref, m, k), 2, ev$b[alive], "*")
      if (isTRUE(ev$dyn)) {
        ga <- vg[, alive, drop = FALSE]
        ga[is.na(ga)] <- pool_g[, alive, drop = FALSE][is.na(ga)]
        if (isTRUE(ev$sexual)) {
          da <- vd[, alive, drop = FALSE]
          da[is.na(da)] <- 0
        }
        r_now <- evolution_growth_cells(pa, e_now, ev$sigma[alive], ev$r_max[alive],
                                        sweep(ga, 2, ev$ve[alive], "+"))
      } else {
        r_now <- evolution_growth(pa, e_now, ev$sigma[alive], ev$r_max[alive],
                                  ev$g[alive] + ev$ve[alive])
      }
    }
  }), st)
  invisible(st)
}

#' reproduction: the offspring `n` of this step, from the abundances at its start (`start`);
#' with a lottery, also the adults that survive (`surv`) and their genetic means
#' @noRd
mc_step_reproduce <- function(st) {
  eval(quote({
    start <- n
    if (lottery) {
      # adults that survive the step, whatever their phenotype; n becomes the offspring
      surv <- matrix(stats::rbinom(m * k, n, 1 - mortality), m, k)
      # extra sampling variance and cells without drift for populations founded this step
      extra_v <- matrix(0, m, k)
      no_drift <- matrix(FALSE, m, k)
      if (!is.null(ev)) {
        z_surv <- za
        if (isTRUE(ev$dyn)) g_surv <- ga
        if (isTRUE(ev$sexual)) d_surv <- da
        # survivors are a sample of the adults: drift of their mean and loss of variance, so
        # that drift per generation approaches that of a Moran model at low mortality
        g0 <- if (isTRUE(ev$dyn)) ga else matrix(ev$g[alive], m, k, byrow = TRUE)
        x <- ifelse(surv > 0 & surv < n, 1 / pmax(surv, 1) - 1 / pmax(n, 1), 0)
        sd_s <- surv > 0 & x > 0 & g0 > 0
        if (isTRUE(ev$loci_on)) {
          loci_step_survive(environment())
        } else if (is.finite(ev$nu) && any(sd_s)) {
          z_surv[sd_s] <- z_surv[sd_s] + stats::rnorm(sum(sd_s), 0,
                                                      sqrt(g0[sd_s] * x[sd_s] / ev$nu))
          if (isTRUE(ev$dyn)) g_surv[sd_s] <- (g_surv * pmax(1 - x / (2 * ev$nu), 0))[sd_s]
          if (isTRUE(ev$sexual)) d_surv[sd_s] <- (d_surv * pmax(1 - x / (2 * ev$nu), 0))[sd_s]
        }
      }
      lam <- r_now * n
    } else {
      # crowding; simulate_evolving_metacommunity() may replace it with competition that
      # depends on the phenotypes in each site
      crowding <- if (is.null(competition_hook)) n %*% at else competition_hook(n, pa)
      lam <- r_now * n / (1 + crowding)
    }
    n <- matrix(stats::rpois(m * k, lam), m, k)
  }), st)
  invisible(st)
}

#' selection toward the local environment and drift in the populations present after
#' reproduction; changes `za`, `pa` and the genetic variances in the state `st`
#' @noRd
mc_step_select <- function(st) {
  eval(quote({
    if (isTRUE(ev$mv)) {
      za <- mv_select(ev, za, pa, e_now, n, alive)
      pa <- mv_phenotype(ev, za, e_now, alive)
    } else if (isTRUE(ev$loci_on)) {
      loci_step_select(environment())
    } else if (!is.null(ev)) {
      # selection toward the local environment and drift, in populations present after
      # reproduction (Lande 1976); no draws without drift
      g <- if (isTRUE(ev$dyn)) ga else matrix(ev$g[alive], m, k, byrow = TRUE)
      width <- g + matrix(ev$sigma[alive]^2 + ev$ve[alive], m, k, byrow = TRUE)
      live <- n > 0
      za[live] <- za[live] + (g[live] / width[live]) * (e_now - pa)[live] / ev$per
      drift <- live & g > 0
      # with a lottery, drift and the dynamics of the variance act on the recruits
      if (!lottery && is.finite(ev$nu) && any(drift)) {
        za[drift] <- za[drift] +
          stats::rnorm(sum(drift), 0, sqrt(g[drift] / (ev$per * ev$nu * n[drift])))
      }
      za[!live] <- NA_real_
      if (isTRUE(ev$dyn) && !lottery) {
        # mutation-drift dynamics of the genetic variance (Lynch & Hill 1986): a loss of
        # 1 / (2 Ne) and a gain of Vm per generation, equilibrium 2 Ne Vm; with
        # `variance_selection`, also the loss to stabilising selection
        mc_variance_step(environment(), live, n, width, 1 / ev$per)
      }
      # phenotype of the offspring, which develop in their site
      pa <- za + sweep(matrix(e_now - ev$ref, m, k), 2, ev$b[alive], "*")
    }
  }), st)
  invisible(st)
}

#' emigration: the emigrants of each population (`em`), those that count toward the flows into
#' open sites (`em_to`) and their totals by species (`tot`); `imm` starts at 0
#' @noRd
mc_step_emigrate <- function(st) {
  eval(quote({
    # phenotypes of the emigrants, for the fitness of immigrants in the speciation engine
    if (isTRUE(ev$imm_sel)) pa_em <- pa
    if (!is.null(disp_prob_time)) p_em <- rep(disp_prob_time[step, alive], each = m)
    p_go <- if (is.null(dispersal_site)) p_em else
      pmin(p_em * rep(dispersal_site[, step], k), 1)
    if (!is.null(kernels)) {
      disp_matrix <- if (is.null(kernel_group)) kernels[[kernel_step[step]]] else
        structure(kernels[kernel_step[step, ]], group = kernel_group[alive])
    }
    if (is.null(pool_from)) {
      em <- em_to <- matrix(stats::rbinom(m * k, n, p_go), m, k)
      open_to <- NULL
    } else {
      # sites not yet open through `pool_from` are not destinations: the emigrants that would
      # go there stay home, and the flows into open sites are unchanged
      open_to <- pool_from[, alive, drop = FALSE] <= step
      share <- if (is.null(disp_matrix)) {
        (matrix(colSums(open_to), m, k, byrow = TRUE) - open_to) / (m - 1)
      } else kernel_product(disp_matrix, open_to)
      em <- matrix(stats::rbinom(m * k, n, p_go * share), m, k)
      em_to <- ifelse(share > 0, em / share, 0)
    }
    tot <- colSums(em)
    imm <- 0L
  }), st)
  invisible(st)
}

#' immigration: the immigrants `imm` of each population from the emigrants of this step, by
#' the kind of landscape and of immigration, in the state `st`
#' @noRd
mc_step_immigrate <- function(st) {
  eval(quote({
    if (is.null(disp_matrix)) {
      # an implicit landscape: the emigrants of the other sites, pooled
      w <- sweep(-em, 2, tot, "+")
      if (!is.null(open_to)) w[!open_to] <- 0
      imm <- matrix(0L, m, k)
      for (j in which(tot > 0)) imm[, j] <- stats::rmultinom(1, tot[j], w[, j])
    } else if (immigration == "multinomial") {
      imm <- place_emigrants(em, disp_matrix, open_to)
    } else if (isTRUE(ev$imm_sel) && isTRUE(ev$mv)) {
      # immigrants filtered by their fitness in the site they reach, several traits
      sel <- mv_immigrants_selected(ev, em_to, za, pa, e_now, disp_matrix, alive)
      imm <- matrix(stats::rpois(m * k, sel$lambda), m, k)
    } else if (isTRUE(ev$imm_sel)) {
      # immigrants filtered by their fitness in the site they reach
      g_src <- if (isTRUE(ev$dyn)) ifelse(is.na(ga), 0, ga) else
        matrix(ev$g[alive], m, k, byrow = TRUE)
      sel <- immigrants_selected(em_to, ifelse(is.na(za), 0, za), ifelse(is.na(pa), 0, pa),
                                 e_now, g_src + matrix(ev$sigma[alive]^2 + ev$ve[alive], m,
                                                       k, byrow = TRUE),
                                 disp_matrix, g_src)
      imm <- matrix(stats::rpois(m * k, sel$lambda), m, k)
    } else {
      imm <- matrix(stats::rpois(m * k, kernel_flow(disp_matrix, em_to)), m, k)
    }
  }), st)
  invisible(st)
}

#' gene flow: residents and immigrants mixed by their numbers, in the state `st`
#' @noRd
mc_step_gene_flow <- function(st) {
  if (isTRUE(st$ev$mv)) {
    eval(quote({
      flowed <- mv_flow(ev, za, pa, n - em + em_to, em_to, imm, disp_matrix, alive,
                        sel = if (isTRUE(ev$imm_sel) && !is.null(disp_matrix)) sel)
      za <- flowed$za
      pa <- flowed$pa
    }), st)
  } else if (!is.null(st$ev)) {
    mc_flow_sources(st)
    if (isTRUE(st$ev$loci_on)) loci_step_flow(st) else mc_flow_merge(st)
  }
  invisible(st)
}

#' what the immigrants of each population carry, with one trait: their mean (`zi`), the
#' variance among the means of their sources (`between`), their phenotype (`p_imm`) and their
#' genetic variance (`vi`), in the state `st`
#' @noRd
mc_flow_sources <- function(st) {
  eval(quote({
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
      zi <- immigrant_mean(em_to, za, disp_matrix, za)
      between <- pmax(immigrant_mean(em_to, za^2, disp_matrix, za^2) - zi^2, 0)
      p_imm <- immigrant_mean(em_to, pa, disp_matrix, pa)
      if (isTRUE(ev$dyn)) {
        g0 <- ifelse(is.na(ga), 0, ga)
        vi <- immigrant_mean(em_to, g0, disp_matrix, g0)
      }
    }
    if (isTRUE(ev$sexual)) {
      d0 <- ifelse(is.na(da), 0, da)
      vi_d <- immigrant_mean(em_to, d0, disp_matrix, d0)
    }
  }), st)
  invisible(st)
}

#' the genetic means, variances and phenotypes of the populations that received immigrants,
#' and of those founded by them, with one trait, in the state `st`
#' @noRd
mc_flow_merge <- function(st) {
  eval(quote({
    pa_now <- pa
    both <- res > 0 & imm_m > 0
    if (isTRUE(ev$dyn)) {
      # variance of the merged population: the mean of the variances of residents and
      # immigrants; the difference between their means is not added, as if selection
      # removed it
      w <- res[both] / (res[both] + imm_m[both])
      ga[both] <- w * ga[both] + (1 - w) * vi[both]
      if (isTRUE(ev$sexual)) da[both] <- w * da[both] + (1 - w) * vi_d[both]
    }
    za[both] <- (res[both] * za[both] + imm_m[both] * zi[both]) / (res[both] + imm_m[both])
    pa_now[both] <- (res[both] * pa[both] + imm_m[both] * p_imm[both]) /
      (res[both] + imm_m[both])
    founded <- res == 0 & imm_m > 0
    if (any(founded) && lottery) {
      # with a lottery, founders are sampled once, when they recruit
      za[founded] <- zi[founded]
      if (isTRUE(ev$dyn)) ga[founded] <- vi[founded]
      if (isTRUE(ev$sexual)) da[founded] <- vi_d[founded]
      extra_v[founded] <- between[founded]
      no_drift[founded] <- !ev$founder
      pa_now[founded] <- p_imm[founded]
    } else if (any(founded)) {
      if (isTRUE(ev$dyn)) {
        # founders carry the variance of their sources, less 1 / (2k) for k founders
        za[founded] <- found(zi, imm_m, founded, alive, g_cells = vi + between)
        ga[founded] <- (vi * pmax(1 - 1 / (2 * imm_m), 0))[founded]
        if (isTRUE(ev$sexual)) da[founded] <- (vi_d * pmax(1 - 1 / (2 * imm_m), 0))[founded]
      } else {
        za[founded] <- found(zi, imm_m, founded, alive, between)
      }
      pa_now[founded] <- p_imm[founded]
    }
    if (isTRUE(ev$dyn)) vg[, alive] <- ga
    if (isTRUE(ev$sexual)) vd[, alive] <- da
    pa <- pa_now
  }), st)
  invisible(st)
}

#' recruitment by lottery: the microsites left by the survivors are filled by a draw among
#' the offspring there; changes `n`, `imm` and the evolutionary state in `st`
#' @noRd
mc_step_recruit <- function(st) {
  eval(quote({
    # the microsites left by the survivors are filled by a draw among the offspring there
    if (!is.null(inside)) n[!inside] <- 0L
    rec <- lottery_draw(n, as.integer(pmax(microsites - rowSums(surv), 0)))
    # immigrants among the recruits, for the colonisation history
    imm <- imm + matrix(0L, m, k)
    from_imm <- rec > 0 & imm > 0
    if (any(from_imm)) {
      imm[from_imm] <- stats::rhyper(sum(from_imm), pmin(imm, n)[from_imm],
                                     pmax(n - imm, 0)[from_imm], rec[from_imm])
    }
    imm[!from_imm] <- 0L
    n_new <- surv + rec
    if (isTRUE(ev$loci_on)) {
      loci_step_recruit(environment())
    } else if (!is.null(ev)) {
      mc_recruit_traits(environment())
    }
    n <- n_new
  }), st)
  invisible(st)
}

#' one update of the genetic variances `ga` (and, with variance_selection = "sexual", the
#' disequilibria `da`) in the cells `cells`, for populations of size `size` after a fraction
#' `f` of a generation; absent populations become NA; changes `ga`, `vg`, `da` and `vd` in
#' the state `st`
#' @noRd
mc_variance_step <- function(st, cells, size, width, f) {
  ev <- st$ev
  alive <- st$alive
  if (isTRUE(ev$sexual)) {
    up <- variance_update_sexual(st$ga, st$da, size, width, ev, alive, f)
    st$da[cells] <- up$d[cells]
    st$da[!cells] <- NA_real_
    st$vd[, alive] <- st$da
    st$ga[cells] <- up$g[cells]
  } else {
    st$ga[cells] <- variance_update(st$ga, size, width, ev, alive, f)[cells]
  }
  st$ga[!cells] <- NA_real_
  st$vg[, alive] <- st$ga
  invisible(st)
}

#' genetic means and variances after a lottery: those of the survivors and of the recruits,
#' weighted by their numbers, with the drift of the recruits; changes `za`, `pa` and `vg`
#' @noRd
mc_recruit_traits <- function(st) {
  eval(quote({
    g <- if (isTRUE(ev$dyn)) ga else matrix(ev$g[alive], m, k, byrow = TRUE)
    g[is.na(g)] <- 0
    v_rec <- g + extra_v
    drift <- rec > 0 & v_rec > 0 & !no_drift
    if (is.finite(ev$nu) && any(drift)) {
      # recruits are a sample of the offspring (and founders of their sources)
      za[drift] <- za[drift] +
        stats::rnorm(sum(drift), 0, sqrt(v_rec[drift] / (ev$nu * rec[drift])))
    }
    has <- n_new > 0
    zs <- ifelse(surv > 0, z_surv, 0)
    zr <- ifelse(rec > 0, za, 0)
    za[has] <- ((surv * zs + rec * zr) / n_new)[has]
    if (isTRUE(ev$dyn)) {
      gs <- ifelse(surv > 0, g_surv, 0)
      ga[has] <- ((surv * gs + rec * g) / n_new)[has]
      if (isTRUE(ev$sexual)) {
        ds <- ifelse(surv > 0, d_surv, 0)
        dr <- ifelse(is.na(da), 0, da)
        da[has] <- ((surv * ds + rec * dr) / n_new)[has]
      }
      # mutation, drift and selection on the variance, in proportion to the recruits
      width <- ga + matrix(ev$sigma[alive]^2 + ev$ve[alive], m, k, byrow = TRUE)
      mc_variance_step(environment(), has, n_new, width, rec / pmax(n_new, 1))
    }
    za[!has] <- NA_real_
    pa <- za + sweep(matrix(e_now - ev$ref, m, k), 2, ev$b[alive], "*")
  }), st)
  invisible(st)
}

#' losses of populations at the end of the step: extirpation, disturbance and the sites
#' outside the pool; changes `n` in the state `st`
#' @noRd
mc_step_losses <- function(st) {
  eval(quote({
    # no draws without extirpation, so that seeds give the same results as before
    if (any(p_ex > 0)) n[stats::runif(m * k) < p_ex] <- 0L
    p_dist <- if (is.matrix(disturbance)) disturbance[, step] else disturbance
    if (any(p_dist > 0)) {
      hit <- stats::runif(m) < p_dist
      if (all(disturbance_sensitivity == 1)) {
        n[hit, ] <- 0L
      } else if (any(hit)) {
        # each population of a disturbed site is lost with the sensitivity of its species
        sub <- n[hit, , drop = FALSE]
        sub[stats::runif(length(sub)) < rep(disturbance_sensitivity[alive], each = sum(hit))] <-
          0L
        n[hit, ] <- sub
      }
    }
    if (!is.null(inside)) n[!inside] <- 0L
  }), st)
  invisible(st)
}

#' writes the genetic means and phenotypes of the alive species back to `z` and `ph` (NA where
#' the population is gone), in the state `st`
#' @noRd
mc_step_store <- function(st) {
  eval(quote({
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
      if (any(fixed_z[alive])) {
        za[, fixed_z[alive]] <- pool_mean[, alive, drop = FALSE][, fixed_z[alive]]
      }
      za[n == 0] <- NA_real_
      pa[n == 0] <- NA_real_
      if (isTRUE(ev$dyn)) {
        va <- vg[, alive, drop = FALSE]
        va[n == 0] <- NA_real_
        vg[, alive] <- va
      }
      if (isTRUE(ev$sexual)) {
        dv <- vd[, alive, drop = FALSE]
        dv[n == 0] <- NA_real_
        vd[, alive] <- dv
      }
      if (isTRUE(ev$loci_on)) loci_step_store(environment())
      z[, alive] <- za
      ph[, alive] <- pa
    }
  }), st)
  invisible(st)
}

#' the colonisation history after this step: arrivals, founders, immigrants and migrant load
#' of the populations of the alive species, in the state `st`
#' @noRd
mc_step_colonisation <- function(st) {
  eval(quote({
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
      cols[new_pop] <- as.integer((if (lottery) n else start + imm)[new_pop])
      propagule[, alive] <- cols
      cols <- immigrants[, alive, drop = FALSE]
      cols[new_pop] <- 0
      immigrants[, alive] <- cols
      cols <- load[, alive, drop = FALSE]
      cols[new_pop] <- 0
      load[, alive] <- cols
    }
    if (any(stay) && (any(imm > 0) || any(added > 0))) mc_migrant_load(environment())
    lost <- was & !now
    if (any(lost)) {
      cols <- arrival[, alive, drop = FALSE]
      cols[lost] <- NA_integer_
      arrival[, alive] <- cols
    }
    present[, alive] <- now
  }), st)
  invisible(st)
}

#' immigrants and seeds into populations already there: their number, and the share of the
#' population they make up, summed over steps (the migrant load), in the state `st`
#' @noRd
mc_migrant_load <- function(st) {
  eval(quote({
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
  }), st)
  invisible(st)
}

#' records the abundances (and the genetic means) at the steps of `record_steps`
#' @noRd
mc_step_record <- function(st) {
  eval(quote({
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

#' the output of mc_dynamics() from the state `st`
#' @noRd
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
  }), st)
  if (isTRUE(st$ev$mv)) {
    mc_finish_evolution_mv(st)
  } else if (!is.null(st$ev)) {
    mc_finish_evolution(st)
  }
  st$out
}

#' adds to the output `out` of the state `st` the evolutionary state with several traits
#' @noRd
mc_finish_evolution_mv <- function(st) {
  eval(quote({
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
                          since_change = since_change(ev$env, n_steps, record_steps, sites),
                          reference = stats::setNames(ev$ref, ev$axes),
                          genetic_var = stats::setNames(ev$G, species),
                          env_var = stats::setNames(ev$E, species),
                          plasticity = matrix(ev$b, k_all, dimnames = list(species, ev$axes)),
                          settings = structure(evolution[setdiff(names(evolution),
                                                                 c("env", "traits"))],
                                               class = "mc_evolution"))
  }), st)
  invisible(st)
}

#' adds to the output `out` of the state `st` the evolutionary state with one trait
#' @noRd
mc_finish_evolution <- function(st) {
  eval(quote({
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
    if (!is.null(vd)) {
      vd[absent] <- NA_real_
      dimnames(vd) <- labels
    }
    out$evolution <- list(genetic = z, plastic = plastic, phenotype = ph, history = z_history,
                          environment = ev$env(n_steps), environment_history = env_history,
                          since_change = since_change(ev$env, n_steps, record_steps, sites),
                          reference = ev$ref, genetic_variance = vg,
                          linkage_disequilibrium = vd,
                          settings = structure(evolution[setdiff(names(evolution),
                                                                 c("env", "traits"))],
                                               class = "mc_evolution"))
    if (isTRUE(ev$loci_on)) {
      out$evolution$allele_frequency <- loci_output(lp, absent, sites, species)
    }
  }), st)
  invisible(st)
}

#' steps since the environment of each site last changed (on any axis), at the end and at the
#' recorded steps: 0 in a step where it changed, and the step itself where it never changed
#' @noRd
since_change <- function(env, n_steps, record_steps, sites) {
  last <- rep(0L, length(sites))
  before <- as.matrix(env(1))
  steps <- sort(unique(as.integer(record_steps[record_steps <= n_steps])))
  at <- matrix(NA_integer_, length(sites), length(steps))
  if (length(steps)) dimnames(at) <- list(sites, step_label(steps))
  for (s in seq_len(n_steps)) {
    now <- as.matrix(env(s))
    if (s > 1) last[rowSums(now != before) > 0] <- s
    before <- now
    if (s %in% steps) at[, step_label(s)] <- s - last
  }
  list(end = stats::setNames(as.integer(n_steps - last), sites), history = if (length(steps)) at)
}

#' recruits of a lottery: in each site, `vacant` microsites are filled by a draw without
#' replacement among the offspring `pool` (sites x species), species by species
#' @noRd
lottery_draw <- function(pool, vacant) {
  out <- pool
  full <- rowSums(pool) > vacant
  if (any(full)) {
    sub <- pool[full, , drop = FALSE]
    left <- as.integer(vacant[full])
    rest <- rowSums(sub)
    res <- matrix(0L, nrow(sub), ncol(sub))
    for (j in seq_len(ncol(sub))) {
      rest <- rest - sub[, j]
      res[, j] <- if (j == ncol(sub)) left else
        stats::rhyper(nrow(sub), sub[, j], rest, left)
      left <- left - res[, j]
    }
    out[full, ] <- res
  }
  out
}

#' immigrants when the emigrants of each species at each site go to the other sites
#' emigration probabilities at each step (steps x all species) from `dispersal_time`: a vector of
#' multipliers of `dispersal`, one per step, or a matrix of probabilities (steps x species);
#' NULL when the probabilities are constant
#' @noRd
dispersal_by_step <- function(dispersal_time, dispersal, n_steps, k) {
  if (is.null(dispersal_time)) return(NULL)
  if (is.matrix(dispersal_time)) return(dispersal_time_matrix(dispersal_time, n_steps, k))
  bad <- !is.numeric(dispersal_time) || length(dispersal_time) < n_steps
  bad <- bad || anyNA(dispersal_time) || any(dispersal_time < 0)
  if (bad || !all(is.finite(dispersal_time))) {
    stop("`dispersal_time` needs one multiplier of `dispersal` (0 or more) per step (",
         n_steps, "), or a matrix of emigration probabilities", call. = FALSE)
  }
  pmin(outer(dispersal_time[seq_len(n_steps)], dispersal), 1)
}

#' emigration probabilities at each step (steps x all species) from a matrix `dispersal_time`
#' with one column, or one per species; stops unless it holds probabilities for every step
#' @noRd
dispersal_time_matrix <- function(dispersal_time, n_steps, k) {
  bad <- nrow(dispersal_time) < n_steps || !ncol(dispersal_time) %in% c(1, k)
  bad <- bad || !is.numeric(dispersal_time) || anyNA(dispersal_time)
  if (bad || any(dispersal_time < 0) || any(dispersal_time > 1)) {
    stop("a matrix `dispersal_time` needs emigration probabilities (0 to 1) with one row ",
         "per step (", n_steps, ") and one column, or one per species (", k, ")",
         call. = FALSE)
  }
  matrix(dispersal_time[seq_len(n_steps), ], n_steps, k)
}

#' dispersal matrices to switch among, checked and made row-stochastic with a zero diagonal;
#' NULL when the matrix of the landscape holds at every step
#' @noRd
kernels_by_step <- function(kernels, kernel_step, n_steps, m, disp_matrix,
                            kernel_group = NULL) {
  if (is.null(kernels) && is.null(kernel_step)) return(NULL)
  if (is.null(kernels) || is.null(kernel_step)) {
    stop("`kernels` and `kernel_step` go together", call. = FALSE)
  }
  if (is.null(disp_matrix)) {
    stop("`kernels` needs a spatially explicit landscape", call. = FALSE)
  }
  if (is.matrix(kernels)) kernels <- list(kernels)
  if (!is.list(kernels) || length(kernels) == 0) {
    stop("`kernels` must be a list of dispersal matrices", call. = FALSE)
  }
  kernels <- lapply(kernels, kernel_matrix, m = m)
  check_kernel_step(kernel_step, kernels, kernel_group, n_steps)
  kernels
}

#' one of `kernels` (a matrix, or a landscape with one), checked and made row-stochastic with
#' a zero diagonal and no names; returns the m x m matrix
#' @noRd
kernel_matrix <- function(x, m) {
  if (inherits(x, "mc_landscape")) x <- x$dispersal
  bad <- !is.matrix(x) || !is.numeric(x) || any(dim(x) != m)
  if (bad || anyNA(x) || any(x < 0)) {
    stop("each of `kernels` must be a non-negative ", m, " x ", m, " matrix",
         call. = FALSE)
  }
  x <- unname(x)
  diag(x) <- 0
  if (any(rowSums(x) <= 0)) {
    stop("each site needs at least one destination in every matrix of `kernels`",
         call. = FALSE)
  }
  # rows already summing to 1 are kept as they are, so that the matrix of the landscape
  # given as a kernel reproduces a run without `kernels` to the last bit
  total <- rowSums(x)
  off <- abs(total - 1) > 1e-12
  x[off, ] <- x[off, ] / total[off]
  x
}

#' stops unless `kernel_step` gives, for each step (and each group of species, with
#' `kernel_group`), the index or the name of one of `kernels`
#' @noRd
check_kernel_step <- function(kernel_step, kernels, kernel_group, n_steps) {
  if (is.character(kernel_step)) {
    if (is.null(names(kernels)) || !all(kernel_step %in% names(kernels))) {
      stop("a character `kernel_step` must hold names of `kernels`", call. = FALSE)
    }
  } else {
    bad <- !is.numeric(kernel_step) || anyNA(kernel_step)
    bad <- bad || any(kernel_step != round(kernel_step)) || any(kernel_step < 1)
    if (bad || any(kernel_step > length(kernels))) {
      stop("`kernel_step` must give, for each step, the index or name of one of `kernels`",
           call. = FALSE)
    }
  }
  if (is.matrix(kernel_step) != !is.null(kernel_group)) {
    stop("a matrix `kernel_step` (steps x groups) goes with `kernel_group`, and a vector ",
         "without it", call. = FALSE)
  }
  if (NROW(kernel_step) < n_steps) {
    stop("`kernel_step` needs one value per step (one row per step with `kernel_group`; ",
         n_steps, " steps)", call. = FALSE)
  }
  invisible(kernel_step)
}

#' the column of a matrix `kernel_step` used by each species, from `kernel_group` (names of the
#' columns or their indices); NULL without groups
#' @noRd
kernel_groups <- function(kernel_group, kernel_step, k) {
  if (is.null(kernel_group)) return(NULL)
  if (!is.matrix(kernel_step)) {
    stop("`kernel_group` needs `kernels` and a matrix `kernel_step` (steps x groups)",
         call. = FALSE)
  }
  if (length(kernel_group) != k || anyNA(kernel_group)) {
    stop("`kernel_group` needs one group per species (", k, ")", call. = FALSE)
  }
  if (is.factor(kernel_group)) kernel_group <- as.character(kernel_group)
  kernel_group_columns(kernel_group, kernel_step)
}

#' the columns of `kernel_step` that `kernel_group` names (a character vector) or indexes
#' (whole numbers); returns their indices as integers
#' @noRd
kernel_group_columns <- function(kernel_group, kernel_step) {
  if (is.character(kernel_group)) {
    if (is.null(colnames(kernel_step)) || !all(kernel_group %in% colnames(kernel_step))) {
      stop("a character `kernel_group` must hold column names of `kernel_step`", call. = FALSE)
    }
    return(match(kernel_group, colnames(kernel_step)))
  }
  bad <- !is.numeric(kernel_group) || any(kernel_group != round(kernel_group))
  if (bad || any(kernel_group < 1) || any(kernel_group > ncol(kernel_step))) {
    stop("`kernel_group` must give, for each species, a column of `kernel_step`", call. = FALSE)
  }
  as.integer(kernel_group)
}

#' expected immigrants (sites x species) from emigrants `x` (sites x species): t(D) %*% x, or,
#' with one kernel per group of species (a list with attribute `group`), each species through
#' the kernel of its group
#' @noRd
kernel_flow <- function(disp_matrix, x) {
  if (!is.list(disp_matrix)) return(crossprod(disp_matrix, x))
  group <- attr(disp_matrix, "group")
  out <- matrix(0, ncol(disp_matrix[[1]]), ncol(x))
  for (g in unique(group)) {
    cols <- which(group == g)
    out[, cols] <- crossprod(disp_matrix[[g]], x[, cols, drop = FALSE])
  }
  out
}

#' the dispersal matrix of species `s`: `disp_matrix` itself or, with one kernel per group of
#' species, that of the group of `s` (`group`: the group of each species)
#' @noRd
species_kernel <- function(disp_matrix, s, group = attr(disp_matrix, "group")) {
  if (is.list(disp_matrix)) disp_matrix[[group[s]]] else disp_matrix
}

#' D %*% x, by group as in kernel_flow()
#' @noRd
kernel_product <- function(disp_matrix, x) {
  if (!is.list(disp_matrix)) return(disp_matrix %*% x)
  group <- attr(disp_matrix, "group")
  out <- matrix(0, nrow(disp_matrix[[1]]), ncol(x))
  for (g in unique(group)) {
    cols <- which(group == g)
    out[, cols] <- disp_matrix[[g]] %*% x[, cols, drop = FALSE]
  }
  out
}

#' multinomially, with the probabilities of the rows of `disp_matrix`
#' @noRd
place_emigrants <- function(em, disp_matrix, open_to = NULL) {
  if (is.list(disp_matrix)) {
    group <- attr(disp_matrix, "group")
    out <- matrix(0L, nrow(em), ncol(em))
    for (g in unique(group)) {
      cols <- which(group == g)
      out[, cols] <- place_emigrants(em[, cols, drop = FALSE], disp_matrix[[g]],
                                     if (!is.null(open_to)) open_to[, cols, drop = FALSE])
    }
    return(out)
  }
  m <- nrow(em)
  immigrants <- matrix(0L, m, ncol(em))
  for (cell in which(em > 0)) {
    from <- (cell - 1L) %% m + 1L
    species <- (cell - 1L) %/% m + 1L
    immigrants[, species] <- immigrants[, species] +
      drop(stats::rmultinom(1, em[cell], if (is.null(open_to)) disp_matrix[from, ] else
        disp_matrix[from, ] * open_to[, species]))
  }
  immigrants
}

#' a sites x species matrix in the order of `sites` and `species`: rows by name when they have
#' names, by position otherwise; columns by species name
#' @noRd
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
