#' Speciation settings
#'
#' Settings for speciation in [simulate_evolving_metacommunity()]. Speciation is checked once
#' per generation (`steps_per_generation` of [evolution_control()]). A species buds when its
#' populations fall into two or more groups: the largest group (most sites, then most
#' individuals) keeps the identity and age of the species, and each other group becomes a new
#' species with its populations, abundances and genetic means. Several new species at once
#' branch off one after another, the most divergent first, at the same time (internal branches
#' of length zero), so that the tree stays binary.
#'
#' Four modes start new species (the initiation of Hubert et al. 2015):
#'
#' * `"geographic"`: each pair of populations of a species has a divergence \eqn{d}{d}, in
#'   generations. With `divergence = "cluster"` (the default, close to gen3sis; Hagen et al.
#'   2021), populations linked by a migration rate of at least `migration_min` form clusters
#'   (chains of links). \eqn{d}{d} grows by 1 per generation between populations of different
#'   clusters; within a cluster it stays as it is, and returns to 0 only after
#'   `reset_generations` generations in a row in the same cluster, so that a single generation
#'   of contact does not erase a long isolation. With `divergence = "gradual"`, \eqn{d}{d} grows
#'   by 1 per generation and gene flow lowers it in proportion to the migrants exchanged: \eqn{d
#'   \leftarrow \max(0, d + 1 - h M)}{d <- max(0, d + 1 - h M)}, with \eqn{M}{M} the expected
#'   migrants between the two populations in the generation and \eqn{h}{h} = `homogenization`.
#'   Populations are grouped by single linkage: two are in the same group when a chain of pairs
#'   with \eqn{d <}{d <} `threshold` links them.
#' * `"ecological"`: each pair of populations counts the generations in a row in which their
#'   genetic means (`evolution`, with genetic variance) differ by at least `threshold` and their
#'   migration rate is below `migration_max`. Populations are grouped by single linkage on pairs
#'   whose count is below `min_generations`. Requiring the difference to last, with reduced gene
#'   flow, keeps drift and founder sampling in small populations from splitting species in a
#'   uniform environment.
#' * `"point"`: point mutation, the mode of neutral theory and of MESS (Overcast et al. 2021).
#'   Each individual becomes, with probability `rate` per generation, the first individual of a
#'   new species in its site. The mother keeps its populations, except where a mutant takes its
#'   last individual. Mutants arrive in the generation of their birth.
#' * `"fission"`: random fission. Each species with two or more populations splits with
#'   probability \eqn{1 - e^{-\lambda}}{1 - exp(-lambda)} per generation (\eqn{\lambda}{lambda}
#'   = `rate`). With `partition = "spatial"`, the populations of the new group are those nearest
#'   (by the dispersal kernel) to a random occupied site, in a random number from 1 to all but
#'   one; with `"random"` (or on a spatially implicit landscape), each population joins it with
#'   probability 1/2, with at least one population on each side.
#'
#' Gene flow between two populations \eqn{a}{a} and \eqn{b}{b} is measured by the migration
#' rate, the larger of the expected migrants from one to the other in the generation over the
#' size of the receiving population, so that it does not depend on abundance. A new population
#' takes the divergences (or counts) of the population that sent it most migrants, and a
#' population lost at any step of a generation founds anew.
#'
#' New species take the genetic means of their populations; with `jump_var > 0` their genetic
#' means (and the optimum of their regional pool) shift by one Gaussian draw with that variance
#' per new species, as in the punctuated change of MESS (Overcast et al. 2021).
#'
#' **Protracted speciation.** With `completion > 0`, each new lineage is incipient until it
#' completes, after `completion` generations (`completion_dist = "fixed"`) or after an
#' exponential time with that mean (`"exponential"`, the rate \eqn{\lambda_2}{lambda2} of the
#' protracted birth-death model; Etienne & Rosindell 2012). Completion only changes what is
#' recognised as a species: incipient lineages are separate lineages in the dynamics. At the
#' end, a lineage that has completed is a good species, and an incipient lineage belongs to the
#' species its parent belonged to at its birth, as in the labels of the PBD. The reconstructed
#' tree has one tip per good species with extant lineages, represented by the good lineage
#' itself or, if it is extinct, by its earliest extant member. Its nodes are times of
#' initiation of lineages, which for a species whose lineages are not monophyletic may predate
#' the initiation of the species itself.
#'
#' @param mode `"geographic"`, `"ecological"`, `"point"` or `"fission"`.
#' @param threshold For `"geographic"` and `"ecological"`: divergence in generations, or
#'   difference of genetic means in units of the environment, that separates groups: one
#'   number greater than 0; `Inf` turns speciation off. With several traits (one per niche axis),
#'   the difference of two
#'   populations in `"ecological"` is the Mahalanobis distance of their genetic means under
#'   the phenotypic matrix of the species, \eqn{\sqrt{\Delta\bar z^\top P^{-1} \Delta\bar
#'   z}}{sqrt(dz' P^-1 dz)}, so `threshold` is in units of phenotypic standard deviations.
#' @param rate For `"point"`: probability per individual per generation, one number between 0
#'   and 1; for `"fission"`: rate per species per generation, one finite number, 0 or more.
#' @param divergence For `"geographic"`: `"cluster"` or `"gradual"`.
#' @param migration_min For `divergence = "cluster"`: migration rate per generation that links
#'   two populations, one number greater than 0 (`Inf` links none). Where clusters form
#'   depends on dispersal, the kernel and abundances: the rate of the strongest link of each
#'   population should straddle it. In the design of
#'   `data-raw/evolving_migration_min.R` (strongest links with median 0.006), species did not
#'   split at all at 0.002 or below and split almost as with no links at 0.02 or above.
#' @param reset_generations For `divergence = "cluster"`: generations in a row in the same
#'   cluster after which divergence returns to 0; one whole number, 1 or more.
#' @param homogenization For `divergence = "gradual"`: reduction of divergence per migrant;
#'   one number, 0 or more.
#' @param min_generations For `"ecological"`: generations in a row of difference and reduced
#'   gene flow before groups split; one whole number, 1 or more.
#' @param migration_max For `"ecological"`: migration rate per generation below which gene flow
#'   counts as reduced; one number, 0 or more (`Inf` for no limit on gene flow).
#' @param partition For `"fission"`: `"spatial"` or `"random"`.
#' @param jump_var Variance of the shift of the genetic mean of each new species: one number, 0
#'   or more. With several traits, one value (independent shifts of that variance on every trait) or
#'   a covariance
#'   matrix with one row and column per trait (for instance proportional to \eqn{G}{G}).
#' @param completion Time to complete speciation, in generations (fixed or mean): one number,
#'   0 or more; 0 for instantaneous speciation.
#' @param completion_dist `"fixed"` or `"exponential"`.
#' @return An `mc_speciation` object.
#' @references Etienne, R. S., & Rosindell, J. (2012). Prolonging the past counteracts the pull
#'   of the present: protracted speciation can explain observed slowdowns in diversification.
#'   *Systematic Biology*, 61, 204--213. \doi{10.1093/sysbio/syr091}
#'
#'   Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M., Rangel, T. F., &
#'   Pellissier, L. (2021). gen3sis: A general engine for eco-evolutionary simulations of the
#'   processes that shape Earth's biodiversity. *PLOS Biology*, 19, e3001340.
#'   \doi{10.1371/journal.pbio.3001340}
#'
#'   Hubert, N., Calcagno, V., Etienne, R. S., & Mouquet, N. (2015). Metacommunity speciation
#'   models and their implications for diversification theory. *Ecology Letters*, 18,
#'   864--881. \doi{10.1111/ele.12458}
#'
#'   Overcast, I., Ruffley, M., Rosindell, J., Harmon, L., Borges, P. A. V., Emerson, B. C.,
#'   Etienne, R. S., Gillespie, R., Krehenwinkel, H., Mahler, D. L., Massol, F., Parent, C. E.,
#'   Patiño, J., Peter, B., Week, B., Wagner, C., Hickerson, M. J., & Rominger, A. (2021). A
#'   unified model of species abundance, genetic diversity, and functional diversity reveals the
#'   mechanisms structuring ecological communities. *Molecular Ecology Resources*, 21,
#'   2782--2800. \doi{10.1111/1755-0998.13514}
#' @family simulate metacommunities
#' @examples
#' speciation_control("geographic", threshold = 50)
#' speciation_control("geographic", threshold = 50, divergence = "gradual", homogenization = 0.5)
#' speciation_control("ecological", threshold = 20, min_generations = 30)
#' speciation_control("point", rate = 1e-4, completion = 50)
#' speciation_control("fission", rate = 0.01, partition = "random")
#' @export
speciation_control <- function(mode = c("geographic", "ecological", "point", "fission"),
                               threshold, rate, divergence = c("cluster", "gradual"),
                               migration_min = 0.005, reset_generations = 10,
                               homogenization = 1, min_generations = 20, migration_max = 0.001,
                               partition = c("spatial", "random"), jump_var = 0,
                               completion = 0, completion_dist = c("fixed", "exponential")) {
  mode <- match.arg(mode)
  divergence <- match.arg(divergence)
  partition <- match.arg(partition)
  completion_dist <- match.arg(completion_dist)
  if (mode %in% c("geographic", "ecological")) {
    if (missing(threshold)) stop("`threshold` is needed", call. = FALSE)
    # Inf turns speciation off
    check_number(threshold, "threshold", lower = 0, open = c(TRUE, FALSE),
                 finite = FALSE)
    rate <- NA_real_
  } else {
    if (missing(rate)) stop("`rate` is needed", call. = FALSE)
    # a probability per individual for "point", a rate per species for "fission"
    check_number(rate, "rate", 0, if (mode == "point") 1 else Inf)
    threshold <- NA_real_
  }
  # infinite migration rates are valid limits: no links (`migration_min`) or no guard on gene
  # flow (`migration_max`)
  check_number(migration_min, "migration_min", lower = 0, open = c(TRUE, FALSE),
               finite = FALSE)
  check_count(reset_generations, "reset_generations")
  check_number(homogenization, "homogenization", lower = 0)
  check_count(min_generations, "min_generations")
  check_number(migration_max, "migration_max", lower = 0, finite = FALSE)
  if (is.matrix(jump_var)) {
    if (!is.numeric(jump_var) || !all(is.finite(jump_var))) {
      stop("a matrix `jump_var` must hold finite numbers", call. = FALSE)
    }
  } else {
    check_number(jump_var, "jump_var", lower = 0)
  }
  check_number(completion, "completion", lower = 0)
  structure(list(mode = mode, threshold = threshold, rate = rate, divergence = divergence,
                 migration_min = migration_min, reset_generations = reset_generations,
                 homogenization = homogenization, min_generations = min_generations,
                 migration_max = migration_max, partition = partition, jump_var = jump_var,
                 completion = completion, completion_dist = completion_dist),
            class = "mc_speciation")
}

#' @export
print.mc_speciation <- function(x, ...) {
  cat("<mc_speciation>", x$mode)
  if (x$mode == "geographic") {
    cat(" | threshold =", format(x$threshold), "| divergence:", x$divergence,
        if (x$divergence == "cluster") {
          paste("| migration_min =", format(x$migration_min), "| reset_generations =",
                format(x$reset_generations))
        } else {
          paste("| homogenization =", format(x$homogenization))
        })
  } else if (x$mode == "ecological") {
    cat(" | threshold =", format(x$threshold), "| min_generations =",
        format(x$min_generations), "| migration_max =", format(x$migration_max))
  } else {
    cat(" | rate =", format(x$rate))
    if (x$mode == "fission") cat(" | partition:", x$partition)
  }
  if (is.matrix(x$jump_var)) {
    cat(" | jump_var:", nrow(x$jump_var), "x", ncol(x$jump_var), "matrix")
  } else if (x$jump_var > 0) {
    cat(" | jump_var =", format(x$jump_var))
  }
  if (x$completion > 0) cat(" | completion:", x$completion_dist, format(x$completion))
  cat("\n")
  invisible(x)
}

#' Simulate a metacommunity with speciation
#'
#' Runs the model of [simulate_metacommunity()], with the evolution of the niche optimum within
#' populations of [evolution_control()], on species that bud as their populations diverge, by
#' point mutation or by fission ([speciation_control()]), and die out when their last
#' population is lost. The tree is not
#' given: it grows out of the dynamics. The community starts from the species of `traits`
#' (usually one ancestor) and is closed: there is no immigration from outside the landscape.
#'
#' Competition follows the tree as it grows and, with `niche_overlap = TRUE`, the phenotypes in
#' each site: \eqn{\alpha_{ii}}{alpha_ii} = `scale` and, between species in site \eqn{s}{s},
#' \deqn{\alpha_{ij}(s) = c \, (q_0 + (q_1 - q_0) e^{-d_{ij} / \delta}) \,
#' \exp\left(-\frac{(p_{is} - p_{js})^2}{2(\sigma_i^2 + \sigma_j^2)}\right),}{alpha_ij(s) = c
#' (q0 + (q1 - q0) exp(-d_ij / delta)) exp(-(p_is - p_js)^2 / (2 (sigma_i^2 + sigma_j^2))),}
#' with \eqn{c}{c} = `scale`, \eqn{q_1}{q1} = `inter_max`, \eqn{q_0}{q0} = `inter_min`,
#' \eqn{d_{ij}}{d_ij} the patristic distance in generations, \eqn{\delta}{delta} =
#' `distance_scale`, \eqn{p}{p} the phenotypes (genetic mean and plastic part) and
#' \eqn{\sigma}{sigma} the niche breadths. The last factor is the overlap of two Gaussian
#' niches, 1 for equal phenotypes, so that species that have diverged in niche compete less.
#' With `inter_max = 1`, sisters just split with the same phenotype compete as members of one
#' species. With `inter_min = 0` and no niche divergence, distant species do not compete at all,
#' and richness is limited only by speciation and extinction. With `competition = "none"`,
#' species do not compete with each other.
#'
#' New species inherit the traits of their mother (`sigma`, `r_max`, `dispersal`,
#' `extirpation`, the settings of evolution) and, in each of their sites, its abundance and
#' genetic mean. Several species at the start need a tree in `traits` (from [species_traits()])
#' with branch lengths in generations, so that their distances are known.
#'
#' Without a split, the simulation is the same as [simulate_metacommunity()] with the same
#' `evolution` and the same seed (with one species, or several with `competition = "none"`).
#'
#' @inheritParams simulate_metacommunity
#' @param traits Traits of the species at the start, from [species_traits()] (one row per
#'   species), or a list of such tables, one per environmental axis of `env`, for a niche of
#'   several axes with one evolving trait per axis (the tree, when needed, goes with the first
#'   table).
#' @param speciation Settings from [speciation_control()].
#' @param evolution Settings from [evolution_control()]; `NULL` keeps the optimum of each
#'   species fixed (`genetic_var = 0`), which rules out ecological speciation.
#' @param competition `"phylo"` or `"none"`.
#' @param inter_max,inter_min,scale,distance_scale Competition between species, relative to
#'   competition within a species: between sisters just split (`inter_max`) and between
#'   distant species (`inter_min`); the overall scale (competition within a species); and the
#'   distance in generations over which the excess over `inter_min` decays by a factor of
#'   \eqn{e}{e}. Each one number: `inter_max` 0 or more, `inter_min` between 0 and `inter_max`,
#'   `scale` greater than 0 and `distance_scale` greater than 0 (`Inf` for no decay).
#' @param niche_overlap Should competition between species also fall with the difference of
#'   their phenotypes in each site? `TRUE` or `FALSE`.
#' @param initial Optional matrix (sites x species) of abundances at the start; by default
#'   Poisson(`seed_mean`) individuals of each species in every site.
#' @param max_species Maximum number of species ever, including extinct ones (one whole number,
#'   1 or more); the simulation stops with an error beyond it.
#' @param env For a list of `traits`: the environmental axes of each site, a matrix or data
#'   frame with one row per site and one column per element of `traits`. The `genetic_var`
#'   and `env_var` of `evolution` are then covariance matrices (see the section Several traits
#'   of [evolution_control()]); competition falls with the overlap of the niches on all axes,
#'   and new species inherit the matrices of their mother. A list of `traits` from
#'   [niche_correlation()] gives the niche a correlation between its axes, in growth,
#'   selection and the overlap of niches; with one matrix per species, new species take the
#'   correlation of their mother, and `traits` in the output keeps the correlation of the
#'   extant species.
#' @param kernel_group Optional routes of dispersal by group of species, as in
#'   [simulate_metacommunity()]: one group per species at the start, given as a column name or
#'   index of a matrix `kernel_step` (steps x groups). New species stay in the group of their
#'   mother. The gene flow that keeps the populations of a species together follows the routes
#'   of its group. `NULL` (the default) for one route for all species.
#' @return An `mc_evolving` object, a list with:
#'   * `abundance`: sites x extant species. With protracted speciation (`completion` of
#'     [speciation_control()]), the species are the good species, each with the individuals
#'     of its extant lineages; otherwise every lineage is a species.
#'   * `tree`: the reconstructed tree of those species, with branch lengths in generations,
#'     nodes at times of initiation and a root edge from the start to the crown, or `NULL`
#'     with fewer than two.
#'   * `tree_full`: the tree of all lineages, extinct ones included.
#'   * `species`: a table with the parent, birth (initiation), end, fate, order of budding,
#'     time of completion and status (`"good"` or `"incipient"` at its end) of every lineage.
#'   * `events`: one row per budding, with the new lineages.
#'   * `lineages` and `lineage_species`: abundances of the extant lineages and the species
#'     each belongs to.
#'   * `richness` and `richness_good`: extant lineages, and extant lineages that have
#'     completed, at each generation (the second can be below the number of good species, some
#'     of which survive only in incipient lineages).
#'   * `genetic`: genetic means of the populations of the extant lineages (sites x lineages,
#'     or sites x lineages x traits with several traits).
#'   * `traits`: of the extant species, with the mean genetic optimum over the populations of
#'     their lineages as `optimum`; with several traits, a list with one such table per axis,
#'     and `env`, the environmental axes.
#'   * `alpha`: among the extant species at the end, without the niche overlap, which varies
#'     among sites.
#'   * The records and colonisation history of [simulate_metacommunity()] for the extant
#'     species: records summed over their extant lineages (lineages extinct by the end are
#'     left out); colonisation, in each site, of the most abundant lineage.
#'   * `series`: one row per lineage alive in each generation, with its `occupancy` (sites),
#'     regional `abundance`, `genetic` mean (weighted by abundance; `genetic_<trait>` for each
#'     trait with several traits), and the sites it
#'     `colonised` and lost (`extirpated`) since the previous generation, out of those where it
#'     was absent (`absent_before`) and present (`present_before`). Sites that a lineage passed
#'     to a new one are not counted as lost, and the transitions are `NA` in the generation of
#'     birth. [persistence()] turns them into rates of colonisation and extirpation.
#'   * `ranges`: the changes of presence of every lineage in every site, one row per change
#'     (`generation`, `lineage`, `site` and `present`), with the starting sites at generation
#'     0 and the sites a new lineage took from its mother at its birth. [range_history()] maps
#'     them on the tree as ranges of areas.
#'   * `evolution`: with `evolution`, the state of the extant species at the end, for
#'     [sample_individuals()] with `source = "evolution"`: genetic means in each site
#'     (`genetic`, weighted over their lineages by abundance; `NA` where absent), the plastic
#'     part of the phenotype (`plastic`), the variance among the lineages of a species in
#'     each site (`within_var`, 0 for species of one lineage), the environment of each site
#'     at the end (`environment`), the same genetic means, variance among lineages and
#'     environment at `record_steps` for the lineages extant at the end (`history`,
#'     `within_history`, `environment_history`), the steps since the environment of each site
#'     last changed (`since_change`, as in [simulate_metacommunity()]; the environment of
#'     this engine is fixed, so they equal the step), the reference environment of plasticity
#'     (`reference`), the genetic and environmental variances and the slope of plasticity of
#'     each species (`genetic_var`, `env_var`, `plasticity`) and the `settings`; otherwise
#'     `NULL`. With several traits, the genetic means, plastic parts and variances among
#'     lineages (trait by trait) are arrays of sites x species x traits (with the recorded
#'     steps as a fourth dimension in the histories), `genetic_var` and `env_var` hold one
#'     matrix per species and `plasticity` is a matrix of species x traits.
#'   * The `landscape`, `dispersal`, `extirpation`, `speciation` and `n_steps` used, and
#'     `dispersal_time`, `kernels` (rescaled), `kernel_step`, `dispersal_site`,
#'     `disturbance_sensitivity` and `kernel_group` (one per species; `NULL` without groups).
#'     New species take the emigration through time, the sensitivity to disturbance and the
#'     group of routes of their mother. The gene flow that keeps populations together follows
#'     the kernel of each step (that of the group of the species, with `kernel_group`);
#'     fission by distance uses the kernel of `landscape`.
#'
#'   [as_metacommunity()] turns it into an `mc_sim` for the other functions of the package.
#' @references Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M., Rangel,
#'   T. F., & Pellissier, L. (2021). gen3sis: A general engine for eco-evolutionary simulations
#'   of the processes that shape Earth's biodiversity. *PLOS Biology*, 19, e3001340.
#'   \doi{10.1371/journal.pbio.3001340}
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' ls <- landscape_grid(6, disp_scale = 0.5)
#' ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
#' sim <- simulate_evolving_metacommunity(
#'   ls, ancestor, dispersal = 0.002, n_steps = 400,
#'   speciation = speciation_control("geographic", threshold = 100),
#'   evolution = evolution_control(genetic_var = 20, ne_ratio = 0.5))
#' sim
#' @export
simulate_evolving_metacommunity <- function(landscape, traits, speciation, evolution = NULL,
                                            dispersal = 0.01, extirpation = 0,
                                            n_steps = 1000, competition = c("phylo", "none"),
                                            inter_max = 1, inter_min = 0.25, scale = 0.05,
                                            distance_scale = 100, niche_overlap = TRUE,
                                            initial = NULL, seed_mean = 0.5,
                                            record_steps = NULL,
                                            immigration = c("poisson", "multinomial"),
                                            max_species = 500, env = NULL, disturbance = 0,
                                            dispersal_time = NULL, kernels = NULL,
                                            kernel_step = NULL, dispersal_site = NULL,
                                            disturbance_sensitivity = 1,
                                            kernel_group = NULL) {
  evolving_check_inputs(landscape, speciation)
  competition <- match.arg(competition)
  if (!is.null(evolution$mutation_var)) {
    stop("`mutation_var` is not yet available in `simulate_evolving_metacommunity()`",
         call. = FALSE)
  }
  if (!is.null(evolution$loci)) {
    stop("`loci` are not yet available in `simulate_evolving_metacommunity()`", call. = FALSE)
  }
  immigration <- match.arg(immigration)
  evolving_check_numbers(n_steps, inter_max, inter_min, scale, distance_scale, niche_overlap,
                         max_species)
  several <- is.list(traits) && !is.data.frame(traits)
  niche <- evolving_niche(traits, env, landscape, several)
  has_evolution <- !is.null(evolution)
  if (!has_evolution) evolution <- evolution_control(0, ne_ratio = Inf)
  evolving_check_evolution(evolution, several, speciation)
  comp <- list(competition = competition, inter_max = inter_max, inter_min = inter_min,
               scale = scale, distance_scale = distance_scale, niche_overlap = niche_overlap)
  # the arguments passed on to the dynamics, evaluated when the state is prepared
  run <- evolving_run(
    landscape, niche, speciation, evolution, comp, n_steps, max_species,
    list(dispersal = dispersal, seed_mean = seed_mean, record_steps = record_steps,
         extirpation = extirpation, initial = initial, immigration = immigration,
         disturbance = disturbance, dispersal_time = dispersal_time, kernels = kernels,
         kernel_step = kernel_step, dispersal_site = dispersal_site,
         disturbance_sensitivity = disturbance_sensitivity, kernel_group = kernel_group))
  evolving_output(run, landscape, speciation, comp, has_evolution)
}

#' stops unless `landscape` and `speciation` are objects of their classes
#' @noRd
evolving_check_inputs <- function(landscape, speciation) {
  if (!inherits(landscape, "mc_landscape")) {
    stop("`landscape` must be an `mc_landscape`, such as one from `landscape_grid()`",
         call. = FALSE)
  }
  if (!inherits(speciation, "mc_speciation")) {
    stop("`speciation` must come from `speciation_control()`", call. = FALSE)
  }
  invisible(landscape)
}

#' stops unless the numbers and the flag of simulate_evolving_metacommunity() are valid
#' @noRd
evolving_check_numbers <- function(n_steps, inter_max, inter_min, scale, distance_scale,
                                   niche_overlap, max_species) {
  check_count(n_steps, "n_steps")
  check_number(inter_max, "inter_max", lower = 0)
  check_number(inter_min, "inter_min", 0, inter_max)
  check_number(scale, "scale", lower = 0, open = TRUE)
  # an infinite scale keeps competition between species at `inter_max`
  check_number(distance_scale, "distance_scale", lower = 0, open = c(TRUE, FALSE),
               finite = FALSE)
  check_flag(niche_overlap, "niche_overlap")
  check_count(max_species, "max_species")
  invisible(n_steps)
}

#' the niche of the species at the start: the `traits` of the first axis and, with several
#' niche axes (one trait per axis, evolving with covariance matrices), the list of tables of
#' traits (`trait_list`), the environment of the sites on each axis (`env`) and the correlation
#' between the axes (`niche_cor`)
#' @noRd
evolving_niche <- function(traits, env, landscape, several) {
  if (!several) {
    if (!is.null(env)) stop("`env` needs a list of traits, one per column", call. = FALSE)
    return(list(traits = as_traits(traits), trait_list = NULL, env = env, several = FALSE))
  }
  niche_cor <- attr(traits, "niche_cor")
  trait_list <- lapply(traits, as_traits)
  env <- evolving_niche_env(env, traits, trait_list, landscape)
  sp1 <- rownames(trait_list[[1]])
  if (!all(vapply(trait_list, function(t) setequal(rownames(t), sp1), logical(1)))) {
    stop("simulate_evolving_metacommunity(): sets of traits have different species",
         call. = FALSE)
  }
  tree_attr <- attr(trait_list[[1]], "tree")
  trait_list <- lapply(trait_list, function(t) t[sp1, , drop = FALSE])
  attr(trait_list[[1]], "tree") <- tree_attr
  list(traits = trait_list[[1]], trait_list = trait_list, env = env, several = TRUE,
       niche_cor = niche_cor)
}

#' the environment of the sites on the axes of the niche, as a matrix with one named column per
#' set of traits, in their order when both are named
#' @noRd
evolving_niche_env <- function(env, traits, trait_list, landscape) {
  if (is.null(env)) env <- landscape$env_axes
  if (is.null(env)) {
    stop("simulate_evolving_metacommunity(): a list of traits needs `env`, one column per ",
         "element", call. = FALSE)
  }
  env <- as.matrix(env)
  bad <- !is.numeric(env) || nrow(env) != length(landscape$env)
  if (bad || ncol(env) != length(trait_list) || anyNA(env)) {
    stop("simulate_evolving_metacommunity(): `env` needs one row per site (",
         length(landscape$env), ") and one column per ",
         "set of traits (", length(trait_list), ")", call. = FALSE)
  }
  # named axes are matched by name
  named <- !is.null(colnames(env)) && !is.null(names(traits))
  if (named) {
    if (!setequal(colnames(env), names(traits)) || anyDuplicated(names(traits))) {
      stop("the names of `traits` and the columns of `env` must be the same axes",
           call. = FALSE)
    }
    env <- env[, names(traits), drop = FALSE]
  }
  colnames(env) <- colnames(env) %||% names(traits) %||% paste0("axis", seq_len(ncol(env)))
  env
}

#' stops unless the settings of `evolution` fit the niche and the mode of speciation
#' @noRd
evolving_check_evolution <- function(evolution, several, speciation) {
  if (!inherits(evolution, "mc_evolution")) {
    stop("simulate_evolving_metacommunity(): `evolution` must come from `evolution_control()`",
         call. = FALSE)
  }
  matrices <- is_covariance(evolution$genetic_var) || is_covariance(evolution$env_var)
  if (!several && matrices) {
    stop("simulate_evolving_metacommunity(): `evolution` with covariance matrices needs a ",
         "niche of several axes (a list of ",
         "`traits` and `env`), one trait per axis", call. = FALSE)
  }
  if (speciation$mode == "ecological" && all(unlist(evolution$genetic_var) == 0)) {
    stop("ecological speciation needs `evolution` with genetic variance", call. = FALSE)
  }
  invisible(evolution)
}

#' a root of the covariance matrix of the jumps of the traits at speciation, with several
#' traits (`jump_var` as a matrix, or one variance for all traits); NULL with one trait
#' @noRd
evolving_jump_root <- function(jump, several, env) {
  if (is.matrix(jump)) {
    if (!several) stop("a matrix `jump_var` needs several traits", call. = FALSE)
    jump <- evolution_matrices(jump, ncol(env), "x", "jump_var")[[1]]
  } else if (several) {
    jump <- diag(jump, ncol(env))
  }
  if (several) evolution_root(jump) else NULL
}

#' stops unless the species at the start, when several, come with an ultrametric tree with
#' branch lengths whose tips are those species
#' @noRd
evolving_check_tree <- function(tree0, traits) {
  k0 <- nrow(traits)
  if (k0 > 1 && is.null(tree0)) {
    stop("several species at the start need a tree in `traits` (from `species_traits()`)",
         call. = FALSE)
  }
  if (k0 > 1) {
    if (is.null(tree0$edge.length) || anyNA(tree0$edge.length)) {
      stop("the tree in `traits` needs branch lengths, in generations", call. = FALSE)
    }
    if (!setequal(tree0$tip.label, rownames(traits))) {
      stop("the tips of the tree in `traits` must be its species", call. = FALSE)
    }
    if (!ape::is.ultrametric(tree0, tol = 1e-6)) {
      stop("the tree in `traits` must be ultrametric", call. = FALSE)
    }
  }
  invisible(tree0)
}

#' runs the simulation of simulate_evolving_metacommunity(), whose arguments are already
#' checked: prepares the state, steps through time with speciation at each generation, and
#' returns the state (`st`), the records of the run (`bk`), the output of the dynamics (`out`)
#' and what the output of the simulation needs
#' @noRd
evolving_run <- function(landscape, niche, speciation, evolution, comp, n_steps, max_species,
                         dyn) {
  traits <- niche$traits
  jump_root <- evolving_jump_root(speciation$jump_var, niche$several, niche$env)
  tree0 <- attr(traits, "tree")
  evolving_check_tree(tree0, traits)
  if (n_steps %% evolution$steps_per_generation != 0) {
    stop("`n_steps` must be a multiple of `steps_per_generation`", call. = FALSE)
  }
  start <- evolving_start(landscape, niche, evolution, comp, n_steps, tree0, dyn)
  st <- start$st
  per <- evolution$steps_per_generation
  bk <- evolving_book(st, start$split0, list(speciation = speciation, comp = comp, per = per,
                                             max_species = max_species,
                                             jump_root = jump_root,
                                             disp = landscape$dispersal))
  for (step in seq_len(n_steps)) {
    mc_step(st, step)
    evolving_after_step(st, bk, step / per)
    if (step %% per == 0) evolving_generation(st, bk, step)
  }
  list(st = st, bk = bk, out = mc_finish(st), tree0 = tree0, niche = niche,
       n_steps = n_steps, generations = n_steps / per, dyn = dyn)
}

#' the state of mc_dynamics() at the start of simulate_evolving_metacommunity() (`st`), with
#' the times at which each pair of species at the start split (`split0`)
#' @noRd
evolving_start <- function(landscape, niche, evolution, comp, n_steps, tree0, dyn) {
  traits <- niche$traits
  k0 <- nrow(traits)
  sites <- paste0("site", seq_along(landscape$env))
  species0 <- rownames(traits)
  evolution$traits <- if (niche$several) {
    stats::setNames(niche$trait_list, colnames(niche$env))
  } else {
    traits
  }
  # a correlation between the axes of the niche (niche_correlation()) goes with the traits
  if (niche$several) attr(evolution$traits, "niche_cor") <- niche$niche_cor
  evolution$env <- if (niche$several) niche$env else landscape$env
  # time at which each pair of species at the start split, in generations before the start
  split0 <- matrix(0, k0, k0, dimnames = list(species0, species0))
  if (k0 > 1) split0 <- -ape::cophenetic.phylo(tree0)[species0, species0] / 2
  alpha0 <- evolving_alpha(split0, 0, comp$competition, comp$inter_max, comp$inter_min,
                           comp$scale, comp$distance_scale)
  # the state of mc_dynamics(), prepared by its own code, with no seeding after the start
  args <- c(list(r = matrix(1, length(sites), k0, dimnames = list(sites, species0)),
                 alpha = alpha0, n_steps = n_steps, seed_steps = numeric(0),
                 disp_matrix = landscape$dispersal, sites = sites, species = species0,
                 pool = NULL, evolution = evolution), dyn)
  st <- do.call(mc_state, args)
  if (comp$niche_overlap && comp$competition == "phylo") {
    st$competition_hook <- niche_crowding(st)
  }
  list(st = st, split0 = split0)
}

#' the records of a run of simulate_evolving_metacommunity(), kept in an environment that the
#' steps of the run change: the table of species, the split times, the speciation events, the
#' divergence among populations, the richness, the series and the changes of range; `cfg`
#' holds the settings of the run
#' @noRd
evolving_book <- function(st, split0, cfg) {
  bk <- new.env(parent = emptyenv())
  species0 <- rownames(split0)
  bk$cfg <- cfg
  bk$sp_table <- data.frame(species = species0, parent = NA_character_, birth = 0,
                            end = NA_real_, fate = "extant", order = 0L, completion = 0,
                            stringsAsFactors = FALSE)
  # times at which each pair of alive species split, by name
  bk$split_time <- split0
  bk$events <- data.frame(time = numeric(0), parent = character(0), new = character(0),
                          stringsAsFactors = FALSE)
  # per species (indexed over all species ever; NULL once extinct), for the geographic and
  # ecological modes: divergence or count among its populations, generations in a row in one
  # cluster, and expected migrants (from row to column) in the current generation
  bk$uses_div <- cfg$speciation$mode %in% c("geographic", "ecological")
  bk$divergence <- bk$together <- bk$flow <- list()
  if (bk$uses_div) {
    k0 <- length(species0)
    bk$divergence <- lapply(seq_len(k0), function(j) matrix(NA_real_, st$m, st$m))
    bk$together <- lapply(seq_len(k0), function(j) matrix(0, st$m, st$m))
    bk$flow <- lapply(seq_len(k0), function(j) matrix(0, st$m, st$m))
  }
  bk$richness <- bk$richness_good <- integer(0)
  # presence of every lineage ever (columns) in each site at the end of the last generation,
  # for the colonisations and extirpations of the series; rows of the series by generation
  bk$before <- evolving_presence(st)
  bk$series <- list()
  # changes of presence of each lineage in each site (generation, lineage index, site, present)
  bk$range_log <- list(range_changes(matrix(FALSE, st$m, st$k_all), bk$before, 0))
  bk
}

#' presence of every lineage ever (columns) in each site, from the state `st`
#' @noRd
evolving_presence <- function(st) {
  b <- matrix(FALSE, st$m, st$k_all)
  b[, st$alive] <- st$n > 0
  b
}

#' a matrix of presences `b` with columns added for the lineages born since it was made
#' @noRd
evolving_grow <- function(b, st) {
  if (ncol(b) < st$k_all) b <- cbind(b, matrix(FALSE, st$m, st$k_all - ncol(b)))
  b
}

#' records in `bk` the extinction at generation `gen` of the species `lost` (indices over all
#' species ever)
#' @noRd
evolving_mark_extinct <- function(bk, lost, gen) {
  bk$sp_table$end[lost] <- gen
  bk$sp_table$fate[lost] <- "extinct"
  keep <- !rownames(bk$split_time) %in% bk$sp_table$species[lost]
  bk$split_time <- bk$split_time[keep, keep, drop = FALSE]
  invisible(bk)
}

#' after each step of the dynamics: species lost at this step leave the state and the records,
#' and the flows of migrants and the populations lost in the step are recorded for the
#' divergence among populations
#' @noRd
evolving_after_step <- function(st, bk, gen) {
  # species lost at this step leave now (mc_step() would drop them at the next one, with no
  # random draws either way), so that the table, richness and output agree
  lost <- evolving_prune(st)
  if (length(lost)) {
    evolving_mark_extinct(bk, lost, gen)
    if (bk$uses_div) bk$divergence[lost] <- bk$together[lost] <- bk$flow[lost] <- list(NULL)
  }
  if (bk$uses_div && st$k > 0) {
    bk$flow <- evolving_flow(st, bk$flow, st$disp_matrix)
    # the lists leave `bk` while they change, so that they are changed in place
    divergence <- bk$divergence
    together <- bk$together
    bk$divergence <- bk$together <- NULL
    # a population lost at any step of the generation founds anew
    for (c in seq_len(st$k)) {
      gone <- st$n[, c] == 0
      if (any(gone)) {
        j <- st$alive[c]
        divergence[[j]][gone, ] <- NA_real_
        divergence[[j]][, gone] <- NA_real_
        together[[j]][gone, ] <- 0
        together[[j]][, gone] <- 0
      }
    }
    bk$divergence <- divergence
    bk$together <- together
  }
  invisible(bk)
}

#' the end of a generation: divergence among populations, speciation, competition among the
#' species alive, and the records of richness, series and ranges, in `st` and `bk`
#' @noRd
evolving_generation <- function(st, bk, step) {
  cfg <- bk$cfg
  gen <- step / cfg$per
  if (st$k == 0) return(evolving_generation_empty(st, bk, gen))
  if (bk$uses_div) {
    upd <- evolving_divergence(st, bk$divergence, bk$together, bk$flow, cfg$speciation)
    bk$divergence <- upd$divergence
    bk$together <- upd$together
  }
  buds <- evolving_initiate(st, bk$divergence, cfg$speciation, cfg$disp)
  if (length(buds)) buds <- evolving_speciate(st, bk, buds, gen)
  if (bk$uses_div) bk$flow <- lapply(bk$flow, function(x) if (!is.null(x)) x * 0)
  # records of this step, taken by mc_step() before the buddings, follow them
  if (length(buds)) evolving_record_again(st, step)
  if (cfg$comp$competition == "phylo" || nrow(bk$events)) {
    now <- st$species[st$alive]
    a <- evolving_alpha(bk$split_time[now, now, drop = FALSE], gen, cfg$comp$competition,
                        cfg$comp$inter_max, cfg$comp$inter_min, cfg$comp$scale,
                        cfg$comp$distance_scale)
    st$at <- t(a)
  }
  bk$richness <- c(bk$richness, length(st$alive))
  bk$richness_good <- c(bk$richness_good, sum(bk$sp_table$completion[st$alive] <= gen))
  bk$series[[length(bk$series) + 1]] <- evolving_series_rows(st, evolving_grow(bk$before, st),
                                                             buds, bk$sp_table, gen)
  after <- evolving_presence(st)
  bk$range_log[[length(bk$range_log) + 1]] <- range_changes(evolving_grow(bk$before, st),
                                                            after, gen)
  bk$before <- after
  invisible(bk)
}

#' the records of a generation that ends with no species alive
#' @noRd
evolving_generation_empty <- function(st, bk, gen) {
  bk$richness <- c(bk$richness, 0L)
  bk$richness_good <- c(bk$richness_good, 0L)
  if (any(bk$before)) {
    before <- evolving_grow(bk$before, st)
    bk$series[[length(bk$series) + 1]] <- evolving_series_rows(st, before, list(),
                                                               bk$sp_table, gen)
    bk$range_log[[length(bk$range_log) + 1]] <- range_changes(before,
                                                              matrix(FALSE, st$m, st$k_all),
                                                              gen)
    bk$before[] <- FALSE
  }
  invisible(bk)
}

#' the buddings of a generation: names the new species, splits them from their mothers in the
#' state `st` and records them in `bk`; returns `buds` with the names and indices of the new
#' species
#' @noRd
evolving_speciate <- function(st, bk, buds, gen) {
  speciation <- bk$cfg$speciation
  qs <- vapply(buds, function(b) length(b$sites), numeric(1))
  if (st$k_all + sum(qs) > bk$cfg$max_species) {
    stop("more than `max_species` (", bk$cfg$max_species, ") species; raise it or lower the ",
         "rate of speciation", call. = FALSE)
  }
  all_names <- evolving_names(st$species, sum(qs))
  last <- cumsum(qs)
  for (i in seq_along(buds)) {
    pos <- (last[i] - qs[i] + 1):last[i]
    buds[[i]]$names <- all_names[pos]
    buds[[i]]$ids <- st$k_all + pos
  }
  # every per-species structure grows once for all the buddings of the generation
  evolving_split_many(st, buds, speciation$jump_var, mutants = speciation$mode == "point",
                      jump_root = bk$cfg$jump_root)
  mothers <- bk$sp_table$species[vapply(buds, function(b) b$j, numeric(1))]
  done <- gen + if (speciation$completion == 0) rep(0, sum(qs)) else
    if (speciation$completion_dist == "fixed") rep(speciation$completion, sum(qs)) else
      stats::rexp(sum(qs), 1 / speciation$completion)
  bk$sp_table <- rbind(bk$sp_table, data.frame(species = all_names,
                                               parent = rep(mothers, qs), birth = gen,
                                               end = NA_real_, fate = "extant",
                                               order = sequence(qs), completion = done,
                                               stringsAsFactors = FALSE))
  bk$events <- rbind(bk$events, data.frame(
    time = gen, parent = mothers,
    new = vapply(buds, function(b) paste(b$names, collapse = ","), character(1)),
    stringsAsFactors = FALSE))
  bk$split_time <- evolving_split_time(bk$split_time, rep(mothers, qs), all_names, gen)
  if (bk$uses_div) evolving_split_divergence(st, bk, buds)
  # mothers whose last individuals became point mutants are gone (no other mode empties
  # a mother, and point mutation keeps no divergence)
  lost <- evolving_prune(st)
  if (length(lost)) evolving_mark_extinct(bk, lost, gen)
  buds
}

#' the divergence among the populations of each group that budded, which becomes a species of
#' its own (the first group keeps the mother), in the records `bk`
#' @noRd
evolving_split_divergence <- function(st, bk, buds) {
  for (b in buds) {
    j <- b$j
    d <- bk$divergence[[j]]
    tg <- bk$together[[j]]
    for (g in seq_along(b$groups)) {
      id <- if (g == 1) j else b$ids[g - 1]
      out <- !(seq_len(st$m) %in% b$groups[[g]])
      dg <- d
      dg[out, ] <- NA_real_
      dg[, out] <- NA_real_
      tgg <- tg
      tgg[out, ] <- 0
      tgg[, out] <- 0
      bk$divergence[[id]] <- dg
      bk$together[[id]] <- tgg
      bk$flow[[id]] <- matrix(0, st$m, st$m)
    }
  }
  invisible(bk)
}

#' writes again the records of a recorded step, taken by mc_step() before the buddings
#' @noRd
evolving_record_again <- function(st, step) {
  if (!is.null(st$history) && step %in% st$record_steps) {
    i <- match(step, st$record_steps)
    st$history[, st$alive, i] <- st$n
    if (isTRUE(st$ev$mv)) {
      st$z_history[, st$alive, , i] <- st$z[, st$alive, ]
    } else if (!is.null(st$z_history)) {
      st$z_history[, st$alive, i] <- st$z[, st$alive]
    }
  }
  invisible(st)
}

#' the output of simulate_evolving_metacommunity() from a run of evolving_run(): the extant
#' lineages are grouped into good species, with their abundances, traits, trees and records
#' @noRd
evolving_output <- function(run, landscape, speciation, comp, has_evolution) {
  st <- run$st
  out <- run$out
  dyn <- run$dyn
  sp_table <- run$bk$sp_table
  extant <- sp_table$species[sp_table$fate == "extant"]
  sp_table$end[sp_table$fate == "extant"] <- run$generations
  sp_table$status <- ifelse(sp_table$completion <= sp_table$end, "good", "incipient")
  tree_full <- evolving_tree(sp_table, run$tree0)
  # extant lineages grouped into good species, each represented in the tree by one lineage
  obs <- evolving_observed(sp_table)
  sp <- evolving_good_species(st, out, obs, extant)
  alpha <- evolving_alpha(run$bk$split_time[obs$rep, obs$rep, drop = FALSE], run$generations,
                          comp$competition, comp$inter_max, comp$inter_min, comp$scale,
                          comp$distance_scale)
  dimnames(alpha) <- list(sp$species, sp$species)
  by <- function(x, sum = FALSE) {
    evolving_by_species(x[, extant, drop = FALSE], sp$lineage_ab, sp$members, sum)
  }
  kernel_group <- dyn$kernel_group
  structure(list(abundance = by(out$abundance, TRUE),
                 tree = evolving_observed_tree(tree_full, obs, run$tree0, run$generations),
                 tree_full = tree_full,
                 species = sp_table, events = run$bk$events, lineages = sp$lineage_ab,
                 lineage_species = obs$group[extant], richness = run$bk$richness,
                 richness_good = run$bk$richness_good, genetic = sp$genetic,
                 traits = evolving_traits_out(st, sp), alpha = alpha, landscape = landscape,
                 dispersal = st$dispersal[sp$idx], extirpation = st$extirpation[sp$idx],
                 disturbance = dyn$disturbance, dispersal_time = dyn$dispersal_time,
                 kernels = st$kernels, kernel_step = dyn$kernel_step,
                 dispersal_site = dyn$dispersal_site,
                 disturbance_sensitivity = st$disturbance_sensitivity[sp$idx],
                 kernel_group = if (is.numeric(kernel_group)) st$kernel_group[sp$idx] else
                   colnames(dyn$kernel_step)[st$kernel_group[sp$idx]],
                 speciation = speciation,
                 evolution = if (has_evolution) evolving_state_out(st, out, sp, run$n_steps),
                 series = evolving_series_out(st, run$bk), ranges = evolving_ranges(st, run$bk),
                 n_steps = run$n_steps,
                 env = if (run$niche$several) as.data.frame(run$niche$env) else NULL,
                 generations = run$generations, history = evolving_history(st, out, sp, by),
                 arrival = by(out$arrival), propagule = by(out$propagule),
                 immigrants = by(out$immigrants), migrant_load = by(out$migrant_load)),
            class = "mc_evolving")
}

#' the good species at the end of a run and their extant lineages: the names of the `species`,
#' their indices in the state (`idx`), the lineages of each one (`members`), the abundances
#' (`lineage_ab`) and the genetic means (`genetic`) of the extant lineages
#' @noRd
evolving_good_species <- function(st, out, obs, extant) {
  species <- names(obs$rep)
  mv <- isTRUE(st$ev$mv)
  genetic <- if (mv) out$evolution$genetic[, extant, , drop = FALSE] else
    out$evolution$genetic[, extant, drop = FALSE]
  list(species = species, idx = match(species, st$species), extant = extant, mv = mv,
       members = split(extant, factor(obs$group[extant], species)),
       lineage_ab = out$abundance[, extant, drop = FALSE], genetic = genetic)
}

#' the tree of the good species at the end of a run: the full tree with one lineage per
#' species, or NULL with fewer than two species
#' @noRd
evolving_observed_tree <- function(tree_full, obs, tree0, generations) {
  species <- names(obs$rep)
  if (length(species) < 2) return(NULL)
  tr <- ape::keep.tip(tree_full, unname(obs$rep))
  tr$tip.label <- species[match(tr$tip.label, obs$rep)]
  # the stem of the reconstructed tree runs from the start to its crown
  if (is.null(tree0)) tr$root.edge <- generations - max(ape::node.depth.edgelength(tr))
  tr
}

#' the traits of the good species at the end of a run, with the mean genetic optimum over the
#' populations: a table, or one table per axis with several traits
#' @noRd
evolving_traits_out <- function(st, sp) {
  species <- sp$species
  if (!sp$mv) {
    optimum <- vapply(sp$members, function(l) mean(sp$genetic[, l], na.rm = TRUE), numeric(1))
    return(data.frame(optimum = unname(optimum), sigma = st$ev$sigma[sp$idx],
                      r_max = st$ev$r_max[sp$idx], row.names = species))
  }
  # one table of traits per axis, with the mean genetic optimum over the populations
  traits_out <- stats::setNames(lapply(seq_len(st$ev$q), function(a) {
    ga <- sp$genetic[, , a]
    ga <- matrix(ga, st$m, dimnames = list(NULL, sp$extant))
    optimum <- vapply(sp$members, function(l) mean(ga[, l], na.rm = TRUE), numeric(1))
    data.frame(optimum = unname(optimum), sigma = st$ev$sigma[sp$idx, a],
               r_max = st$ev$r_axis[sp$idx, a], row.names = species)
  }), st$ev$axes)
  # the correlation between the axes: the shared matrix, or that of each extant species
  if (is.list(st$ev$niche_cor)) {
    attr(traits_out, "niche_cor") <- stats::setNames(
      lapply(st$ev$sel[sp$idx], function(s) s$cor), species)
  } else if (!is.null(st$ev$niche_cor)) {
    attr(traits_out, "niche_cor") <- st$ev$niche_cor
  }
  traits_out
}

#' the abundances of the good species at the recorded steps (sites x species x steps), or NULL
#' without records
#' @noRd
evolving_history <- function(st, out, sp, by) {
  if (is.null(out$history)) return(NULL)
  history <- array(0L, c(st$m, length(sp$species), dim(out$history)[3]),
                   dimnames = list(rownames(out$history), sp$species,
                                   dimnames(out$history)[[3]]))
  for (i in seq_len(dim(out$history)[3])) {
    history[, , i] <- by(matrix(out$history[, , i], st$m,
                                dimnames = dimnames(out$history)[1:2]), TRUE)
  }
  history
}

#' the series of a run as one data frame, with one row per lineage and generation
#' @noRd
evolving_series_out <- function(st, bk) {
  series_df <- if (length(bk$series)) {
    do.call(rbind, bk$series)
  } else {
    evolving_series_rows(st, bk$before, list(), bk$sp_table, 0)[0, ]
  }
  rownames(series_df) <- NULL
  series_df
}

#' the changes of range of a run as a data frame of generation, lineage, site and presence
#' @noRd
evolving_ranges <- function(st, bk) {
  range_log <- do.call(rbind, bk$range_log)
  data.frame(generation = range_log[, 1], lineage = st$species[range_log[, 2]],
             site = as.integer(range_log[, 3]), present = range_log[, 4] == 1,
             stringsAsFactors = FALSE)
}

#' the evolutionary state by species, for sample_individuals(): genetic means weighted by the
#' abundance of the lineages, the plastic part and the variances among individuals
#' @noRd
evolving_state_out <- function(st, out, sp, n_steps) {
  if (sp$mv) {
    return(evolving_state_mv(st, out, sp$lineage_ab, sp$genetic, sp$members, sp$species,
                             sp$extant, sp$idx, n_steps))
  }
  species <- sp$species
  now <- evolving_species_state(sp$lineage_ab, sp$genetic, sp$members, species)
  wmean <- now$genetic
  plastic <- outer(st$ev$env(n_steps) - st$ev$ref, st$ev$b[sp$idx])
  plastic[is.na(wmean)] <- NA_real_
  dimnames(plastic) <- dimnames(wmean)
  hist <- evolving_state_history(st, out, sp)
  list(genetic = wmean, plastic = plastic, within_var = now$within,
       history = hist$z, within_history = hist$within,
       environment = st$ev$env(n_steps), environment_history = hist$env,
       since_change = out$evolution$since_change,
       reference = out$evolution$reference,
       genetic_var = stats::setNames(st$ev$g[sp$idx], species),
       env_var = stats::setNames(st$ev$ve[sp$idx], species),
       plasticity = stats::setNames(st$ev$b[sp$idx], species),
       settings = out$evolution$settings)
}

#' the evolutionary state by species at the recorded steps, for the lineages extant at the
#' end: the genetic means (`z`), the variance among lineages (`within`) and the environment
#' (`env`); NULL each without records
#' @noRd
evolving_state_history <- function(st, out, sp) {
  z_hist <- within_hist <- env_hist <- NULL
  if (!is.null(out$history) && !is.null(out$evolution$history)) {
    dims <- c(st$m, length(sp$species), dim(out$history)[3])
    names3 <- list(rownames(sp$lineage_ab), sp$species, dimnames(out$history)[[3]])
    z_hist <- within_hist <- array(NA_real_, dims, dimnames = names3)
    for (i in seq_len(dims[3])) {
      slice <- function(a) {
        matrix(a[, sp$extant, i], st$m, dimnames = list(rownames(sp$lineage_ab), sp$extant))
      }
      at <- evolving_species_state(slice(out$history), slice(out$evolution$history),
                                   sp$members, sp$species)
      z_hist[, , i] <- at$genetic
      within_hist[, , i] <- at$within
    }
    env_hist <- out$evolution$environment_history
  }
  list(z = z_hist, within = within_hist, env = env_hist)
}

#' @export
print.mc_evolving <- function(x, ...) {
  fates <- table(factor(x$species$fate, c("extant", "extinct")))
  cat("<mc_evolving>", nrow(x$abundance), "sites,", x$generations, "generations\n")
  cat("  lineages ever:", nrow(x$species), "| extant:", fates[["extant"]], "| extinct:",
      fates[["extinct"]], "| buddings:", nrow(x$events), "\n")
  if (x$speciation$completion > 0) cat("  good species with extant lineages:",
                                       ncol(x$abundance), "\n")
  cat("  speciation:", format(x$speciation$mode), "\n")
  invisible(x)
}

#' Convert to a simulated metacommunity
#'
#' Turns the output of [simulate_evolving_metacommunity()] into an `mc_sim` object, as from
#' [simulate_metacommunity()], with the extant species, their reconstructed tree and their
#' traits (the mean genetic optimum over their populations), for [as_mc_data()],
#' [summary_stats()] and the other functions of the package.
#'
#' @param x An `mc_evolving` object.
#' @param ... Unused.
#' @return An `mc_sim` object.
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' ls <- landscape_grid(5, disp_scale = 0.5)
#' ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
#' sim <- simulate_evolving_metacommunity(ls, ancestor, dispersal = 0.002, n_steps = 300,
#'                                        speciation = speciation_control("geographic", 80))
#' as_metacommunity(sim)
#' @export
as_metacommunity <- function(x, ...) UseMethod("as_metacommunity")

#' @export
as_metacommunity.mc_evolving <- function(x, ...) {
  traits <- x$traits
  attr(traits, "tree") <- x$tree
  rownames(x$abundance) <- paste0("site", seq_len(nrow(x$abundance)))
  fix <- function(y) {
    if (!is.null(y)) rownames(y) <- rownames(x$abundance)
    y
  }
  structure(list(abundance = x$abundance, history = x$history, landscape = x$landscape,
                 env = x$env, traits = traits, alpha = x$alpha, dispersal = x$dispersal,
                 extirpation = x$extirpation, pool = NULL, tree = x$tree,
                 arrival = fix(x$arrival), propagule = fix(x$propagule),
                 immigrants = fix(x$immigrants), migrant_load = fix(x$migrant_load),
                 evolution = x$evolution, n_steps = x$n_steps), class = "mc_sim")
}

#' competition among species from the times at which each pair split (generations; 0 on the
#' diagonal is ignored) at generation `now`
#' @noRd
evolving_alpha <- function(split_time, now, competition, inter_max, inter_min, scale,
                           distance_scale) {
  k <- nrow(split_time)
  a <- if (competition == "none") matrix(0, k, k) else
    inter_min + (inter_max - inter_min) * exp(-2 * (now - split_time) / distance_scale)
  diag(a) <- 1
  dimnames(a) <- dimnames(split_time)
  scale * a
}

#' crowding in each site with competition scaled by the overlap of the Gaussian niches of the
#' phenotypes p there (sites x alive species); the phylogenetic part is the alpha of the state
#' @noRd
niche_crowding <- function(st) {
  if (isTRUE(st$ev$mv) && isTRUE(st$ev$correlated)) return(niche_crowding_correlated(st))
  if (isTRUE(st$ev$mv)) {
    return(function(n, p) {
      sigma2 <- st$ev$sigma[st$alive, , drop = FALSE]^2
      out <- matrix(0, nrow(n), ncol(n))
      for (s in which(rowSums(n) > 0)) {
        here <- which(n[s, ] > 0)
        # the overlap of two Gaussian niches is the product over axes
        log_overlap <- 0
        for (a in seq_len(st$ev$q)) {
          w <- outer(sigma2[here, a], sigma2[here, a], "+")
          d2 <- outer(p[s, here, a], p[s, here, a], "-")^2
          # sigma = Inf gives no penalty; sigma = 0 on both sides gives overlap 1 at the same
          # phenotype and 0 elsewhere
          log_overlap <- log_overlap - ifelse(is.finite(w), niche_penalty(w, d2), 0)
        }
        out[s, here] <- n[s, here] %*% (st$at[here, here, drop = FALSE] * exp(log_overlap))
      }
      out
    })
  }
  function(n, p) {
    sigma2 <- st$ev$sigma[st$alive]^2
    out <- matrix(0, nrow(n), ncol(n))
    # only species present in a site need its crowding (the others have no one to grow)
    for (s in which(rowSums(n) > 0)) {
      here <- which(n[s, ] > 0)
      w <- outer(sigma2[here], sigma2[here], "+")
      d2 <- outer(p[s, here], p[s, here], "-")^2
      overlap <- exp(-niche_penalty(w, d2))
      out[s, here] <- n[s, here] %*% (st$at[here, here, drop = FALSE] * overlap)
    }
    out
  }
}

#' the same with a correlation between the axes of the niches (niche_correlation()): the overlap
#' of the niches of two species is exp(-d' M d / 2), with d the difference of their phenotypes
#' and M = (Omega_i + Omega_j)^-1. The breadths and correlations of a species never change, so M
#' is kept for each pair of species (indexed over all species ever)
#' @noRd
niche_crowding_correlated <- function(st) {
  cache <- new.env(parent = emptyenv())
  precision <- function(l) {
    key <- paste0("precision", l)
    if (is.null(cache[[key]])) {
      cor <- st$ev$sel[[l]]$cor %||% diag(st$ev$q)
      cache[[key]] <- niche_precision(st$ev$sigma[l, ], cor)
    }
    cache[[key]]
  }
  overlap_matrix <- function(i, j) {
    key <- paste(i, j)
    if (is.null(cache[[key]])) {
      cache[[key]] <- niche_overlap_matrix(precision(i), precision(j))
    }
    cache[[key]]
  }
  function(n, p) {
    alive <- st$alive
    at <- st$at
    out <- sweep(n, 2, diag(at), "*")
    here <- n > 0
    for (a in seq_len(ncol(n) - 1L)) {
      for (b in (a + 1L):ncol(n)) {
        both <- here[, a] & here[, b]
        if (!any(both)) next
        d <- matrix(p[both, a, ] - p[both, b, ], sum(both))
        overlap <- exp(-0.5 * rowSums((d %*% overlap_matrix(alive[a], alive[b])) * d))
        out[both, b] <- out[both, b] + n[both, a] * at[a, b] * overlap
        out[both, a] <- out[both, a] + n[both, b] * at[b, a] * overlap
      }
    }
    out
  }
}

#' split times (by name, among alive species) after the species `names_new` bud from `mothers`
#' (one per new species) at generation `gen`: a new species shares the history of its mother up
#' to `gen`, and new species of the same mother split from each other at `gen`
#' @noRd
evolving_split_time <- function(split_time, mothers, names_new, gen) {
  k <- nrow(split_time)
  anchor <- c(seq_len(k), match(mothers, rownames(split_time)))
  out <- split_time[anchor, anchor, drop = FALSE]
  new <- k + seq_along(names_new)
  out[cbind(new, anchor[new])] <- gen
  out[cbind(anchor[new], new)] <- gen
  same <- outer(anchor[new], anchor[new], "==")
  out[new, new][same] <- gen
  diag(out) <- 0
  dimnames(out) <- list(c(rownames(split_time), names_new), c(rownames(split_time), names_new))
  out
}

#' expected migrants of each species from each site (rows) to each other site (columns) in this
#' step, added to `flow` (indexed by species over all species ever)
#' @noRd
evolving_flow <- function(st, flow, disp_matrix) {
  em <- st$em
  if (is.null(em)) return(flow)
  for (c in seq_len(st$k)) {
    j <- st$alive[c]
    e <- em[, c]
    if (!any(e > 0)) next
    # implicit landscapes: emigrants are pooled and reach each other site with weight equal
    # to the emigrants of the others, so a site sends e / (m - 1) to each other site
    f <- if (!is.null(disp_matrix)) e * species_kernel(disp_matrix, j, st$kernel_group) else {
      out <- matrix(e / (length(e) - 1), length(e), length(e))
      diag(out) <- 0
      out
    }
    # with immigrant selection, only the immigrants that survive in the site they reach
    if (isTRUE(st$ev$imm_sel)) f <- f * evolving_immigrant_fitness(st, c, j)
    flow[[j]] <- flow[[j]] + f
  }
  flow
}

#' relative fitness (sources x sites) of the emigrants of the alive column c (species j) in
#' each site, as in immigrants_selected() and mv_immigrants_selected()
#' @noRd
evolving_immigrant_fitness <- function(st, c, j) {
  ev <- st$ev
  if (isTRUE(ev$mv)) {
    p <- matrix(st$pa_em[, c, ], st$m)
    p[is.na(p)] <- 0
    return(mv_fitness(ev, p, st$e_now, j))
  }
  p <- st$pa_em[, c]
  p[is.na(p)] <- 0
  width <- ev$g[j] + ev$sigma[j]^2 + ev$ve[j]
  if (!is.finite(width)) return(matrix(1, st$m, st$m))
  exp(-outer(p, st$e_now, "-")^2 / (2 * width))
}

#' migration rate between pairs of populations: the larger of the migrants each sends to the
#' other over the size of the receiver
#' @noRd
migration_rate <- function(f, n) {
  into <- sweep(f, 2, n, "/")
  into[!is.finite(into)] <- 0
  pmax(into, t(into))
}

#' divergence (geographic) or counts (ecological) among the populations of each species after
#' one generation, with the generations in a row in one cluster
#' @noRd
evolving_divergence <- function(st, divergence, together, flow, speciation) {
  # the matrices of each species change in this one frame, where they are copied only once
  for (c in seq_len(st$k)) {
    j <- st$alive[c]
    here <- st$n[, c] > 0
    d <- divergence[[j]]
    tg <- together[[j]]
    f <- flow[[j]]
    known <- !is.na(diag(d))
    # new populations take the divergences of the population that sent them most migrants
    for (a in which(here & !known)) {
      src <- which(known & here)
      if (length(src)) {
        b <- divergence_source(st, c, j, src, a, f)
        d[a, ] <- d[b, ]
        d[, a] <- d[, b]
        d[a, b] <- d[b, a] <- 0
        tg[a, ] <- tg[b, ]
        tg[, a] <- tg[, b]
        tg[a, b] <- tg[b, a] <- 0
      } else {
        d[a, ] <- d[, a] <- NA_real_
        d[a, here & !is.na(diag(d))] <- 0
        d[here & !is.na(diag(d)), a] <- 0
      }
      d[a, a] <- 0
      known[a] <- TRUE
    }
    d[!here, ] <- NA_real_
    d[, !here] <- NA_real_
    idx <- which(here)
    if (length(idx) > 1) {
      upd <- divergence_pairs(st, c, j, idx, d, tg, f, speciation)
      d[idx, idx] <- upd$sub
      if (!is.null(upd$tsub)) tg[idx, idx] <- upd$tsub
    } else if (length(idx) == 1) {
      d[idx, idx] <- 0
    }
    divergence[[j]] <- d
    together[[j]] <- tg
  }
  list(divergence = divergence, together = together)
}

#' the population, among those in `src`, from which the new population `a` of species j (alive
#' column c) takes its divergences: the one that sent it most migrants (`f`); without recorded
#' flow, the largest one weighted by the kernel toward `a`
#' @noRd
divergence_source <- function(st, c, j, src, a, f) {
  weight <- f[src, a]
  if (!any(weight > 0)) {
    kern <- species_kernel(st$disp_matrix, j, st$kernel_group)
    toward <- if (is.null(kern)) 1 else kern[src, a]
    weight <- st$n[src, c] * toward
  }
  src[which.max(weight)]
}

#' the divergences (`sub`) and, when they change, the generations together (`tsub`) among the
#' populations `idx` (two or more) of species j (alive column c) after one generation, by the
#' mode of speciation
#' @noRd
divergence_pairs <- function(st, c, j, idx, d, tg, f, speciation) {
  out <- list()
  sub <- d[idx, idx]
  rate <- migration_rate(f[idx, idx, drop = FALSE], st$n[idx, c])
  if (speciation$mode == "ecological") {
    far <- if (isTRUE(st$ev$mv)) {
      mv_distance(st$ev, j, matrix(st$z[idx, j, ], length(idx)))
    } else {
      zj <- st$z[idx, j]
      abs(outer(zj, zj, "-"))
    }
    apart <- far >= speciation$threshold & rate < speciation$migration_max
    sub[] <- ifelse(apart, sub + 1, 0)
  } else if (speciation$divergence == "cluster") {
    cl <- link_groups(rate >= speciation$migration_min)
    same <- outer(cl, cl, "==")
    tsub <- tg[idx, idx]
    tsub[] <- ifelse(same, tsub + 1, 0)
    kept <- ifelse(tsub >= speciation$reset_generations, 0, sub)
    sub[] <- ifelse(same, kept, sub + 1)
    out$tsub <- tsub
  } else {
    M <- f[idx, idx, drop = FALSE] + t(f[idx, idx, drop = FALSE])
    sub[] <- pmax(0, sub + 1 - speciation$homogenization * M)
  }
  diag(sub) <- 0
  out$sub <- sub
  out
}

#' groups of the populations of species j (column index over all species ever) by single
#' linkage on pairs below the limit; a list of sites
#' @noRd
evolving_groups <- function(st, j, d, speciation) {
  if (!is.finite(speciation$threshold)) return(list())
  c <- match(j, st$alive)
  idx <- which(st$n[, c] > 0)
  if (length(idx) < 2) return(list())
  limit <- if (speciation$mode == "geographic") speciation$threshold else
    speciation$min_generations
  sub <- d[idx, idx]
  sub[is.na(sub)] <- 0
  cl <- link_groups(sub < limit)
  if (length(unique(cl)) < 2) return(list())
  unname(split(idx, cl))
}

#' groups in the order of budding: first the one that keeps the species (most sites, then most
#' individuals), then the others from the most to the least divergent from it
#' @noRd
evolving_order <- function(st, j, groups, d, speciation) {
  c <- match(j, st$alive)
  size <- vapply(groups, length, numeric(1))
  total <- vapply(groups, function(s) sum(st$n[s, c]), numeric(1))
  keep <- order(-size, -total)[1]
  away <- function(s) {
    if (speciation$mode == "geographic") {
      mean(d[s, groups[[keep]], drop = FALSE], na.rm = TRUE)
    } else if (isTRUE(st$ev$mv)) {
      centre <- function(rows) colMeans(matrix(st$z[rows, j, ], length(rows)), na.rm = TRUE)
      mv_distance(st$ev, j, rbind(centre(s), centre(groups[[keep]])))[1, 2]
    } else {
      abs(mean(st$z[s, j], na.rm = TRUE) - mean(st$z[groups[[keep]], j], na.rm = TRUE))
    }
  }
  others <- setdiff(seq_along(groups), keep)
  dist <- vapply(groups[others], away, numeric(1))
  c(groups[keep], groups[others[order(-dist)]])
}

#' connected components of a symmetric logical adjacency matrix
#' @noRd
link_groups <- function(adj) {
  n <- nrow(adj)
  group <- integer(n)
  g <- 0L
  for (i in seq_len(n)) {
    if (group[i] > 0) next
    g <- g + 1L
    queue <- i
    group[i] <- g
    while (length(queue)) {
      v <- queue[1]
      queue <- queue[-1]
      nb <- which((adj[v, ] | adj[, v]) & group == 0)
      group[nb] <- g
      queue <- c(queue, nb)
    }
  }
  group
}

#' new species to start this generation: a list with, for each budding mother, `j` (index over
#' all species ever), `sites` (a list of the sites of each new species), `amounts` (individuals
#' each new species takes in its sites, or NULL for whole populations) and `groups` (for the
#' geographic and ecological modes, the groups of populations with the one kept first)
#' @noRd
evolving_initiate <- function(st, divergence, speciation, disp_matrix) {
  mode <- speciation$mode
  if (mode %in% c("geographic", "ecological")) {
    return(initiate_groups(st, divergence, speciation))
  }
  if (mode == "point") return(initiate_point(st, speciation))
  initiate_fission(st, speciation, disp_matrix)
}

#' buddings of the geographic and ecological modes: the species whose populations form two or
#' more groups, with the group that keeps the species first
#' @noRd
initiate_groups <- function(st, divergence, speciation) {
  out <- list()
  for (j in st$alive) {
    groups <- evolving_groups(st, j, divergence[[j]], speciation)
    if (length(groups) < 2) next
    groups <- evolving_order(st, j, groups, divergence[[j]], speciation)
    out[[length(out) + 1]] <- list(j = j, sites = groups[-1], amounts = NULL,
                                   groups = groups)
  }
  out
}

#' buddings by point mutation: each individual starts a species with the rate of speciation,
#' alone in its site
#' @noRd
initiate_point <- function(st, speciation) {
  out <- list()
  if (speciation$rate == 0) return(out)
  mut <- matrix(stats::rbinom(st$m * st$k, st$n, speciation$rate), st$m, st$k)
  for (c in which(colSums(mut) > 0)) {
    s <- rep(seq_len(st$m), mut[, c])
    out[[length(out) + 1]] <- list(j = st$alive[c], sites = as.list(s),
                                   amounts = as.list(rep(1L, length(s))))
  }
  out
}

#' buddings by fission: each species of two or more populations splits in two groups of sites
#' with the rate of speciation
#' @noRd
initiate_fission <- function(st, speciation, disp_matrix) {
  out <- list()
  p <- 1 - exp(-speciation$rate)
  if (p == 0) return(out)
  for (c in seq_len(st$k)) {
    idx <- which(st$n[, c] > 0)
    if (length(idx) < 2 || stats::runif(1) >= p) next
    j <- st$alive[c]
    groups <- evolving_order(st, j, evolving_partition(idx, speciation$partition,
                                                       disp_matrix), NULL, speciation)
    out[[length(out) + 1]] <- list(j = j, sites = groups[-1], amounts = NULL)
  }
  out
}

#' two groups of the occupied sites `idx` for fission: the sites nearest (by the kernel) to a
#' random one, in a random number, or a random partition
#' @noRd
evolving_partition <- function(idx, partition, disp_matrix) {
  q <- length(idx)
  if (partition == "spatial" && !is.null(disp_matrix)) {
    seed <- idx[sample.int(q, 1)]
    others <- setdiff(idx, seed)
    near <- c(seed, others[order(-disp_matrix[seed, others])])
    size <- sample.int(q - 1, 1)
    return(list(sort(near[-seq_len(size)]), sort(near[seq_len(size)])))
  }
  repeat {
    side <- stats::runif(q) < 0.5
    if (any(side) && !all(side)) break
  }
  list(idx[!side], idx[side])
}

#' several mothers bud at once: each element of `buds` has the mother `j` (index over all
#' species ever), the `sites` and `amounts` of each new species (amounts NULL for whole
#' populations), and their `names` and `ids`. Every per-species structure of the state grows once;
#' new species take the genetic means and history of the mother in their sites, shifted by one
#' Gaussian draw of variance `jump_var` each. Point mutants (`mutants = TRUE`) arrive now, from a
#' propagule of their own individuals.
#' @noRd
evolving_split_many <- function(st, buds, jump_var = 0, mutants = FALSE, jump_root = NULL) {
  q <- sum(vapply(buds, function(b) length(b$sites), numeric(1)))
  m <- st$m
  mothers <- rep(vapply(buds, function(b) b$j, numeric(1)),
                 vapply(buds, function(b) length(b$sites), numeric(1)))
  names_new <- unlist(lapply(buds, function(b) b$names))
  cols_j <- match(vapply(buds, function(b) b$j, numeric(1)), st$alive)
  evolving_grow_state(st, mothers, names_new, q)
  evolving_grow_species(st, mothers, names_new)
  # local copies in the frame `w`, modified in place and written back once
  fields <- c("arrival", "propagule", "immigrants", "load", "present", "z", "ph", "pool_mean")
  w <- list2env(mget(fields, envir = st), parent = parent.env(environment()))
  list2env(list(optimum = st$ev$optimum, n = st$n, new_n = matrix(0L, m, q), col = 0, m = m,
                mv = isTRUE(st$ev$mv), mutants = mutants, jump_var = jump_var,
                jump_root = jump_root, step = st$step), envir = w)
  for (bi in seq_along(buds)) {
    w$b <- buds[[bi]]
    w$c_j <- cols_j[bi]
    evolving_bud(w)
  }
  for (nm in fields) assign(nm, w[[nm]], envir = st)
  st$ev$optimum <- w$optimum
  st$n <- w$n
  new_ids <- unlist(lapply(buds, function(b) b$ids))
  st$n <- cbind(st$n, w$new_n)
  # emigrants of this step stay aligned with the alive species (none for the new ones)
  if (!is.null(st$em)) st$em <- cbind(st$em, matrix(0, m, q))
  # and so do their phenotypes, kept with immigrant selection
  if (!is.null(st$pa_em)) st$pa_em <- grow_species_array(st$pa_em, NA_real_, q)
  st$alive <- c(st$alive, new_ids)
  st$k_all <- st$k_all + q
  st$k <- length(st$alive)
  st$p_em <- rep(st$dispersal[st$alive], each = m)
  st$p_ex <- rep(st$extirpation[st$alive], each = m)
  st$at <- matrix(0, st$k, st$k)
  invisible(st)
}

#' a sites x species (x traits) array of the state with `q` more species, filled with `fill`
#' @noRd
grow_species_array <- function(x, fill, q) {
  if (length(dim(x)) == 2) return(cbind(x, matrix(fill, nrow(x), q)))
  d <- dim(x)
  out <- array(fill, c(d[1], d[2] + q, d[3]))
  out[, seq_len(d[2]), ] <- x
  out
}

#' a recorded history (sites x species x steps, or sites x species x traits x steps) with the
#' species `names_new` added, filled with `fill`; NULL without records
#' @noRd
grow_species_history <- function(x, fill, names_new) {
  if (is.null(x)) return(NULL)
  d <- dim(x)
  nd <- length(d)
  out <- array(fill, replace(d, 2, d[2] + length(names_new)))
  if (nd == 3) out[, seq_len(d[2]), ] <- x else out[, seq_len(d[2]), , ] <- x
  dn <- dimnames(x)
  dn[[2]] <- c(dn[[2]], names_new)
  dimnames(out) <- dn
  out
}

#' grows by `q` species the arrays of sites x species of the state `st`: the colonisation
#' history, the evolutionary state and the recorded histories
#' @noRd
evolving_grow_state <- function(st, mothers, names_new, q) {
  grow_cols <- function(x, fill) cbind(x, matrix(fill, nrow(x), q))
  for (nm in c("arrival", "propagule")) assign(nm, grow_cols(get(nm, st), NA_integer_), st)
  for (nm in c("immigrants", "load")) assign(nm, grow_cols(get(nm, st), 0), st)
  st$present <- grow_cols(st$present, FALSE)
  st$z <- grow_species_array(st$z, NA_real_, q)
  st$ph <- grow_species_array(st$ph, NA_real_, q)
  st$pool_mean <- grow_species_array(st$pool_mean, NA_real_, q)
  st$history <- grow_species_history(st$history, 0L, names_new)
  st$z_history <- grow_species_history(st$z_history, NA_real_, names_new)
  invisible(st)
}

#' adds to the per-species settings of the state `st` the new species `names_new`, each with
#' the values of its mother (`mothers`: indices over all species ever)
#' @noRd
evolving_grow_species <- function(st, mothers, names_new) {
  grow_vec <- function(x) c(x, x[mothers])
  grow_rows <- function(x) rbind(x, x[mothers, , drop = FALSE])
  st$dispersal <- grow_vec(st$dispersal)
  st$extirpation <- grow_vec(st$extirpation)
  st$disturbance_sensitivity <- grow_vec(st$disturbance_sensitivity)
  # routes of dispersal by group of species: new species stay in the group of their mother
  st$kernel_group <- grow_vec(st$kernel_group)
  # emigration through time (steps x species): new species take that of their mother
  if (!is.null(st$disp_prob_time)) {
    st$disp_prob_time <- cbind(st$disp_prob_time, st$disp_prob_time[, mothers, drop = FALSE])
  }
  if (isTRUE(st$ev$mv)) {
    for (nm in c("optimum", "sigma", "b", "r_axis")) st$ev[[nm]] <- grow_rows(st$ev[[nm]])
    for (nm in c("r_max", "fixed")) st$ev[[nm]] <- grow_vec(st$ev[[nm]])
    for (nm in c("G", "E", "sel", "L", "whiten")) {
      st$ev[[nm]] <- c(st$ev[[nm]], st$ev[[nm]][mothers])
    }
  } else {
    for (nm in c("g", "b", "ve", "optimum", "sigma", "r_max")) {
      st$ev[[nm]] <- grow_vec(st$ev[[nm]])
    }
  }
  st$fixed_z <- grow_vec(st$fixed_z)
  st$species <- c(st$species, names_new)
  invisible(st)
}

#' one mother buds, in the working frame `w` of evolving_split_many(): each new species of the
#' bud `b` takes its individuals and the history of the mother (alive column `c_j`) in its
#' sites, and the mother keeps what is left
#' @noRd
evolving_bud <- function(w) {
  w$j <- w$b$j
  w$taken <- integer(w$m)
  for (g in seq_along(w$b$sites)) {
    w$g <- g
    evolving_bud_take(w)
    evolving_bud_jump(w)
  }
  eval(quote({
    # the mother keeps what is left; where nothing is left, its population is gone
    n[, c_j] <- n[, c_j] - as.integer(taken)
    emptied <- which(taken > 0 & n[, c_j] == 0)
    arrival[emptied, j] <- NA
    propagule[emptied, j] <- NA
    if (mv) {
      z[emptied, j, ] <- NA
      ph[emptied, j, ] <- NA
    } else {
      z[emptied, j] <- NA
      ph[emptied, j] <- NA
    }
    immigrants[emptied, j] <- 0
    load[emptied, j] <- 0
    present[emptied, j] <- FALSE
  }), w)
  invisible(w)
}

#' the new species `g` of the bud `b` takes its individuals, and the colonisation history and
#' the genetic means of the mother in its sites, in the working frame `w`
#' @noRd
evolving_bud_take <- function(w) {
  eval(quote({
    s <- b$sites[[g]]
    id <- b$ids[g]
    col <- col + 1
    a <- if (is.null(b$amounts)) n[s, c_j] else b$amounts[[g]]
    new_n[s, col] <- a
    taken[s] <- taken[s] + a
    arrival[s, id] <- arrival[s, j]
    propagule[s, id] <- propagule[s, j]
    immigrants[s, id] <- immigrants[s, j]
    load[s, id] <- load[s, j]
    present[s, id] <- present[s, j]
    if (mv) {
      pool_mean[, id, ] <- pool_mean[, j, ]
      z[s, id, ] <- z[s, j, ]
      ph[s, id, ] <- ph[s, j, ]
    } else {
      pool_mean[, id] <- pool_mean[, j]
      z[s, id] <- z[s, j]
      ph[s, id] <- ph[s, j]
    }
    if (mutants) {
      arrival[s, id] <- as.integer(step)
      propagule[s, id] <- as.integer(a)
      immigrants[s, id] <- 0
      load[s, id] <- 0
    }
  }), w)
  invisible(w)
}

#' the jump of the traits of the new species `id` at its birth: one Gaussian draw of variance
#' `jump_var` (one multivariate draw with several traits), added to its genetic means,
#' phenotypes and optimum in the working frame `w`
#' @noRd
evolving_bud_jump <- function(w) {
  eval(quote({
    if (mv && !is.null(jump_root) && any(jump_root != 0)) {
      # one multivariate draw per new species, added to every trait of its populations
      shift <- drop(jump_root %*% stats::rnorm(ncol(jump_root)))
      for (tr in seq_along(shift)) {
        z[s, id, tr] <- z[s, id, tr] + shift[tr]
        ph[s, id, tr] <- ph[s, id, tr] + shift[tr]
        pool_mean[, id, tr] <- pool_mean[, id, tr] + shift[tr]
      }
      optimum[id, ] <- optimum[id, ] + shift
    } else if (!mv && jump_var > 0) {
      shift <- stats::rnorm(1, 0, sqrt(jump_var))
      z[s, id] <- z[s, id] + shift
      ph[s, id] <- ph[s, id] + shift
      pool_mean[, id] <- pool_mean[, id] + shift
      optimum[id] <- optimum[id] + shift
    }
  }), w)
  invisible(w)
}

#' the full tree (with extinct species) from the table of species: each species is a lineage
#' from its birth to its end, from which the species that budded from it branch off at their
#' birth, in their order; species at the start hang from the tips of `tree0` (or from the root,
#' alone). Built from the youngest lineages to the oldest, without recursion, since a lineage
#' may have many buds.
#' @noRd
evolving_tree <- function(sp, tree0) {
  buds <- split(seq_len(nrow(sp)), factor(sp$parent, sp$species))
  inner <- character(nrow(sp))
  stem <- numeric(nrow(sp))
  for (i in rev(seq_len(nrow(sp)))) {
    b <- buds[[i]]
    b <- b[order(sp$birth[b], sp$order[b])]
    cur <- sp$species[i]
    times <- c(sp$birth[i], sp$birth[b])
    len <- sp$end[i] - times[length(times)]
    for (k in rev(seq_along(b))) {
      cur <- sprintf("(%s:%.10g,%s:%.10g)", inner[b[k]], stem[b[k]], cur, len)
      len <- times[k + 1] - times[k]
      inner[b[k]] <- ""
    }
    inner[i] <- cur
    stem[i] <- len
  }
  roots <- which(is.na(sp$parent))
  if (length(roots) == 1) {
    tr <- ape::read.tree(text = sprintf("(%s:%.10g);", inner[roots], stem[roots]))
    return(if (length(tr$tip.label) > 1) ape::collapse.singles(tr, root.edge = TRUE) else tr)
  }
  ntip <- length(tree0$tip.label)
  node_newick <- function(v, len) {
    if (v <= ntip) {
      i <- match(tree0$tip.label[v], sp$species)
      return(sprintf("%s:%.10g", inner[i], stem[i] + len))
    }
    kids <- tree0$edge[tree0$edge[, 1] == v, 2]
    lens <- tree0$edge.length[tree0$edge[, 1] == v]
    sprintf("(%s):%.10g", paste(mapply(node_newick, kids, lens), collapse = ","), len)
  }
  text <- sub(":0$", "", node_newick(ntip + 1L, 0))
  ape::read.tree(text = paste0(text, ";"))
}

#' the evolutionary state by species with several traits: the genetic means and variance among
#' lineages of evolving_species_state(), trait by trait, as arrays of sites x species x traits,
#' at the end and at the recorded steps, with the plastic part and the matrices of each species
#' @noRd
evolving_state_mv <- function(st, out, lineage_ab, genetic, members, species, extant, sp_idx,
                              n_steps) {
  ev <- st$ev
  q <- ev$q
  sites <- rownames(lineage_ab)
  by_trait <- function(ab, z) {
    parts <- lapply(seq_len(q), function(a) {
      evolving_species_state(ab, matrix(z[, , a], st$m, dimnames = list(sites, extant)),
                             members, species)
    })
    shape <- function(what) {
      array(vapply(parts, function(p) p[[what]], matrix(0, st$m, length(species))),
            c(st$m, length(species), q), dimnames = list(sites, species, ev$axes))
    }
    list(genetic = shape("genetic"), within = shape("within"))
  }
  now <- by_trait(lineage_ab, genetic)
  e_end <- ev$env(n_steps)
  plastic <- now$genetic
  for (a in seq_len(q)) plastic[, , a] <- outer(e_end[, a] - ev$ref[a], ev$b[sp_idx, a])
  plastic[is.na(now$genetic)] <- NA_real_
  z_hist <- within_hist <- env_hist <- NULL
  if (!is.null(out$history) && !is.null(out$evolution$history)) {
    steps <- dimnames(out$history)[[3]]
    dims <- c(st$m, length(species), q, length(steps))
    names4 <- list(sites, species, ev$axes, steps)
    z_hist <- within_hist <- array(NA_real_, dims, dimnames = names4)
    for (i in seq_along(steps)) {
      ab <- matrix(out$history[, extant, i], st$m, dimnames = list(sites, extant))
      z <- out$evolution$history[, extant, , i, drop = FALSE]
      z <- array(z, c(st$m, length(extant), q))
      at <- by_trait(ab, z)
      z_hist[, , , i] <- at$genetic
      within_hist[, , , i] <- at$within
    }
    env_hist <- out$evolution$environment_history
    dimnames(env_hist)[[1]] <- sites
  }
  environment <- matrix(e_end, st$m, dimnames = list(sites, ev$axes))
  list(genetic = now$genetic, plastic = plastic, within_var = now$within, history = z_hist,
       within_history = within_hist, environment = environment,
       environment_history = env_hist, since_change = out$evolution$since_change,
       reference = stats::setNames(ev$ref, ev$axes),
       genetic_var = stats::setNames(ev$G[sp_idx], species),
       env_var = stats::setNames(ev$E[sp_idx], species),
       plasticity = matrix(ev$b[sp_idx, , drop = FALSE], length(species), length(ev$axes),
                           dimnames = list(species, ev$axes)),
       settings = out$evolution$settings)
}

#' genetic means of species (weighted by the abundance of their lineages) and the variance among
#' their lineages within each site (0 for species of one lineage), from lineage abundances and
#' genetic means (sites x extant lineages)
#' @noRd
evolving_species_state <- function(ab, z, members, species) {
  wmean <- vapply(members, function(l) {
    w <- ab[, l, drop = FALSE]
    zz <- z[, l, drop = FALSE]
    w[is.na(zz)] <- 0
    zz[w == 0] <- 0
    ifelse(rowSums(w) > 0, rowSums(w * zz) / rowSums(w), NA_real_)
  }, numeric(nrow(z)))
  wmean <- matrix(wmean, nrow(z), length(species), dimnames = list(rownames(ab), species))
  within <- vapply(seq_along(members), function(i) {
    l <- members[[i]]
    w <- ab[, l, drop = FALSE]
    dev <- sweep(z[, l, drop = FALSE], 1, wmean[, i])
    w[is.na(dev)] <- 0
    dev[w == 0] <- 0
    ifelse(rowSums(w) > 0, rowSums(w * dev^2) / rowSums(w), NA_real_)
  }, numeric(nrow(z)))
  list(genetic = wmean,
       within = matrix(within, nrow(z), length(species), dimnames = dimnames(wmean)))
}

#' changes of presence between two matrices of sites x lineages (all lineages ever), as rows of
#' generation, lineage, site and new presence (1 or 0)
#' @noRd
range_changes <- function(was, now, gen) {
  at <- which(was != now, arr.ind = TRUE)
  if (!length(at)) return(matrix(numeric(0), 0, 4))
  cbind(gen, at[, 2], at[, 1], now[at])
}

#' rows of the series for generation `gen`: for each lineage alive at its start, its occupancy,
#' abundance and mean genetic optimum (weighted by abundance) at its end, and the sites it
#' colonised and lost since the last generation, out of those it did not and did occupy then.
#' Lineages that died out in this generation get a row with their last losses. A site a mother
#' passed to its buds in this generation counts neither as a loss nor as exposure; a lineage
#' born in this generation has no transitions yet.
#' @noRd
evolving_series_rows <- function(st, before, buds, sp_table, gen) {
  dead <- setdiff(which(colSums(before) > 0), st$alive)
  ids <- c(st$alive, dead)
  n <- cbind(st$n, matrix(0, st$m, length(dead)))
  now <- n > 0
  was <- before[, ids, drop = FALSE]
  skip <- matrix(FALSE, st$m, length(ids))
  # a mother was alive at the start of the generation, so it is among `ids`
  for (b in buds) {
    c_j <- match(b$j, ids)
    given <- unique(unlist(b$sites))
    given <- given[!now[given, c_j]]
    skip[given, c_j] <- TRUE
  }
  newborn <- sp_table$birth[ids] == gen & gen > 0
  abundance <- colSums(n)
  # mean genetic optimum of each lineage, weighted by abundance (one column per trait)
  lineage_mean <- function(z) {
    zz <- cbind(z, matrix(0, st$m, length(dead)))
    zz[!now] <- 0
    ifelse(abundance > 0, colSums(n * zz) / abundance, NA_real_)
  }
  genetic <- if (isTRUE(st$ev$mv)) {
    g <- lapply(seq_len(st$ev$q), function(a) lineage_mean(st$z[, st$alive, a]))
    stats::setNames(g, paste0("genetic_", st$ev$axes))
  } else {
    list(genetic = lineage_mean(st$z[, st$alive, drop = FALSE]))
  }
  keep <- !skip
  out <- data.frame(generation = rep(gen, length(ids)), lineage = st$species[ids],
                    occupancy = colSums(now), abundance = abundance, genetic,
                    colonised = colSums(now & !was & keep),
                    extirpated = colSums(!now & was & keep),
                    present_before = colSums(was & keep), absent_before = colSums(!was & keep),
                    stringsAsFactors = FALSE)
  out[newborn, c("colonised", "extirpated", "present_before", "absent_before")] <- NA
  out
}

#' extant lineages grouped into good species: a lineage that completed by its end (`status`
#' "good") is its own species; an incipient one belongs to the species its parent belonged to at
#' its birth, as in the labels of the PBD; `group` maps lineages to species and
#' `rep` names, for each species (in the order of the table), the lineage that represents it:
#' itself if extant, otherwise its earliest extant member
#' @noRd
evolving_observed <- function(sp) {
  good <- sp$status == "good"
  extant <- which(sp$fate == "extant")
  owner <- vapply(extant, function(i) {
    if (good[i]) return(i)
    # as in the PBD: the species of the parent at the birth of the lineage
    t <- sp$birth[i]
    p <- match(sp$parent[i], sp$species)
    while (sp$completion[p] > t) {
      t <- sp$birth[p]
      p <- match(sp$parent[p], sp$species)
    }
    p
  }, numeric(1))
  group <- stats::setNames(sp$species[owner], sp$species[extant])
  ids <- sort(unique(owner))
  rep <- vapply(ids, function(o) {
    if (sp$fate[o] == "extant") return(sp$species[o])
    m <- extant[owner == o]
    sp$species[m[order(sp$birth[m])][1]]
  }, character(1))
  list(group = group, rep = stats::setNames(rep, sp$species[ids]))
}

#' a record of the extant lineages (sites x lineages) by good species: summed, or taken from the
#' most abundant lineage of the species in each site
#' @noRd
evolving_by_species <- function(x, abundance, members, sum = FALSE) {
  out <- vapply(members, function(l) {
    if (length(l) == 1) return(as.numeric(x[, l]))
    if (sum) return(rowSums(x[, l, drop = FALSE]))
    top <- l[max.col(abundance[, l, drop = FALSE], ties.method = "first")]
    as.numeric(x[cbind(seq_len(nrow(x)), match(top, colnames(x)))])
  }, numeric(nrow(x)))
  out <- matrix(out, nrow(x), length(members), dimnames = list(rownames(x), names(members)))
  if (is.integer(x)) storage.mode(out) <- "integer"
  out
}

# the state of mc_dynamics() after its preparation: a frame with the same arguments
mc_state <- mc_dynamics
body(mc_state) <- quote({
  mc_setup(environment())
  environment()
})

#' names for q new species that no species has used
#' @noRd
evolving_names <- function(used, q) {
  out <- character(0)
  i <- length(used)
  while (length(out) < q) {
    i <- i + 1
    nm <- paste0("sp", i)
    if (!nm %in% c(used, out)) out <- c(out, nm)
  }
  out
}

#' drops from the state the alive species with no individual left, as mc_step() does at the
#' start of the next step; returns their indices over all species ever
#' @noRd
evolving_prune <- function(st) {
  if (st$k == 0) return(integer(0))
  empty <- colSums(st$n) == 0
  if (!any(empty)) return(integer(0))
  lost <- st$alive[empty]
  keep <- !empty
  st$alive <- st$alive[keep]
  st$n <- st$n[, keep, drop = FALSE]
  st$at <- st$at[keep, keep, drop = FALSE]
  if (!is.null(st$em)) st$em <- st$em[, keep, drop = FALSE]
  if (!is.null(st$pa_em)) {
    st$pa_em <- if (length(dim(st$pa_em)) == 3) st$pa_em[, keep, , drop = FALSE] else
      st$pa_em[, keep, drop = FALSE]
  }
  st$k <- length(st$alive)
  st$p_em <- rep(st$dispersal[st$alive], each = st$m)
  st$p_ex <- rep(st$extirpation[st$alive], each = st$m)
  lost
}
