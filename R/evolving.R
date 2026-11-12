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
#'   difference of genetic means in units of the environment, that separates groups; `Inf`
#'   turns speciation off.
#' @param rate For `"point"`: probability per individual per generation; for `"fission"`: rate
#'   per species per generation.
#' @param divergence For `"geographic"`: `"cluster"` or `"gradual"`.
#' @param migration_min For `divergence = "cluster"`: migration rate per generation that links
#'   two populations. Where clusters form depends on dispersal, the kernel and abundances: the
#'   rate of the strongest link of each population should straddle it. In the design of
#'   `data-raw/evolving_migration_min.R` (strongest links with median 0.006), species did not
#'   split at all at 0.002 or below and split almost as with no links at 0.02 or above.
#' @param reset_generations For `divergence = "cluster"`: generations in a row in the same
#'   cluster after which divergence returns to 0.
#' @param homogenization For `divergence = "gradual"`: reduction of divergence per migrant.
#' @param min_generations For `"ecological"`: generations in a row of difference and reduced
#'   gene flow before groups split.
#' @param migration_max For `"ecological"`: migration rate per generation below which gene flow
#'   counts as reduced.
#' @param partition For `"fission"`: `"spatial"` or `"random"`.
#' @param jump_var Variance of the shift of the genetic mean of each new species.
#' @param completion Time to complete speciation, in generations (fixed or mean); 0 for
#'   instantaneous speciation.
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
  one <- function(x) is.numeric(x) && length(x) == 1
  if (mode %in% c("geographic", "ecological")) {
    if (missing(threshold)) stop("`threshold` is needed", call. = FALSE)
    stopifnot(one(threshold), threshold > 0)
    rate <- NA_real_
  } else {
    if (missing(rate)) stop("`rate` is needed", call. = FALSE)
    stopifnot(one(rate), rate >= 0, mode == "fission" || rate <= 1)
    threshold <- NA_real_
  }
  stopifnot(one(migration_min), migration_min > 0, one(reset_generations),
            reset_generations >= 1, one(homogenization), homogenization >= 0,
            one(min_generations), min_generations >= 1, one(migration_max), migration_max >= 0,
            one(jump_var), jump_var >= 0, one(completion), completion >= 0)
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
  if (x$jump_var > 0) cat(" | jump_var =", format(x$jump_var))
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
#'   species).
#' @param speciation Settings from [speciation_control()].
#' @param evolution Settings from [evolution_control()]; `NULL` keeps the optimum of each
#'   species fixed (`genetic_var = 0`), which rules out ecological speciation.
#' @param competition `"phylo"` or `"none"`.
#' @param inter_max,inter_min,scale,distance_scale Competition between species, relative to
#'   competition within a species: between sisters just split (`inter_max`) and between
#'   distant species (`inter_min`); the overall scale (competition within a species); and the
#'   distance in generations over which the excess over `inter_min` decays by a factor of
#'   \eqn{e}{e}.
#' @param niche_overlap Should competition between species also fall with the difference of
#'   their phenotypes in each site?
#' @param initial Optional matrix (sites x species) of abundances at the start; by default
#'   Poisson(`seed_mean`) individuals of each species in every site.
#' @param max_species Maximum number of species ever, including extinct ones; the simulation
#'   stops with an error beyond it.
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
#'   * `genetic`: genetic means of the populations of the extant lineages.
#'   * `traits`: of the extant species, with the mean genetic optimum over the populations of
#'     their lineages as `optimum`.
#'   * `alpha`: among the extant species at the end, without the niche overlap, which varies
#'     among sites.
#'   * The records and colonisation history of [simulate_metacommunity()] for the extant
#'     species: records summed over their extant lineages (lineages extinct by the end are
#'     left out); colonisation, in each site, of the most abundant lineage.
#'   * `series`: one row per lineage alive in each generation, with its `occupancy` (sites),
#'     regional `abundance`, `genetic` mean (weighted by abundance), and the sites it
#'     `colonised` and lost (`extirpated`) since the previous generation, out of those where it
#'     was absent (`absent_before`) and present (`present_before`). Sites that a lineage passed
#'     to a new one are not counted as lost, and the transitions are `NA` in the generation of
#'     birth. [persistence()] turns them into rates of colonisation and extirpation.
#'   * `evolution`: with `evolution`, the state of the extant species at the end, for
#'     [sample_individuals()] with `source = "evolution"`: genetic means in each site
#'     (`genetic`, weighted over their lineages by abundance; `NA` where absent), the plastic
#'     part of the phenotype (`plastic`), the variance among the lineages of a species in
#'     each site (`within_var`, 0 for species of one lineage), the environment of each site
#'     at the end (`environment`), the reference environment of plasticity (`reference`), the genetic and environmental variances and the slope of plasticity of
#'     each species (`genetic_var`, `env_var`, `plasticity`) and the `settings`; otherwise
#'     `NULL`.
#'   * The `landscape`, `dispersal`, `extirpation`, `speciation` and `n_steps` used.
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
                                            max_species = 500) {
  stopifnot(inherits(landscape, "mc_landscape"))
  if (!inherits(speciation, "mc_speciation")) {
    stop("`speciation` must come from `speciation_control()`", call. = FALSE)
  }
  competition <- match.arg(competition)
  immigration <- match.arg(immigration)
  stopifnot(inter_max >= 0, inter_min >= 0, inter_min <= inter_max, scale > 0,
            distance_scale > 0, max_species >= 1, is.logical(niche_overlap))
  traits <- as_traits(traits)
  k0 <- nrow(traits)
  has_evolution <- !is.null(evolution)
  if (!has_evolution) evolution <- evolution_control(0, ne_ratio = Inf)
  if (!inherits(evolution, "mc_evolution")) {
    stop("`evolution` must come from `evolution_control()`", call. = FALSE)
  }
  if (speciation$mode == "ecological" && all(evolution$genetic_var == 0)) {
    stop("ecological speciation needs `evolution` with genetic variance", call. = FALSE)
  }
  tree0 <- attr(traits, "tree")
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
  if (n_steps %% evolution$steps_per_generation != 0) {
    stop("`n_steps` must be a multiple of `steps_per_generation`", call. = FALSE)
  }
  sites <- paste0("site", seq_along(landscape$env))
  species0 <- rownames(traits)
  evolution$traits <- traits
  evolution$env <- landscape$env
  # time at which each pair of species at the start split, in generations before the start
  split0 <- matrix(0, k0, k0, dimnames = list(species0, species0))
  if (k0 > 1) split0 <- -ape::cophenetic.phylo(tree0)[species0, species0] / 2
  alpha0 <- evolving_alpha(split0, 0, competition, inter_max, inter_min, scale,
                           distance_scale)

  # the state of mc_dynamics(), prepared by its own code, with no seeding after the start
  args <- list(r = matrix(1, length(sites), k0, dimnames = list(sites, species0)),
               alpha = alpha0, dispersal = dispersal, n_steps = n_steps,
               seed_steps = numeric(0), seed_mean = seed_mean,
               disp_matrix = landscape$dispersal, record_steps = record_steps, sites = sites,
               species = species0, extirpation = extirpation, pool = NULL, initial = initial,
               immigration = immigration, evolution = evolution)
  st <- do.call(mc_state, args)
  if (niche_overlap && competition == "phylo") st$competition_hook <- niche_crowding(st)

  per <- evolution$steps_per_generation
  generations <- n_steps / per
  sp_table <- data.frame(species = species0, parent = NA_character_, birth = 0,
                         end = NA_real_, fate = "extant", order = 0L, completion = 0,
                         stringsAsFactors = FALSE)
  # times at which each pair of alive species split, by name
  split_time <- split0
  events <- data.frame(time = numeric(0), parent = character(0), new = character(0),
                       stringsAsFactors = FALSE)
  # per species (indexed over all species ever; NULL once extinct), for the geographic and
  # ecological modes: divergence or count among its populations, generations in a row in one
  # cluster, and expected migrants (from row to column) in the current generation
  uses_div <- speciation$mode %in% c("geographic", "ecological")
  empty <- function() matrix(NA_real_, st$m, st$m)
  divergence <- together <- flow <- list()
  if (uses_div) {
    divergence <- lapply(seq_len(k0), function(j) empty())
    together <- lapply(seq_len(k0), function(j) matrix(0, st$m, st$m))
    flow <- lapply(seq_len(k0), function(j) matrix(0, st$m, st$m))
  }
  richness <- richness_good <- integer(0)
  # presence of every lineage at the end of the last generation (sites x all lineages ever), for
  # the colonisations and extirpations of the series; rows of the series by generation
  # presence of every lineage ever (columns) in each site at the end of the last generation
  grow <- function(b) {
    if (ncol(b) < st$k_all) b <- cbind(b, matrix(FALSE, st$m, st$k_all - ncol(b)))
    b
  }
  now_present <- function() {
    b <- matrix(FALSE, st$m, st$k_all)
    b[, st$alive] <- st$n > 0
    b
  }
  before <- now_present()
  series <- list()
  for (step in seq_len(n_steps)) {
    mc_step(st, step)
    gen <- step / per
    # species lost at this step leave now (mc_step() would drop them at the next one, with no
    # random draws either way), so that the table, richness and output agree
    lost <- evolving_prune(st)
    if (length(lost)) {
      sp_table$end[lost] <- gen
      sp_table$fate[lost] <- "extinct"
      keep <- !rownames(split_time) %in% sp_table$species[lost]
      split_time <- split_time[keep, keep, drop = FALSE]
      if (uses_div) divergence[lost] <- together[lost] <- flow[lost] <- list(NULL)
    }
    if (uses_div && st$k > 0) {
      flow <- evolving_flow(st, flow, landscape$dispersal)
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
    }
    if (step %% per != 0) next
    if (st$k == 0) {
      richness <- c(richness, 0L)
      richness_good <- c(richness_good, 0L)
      if (any(before)) {
        series[[length(series) + 1]] <- evolving_series_rows(st, grow(before), list(),
                                                             sp_table, gen)
        before[] <- FALSE
      }
      next
    }
    if (uses_div) {
      upd <- evolving_divergence(st, divergence, together, flow, speciation)
      divergence <- upd$divergence
      together <- upd$together
    }
    buds <- evolving_initiate(st, divergence, speciation, landscape$dispersal)
    if (length(buds)) {
      qs <- vapply(buds, function(b) length(b$sites), numeric(1))
      if (st$k_all + sum(qs) > max_species) {
        stop("more than `max_species` (", max_species, ") species; raise it or lower the ",
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
      evolving_split_many(st, buds, speciation$jump_var, mutants = speciation$mode == "point")
      mothers <- sp_table$species[vapply(buds, function(b) b$j, numeric(1))]
      done <- gen + if (speciation$completion == 0) rep(0, sum(qs)) else
        if (speciation$completion_dist == "fixed") rep(speciation$completion, sum(qs)) else
          stats::rexp(sum(qs), 1 / speciation$completion)
      sp_table <- rbind(sp_table, data.frame(species = all_names,
                                             parent = rep(mothers, qs), birth = gen,
                                             end = NA_real_, fate = "extant",
                                             order = sequence(qs), completion = done,
                                             stringsAsFactors = FALSE))
      events <- rbind(events, data.frame(
        time = gen, parent = mothers,
        new = vapply(buds, function(b) paste(b$names, collapse = ","), character(1)),
        stringsAsFactors = FALSE))
      split_time <- evolving_split_time(split_time, rep(mothers, qs), all_names, gen)
      for (b in buds) {
        j <- b$j
        if (uses_div) {
          d <- divergence[[j]]
          tg <- together[[j]]
          for (g in seq_along(b$groups)) {
            id <- if (g == 1) j else b$ids[g - 1]
            out <- !(seq_len(st$m) %in% b$groups[[g]])
            dg <- d
            dg[out, ] <- NA_real_
            dg[, out] <- NA_real_
            tgg <- tg
            tgg[out, ] <- 0
            tgg[, out] <- 0
            divergence[[id]] <- dg
            together[[id]] <- tgg
            flow[[id]] <- matrix(0, st$m, st$m)
          }
        }
      }
      # mothers whose last individuals became point mutants are gone (no other mode empties
      # a mother, and point mutation keeps no divergence)
      lost <- evolving_prune(st)
      if (length(lost)) {
        sp_table$end[lost] <- gen
        sp_table$fate[lost] <- "extinct"
        keep <- !rownames(split_time) %in% sp_table$species[lost]
        split_time <- split_time[keep, keep, drop = FALSE]
      }
    }
    if (uses_div) flow <- lapply(flow, function(x) if (!is.null(x)) x * 0)
    # records of this step, taken by mc_step() before the buddings, follow them
    if (length(buds) && !is.null(st$history) && step %in% st$record_steps) {
      i <- match(step, st$record_steps)
      st$history[, st$alive, i] <- st$n
      if (!is.null(st$z_history)) st$z_history[, st$alive, i] <- st$z[, st$alive]
    }
    if (competition == "phylo" || nrow(events)) {
      now <- st$species[st$alive]
      a <- evolving_alpha(split_time[now, now, drop = FALSE], gen, competition, inter_max,
                          inter_min, scale, distance_scale)
      st$at <- t(a)
    }
    richness <- c(richness, length(st$alive))
    richness_good <- c(richness_good, sum(sp_table$completion[st$alive] <= gen))
    series[[length(series) + 1]] <- evolving_series_rows(st, grow(before), buds, sp_table,
                                                         gen)
    before <- now_present()
  }
  out <- mc_finish(st)
  extant <- sp_table$species[sp_table$fate == "extant"]
  sp_table$end[sp_table$fate == "extant"] <- generations
  sp_table$status <- ifelse(sp_table$completion <= sp_table$end, "good", "incipient")
  tree_full <- evolving_tree(sp_table, tree0)
  # extant lineages grouped into good species, each represented in the tree by one lineage
  obs <- evolving_observed(sp_table)
  species <- names(obs$rep)
  tree <- if (length(species) >= 2) {
    tr <- ape::keep.tip(tree_full, unname(obs$rep))
    tr$tip.label <- species[match(tr$tip.label, obs$rep)]
    # the stem of the reconstructed tree runs from the start to its crown
    if (is.null(tree0)) tr$root.edge <- generations - max(ape::node.depth.edgelength(tr))
    tr
  } else {
    NULL
  }
  lineage_ab <- out$abundance[, extant, drop = FALSE]
  genetic <- out$evolution$genetic[, extant, drop = FALSE]
  members <- split(extant, factor(obs$group[extant], species))
  optimum <- vapply(members, function(l) mean(genetic[, l], na.rm = TRUE), numeric(1))
  sp_idx <- match(species, st$species)
  traits_out <- data.frame(optimum = unname(optimum), sigma = st$ev$sigma[sp_idx],
                           r_max = st$ev$r_max[sp_idx], row.names = species)
  alpha <- evolving_alpha(split_time[obs$rep, obs$rep, drop = FALSE], generations,
                          competition, inter_max, inter_min, scale, distance_scale)
  dimnames(alpha) <- list(species, species)
  by <- function(x, sum = FALSE) evolving_by_species(x[, extant, drop = FALSE], lineage_ab,
                                                     members, sum)
  history <- NULL
  if (!is.null(out$history)) {
    history <- array(0L, c(st$m, length(species), dim(out$history)[3]),
                     dimnames = list(rownames(out$history), species,
                                     dimnames(out$history)[[3]]))
    for (i in seq_len(dim(out$history)[3])) history[, , i] <- by(out$history[, , i], TRUE)
  }
  series_df <- if (length(series)) do.call(rbind, series) else evolving_series_rows(st, before,
                                                                                   list(), sp_table, 0)[0, ]
  rownames(series_df) <- NULL
  # the evolutionary state by species, for sample_individuals(): genetic means weighted by the
  # abundance of the lineages, the plastic part and the variances among individuals
  ev_out <- NULL
  if (has_evolution) {
    wmean <- vapply(members, function(l) {
      w <- lineage_ab[, l, drop = FALSE]
      zz <- genetic[, l, drop = FALSE]
      zz[w == 0] <- 0
      ifelse(rowSums(w) > 0, rowSums(w * zz) / rowSums(w), NA_real_)
    }, numeric(st$m))
    wmean <- matrix(wmean, st$m, length(species), dimnames = list(rownames(lineage_ab), species))
    # variance among the lineages of a species within each site (weighted by abundance), 0 for
    # species of one lineage
    within <- vapply(seq_along(members), function(i) {
      l <- members[[i]]
      w <- lineage_ab[, l, drop = FALSE]
      dev <- sweep(genetic[, l, drop = FALSE], 1, wmean[, i])
      dev[w == 0] <- 0
      ifelse(rowSums(w) > 0, rowSums(w * dev^2) / rowSums(w), NA_real_)
    }, numeric(st$m))
    within <- matrix(within, st$m, length(species), dimnames = dimnames(wmean))
    plastic <- outer(st$ev$env(n_steps) - st$ev$ref, st$ev$b[sp_idx])
    plastic[is.na(wmean)] <- NA_real_
    dimnames(plastic) <- dimnames(wmean)
    ev_out <- list(genetic = wmean, plastic = plastic, within_var = within,
                   environment = st$ev$env(n_steps), reference = out$evolution$reference,
                   genetic_var = stats::setNames(st$ev$g[sp_idx], species),
                   env_var = stats::setNames(st$ev$ve[sp_idx], species),
                   plasticity = stats::setNames(st$ev$b[sp_idx], species),
                   settings = out$evolution$settings)
  }
  structure(list(abundance = by(out$abundance, TRUE), tree = tree, tree_full = tree_full,
                 species = sp_table, events = events, lineages = lineage_ab,
                 lineage_species = obs$group[extant], richness = richness,
                 richness_good = richness_good, genetic = genetic, traits = traits_out,
                 alpha = alpha, landscape = landscape, dispersal = st$dispersal[sp_idx],
                 extirpation = st$extirpation[sp_idx], speciation = speciation,
                 evolution = ev_out, series = series_df, n_steps = n_steps,
                 generations = generations, history = history,
                 arrival = by(out$arrival), propagule = by(out$propagule),
                 immigrants = by(out$immigrants), migrant_load = by(out$migrant_load)),
            class = "mc_evolving")
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
                 env = NULL, traits = traits, alpha = x$alpha, dispersal = x$dispersal,
                 extirpation = x$extirpation, pool = NULL, tree = x$tree,
                 arrival = fix(x$arrival), propagule = fix(x$propagule),
                 immigrants = fix(x$immigrants), migrant_load = fix(x$migrant_load),
                 evolution = x$evolution, n_steps = x$n_steps), class = "mc_sim")
}

# competition among species from the times at which each pair split (generations; 0 on the
# diagonal is ignored) at generation `now`
evolving_alpha <- function(split_time, now, competition, inter_max, inter_min, scale,
                           distance_scale) {
  k <- nrow(split_time)
  a <- if (competition == "none") matrix(0, k, k) else
    inter_min + (inter_max - inter_min) * exp(-2 * (now - split_time) / distance_scale)
  diag(a) <- 1
  dimnames(a) <- dimnames(split_time)
  scale * a
}

# crowding in each site with competition scaled by the overlap of the Gaussian niches of the
# phenotypes p there (sites x alive species); the phylogenetic part is the alpha of the state
niche_crowding <- function(st) {
  function(n, p) {
    sigma2 <- st$ev$sigma[st$alive]^2
    out <- matrix(0, nrow(n), ncol(n))
    # only species present in a site need its crowding (the others have no one to grow)
    for (s in which(rowSums(n) > 0)) {
      here <- which(n[s, ] > 0)
      overlap <- exp(-outer(p[s, here], p[s, here], "-")^2 /
                       (2 * outer(sigma2[here], sigma2[here], "+")))
      out[s, here] <- n[s, here] %*% (st$at[here, here, drop = FALSE] * overlap)
    }
    out
  }
}

# split times (by name, among alive species) after the species `names_new` bud from `mothers`
# (one per new species) at generation `gen`: a new species shares the history of its mother up
# to `gen`, and new species of the same mother split from each other at `gen`
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

# expected migrants of each species from each site (rows) to each other site (columns) in this
# step, added to `flow` (indexed by species over all species ever)
evolving_flow <- function(st, flow, disp_matrix) {
  em <- st$em
  if (is.null(em)) return(flow)
  for (c in seq_len(st$k)) {
    j <- st$alive[c]
    e <- em[, c]
    if (!any(e > 0)) next
    # implicit landscapes: emigrants are pooled and reach each other site with weight equal
    # to the emigrants of the others, so a site sends e / (m - 1) to each other site
    f <- if (!is.null(disp_matrix)) e * disp_matrix else {
      out <- matrix(e / (length(e) - 1), length(e), length(e))
      diag(out) <- 0
      out
    }
    flow[[j]] <- flow[[j]] + f
  }
  flow
}

# migration rate between pairs of populations: the larger of the migrants each sends to the
# other over the size of the receiver
migration_rate <- function(f, n) {
  into <- sweep(f, 2, n, "/")
  into[!is.finite(into)] <- 0
  pmax(into, t(into))
}

# divergence (geographic) or counts (ecological) among the populations of each species after
# one generation, with the generations in a row in one cluster
evolving_divergence <- function(st, divergence, together, flow, speciation) {
  for (c in seq_len(st$k)) {
    j <- st$alive[c]
    here <- st$n[, c] > 0
    d <- divergence[[j]]
    tg <- together[[j]]
    f <- flow[[j]]
    known <- !is.na(diag(d))
    # new populations take the divergences of the population that sent them most migrants;
    # without recorded flow, the largest population weighted by the kernel toward them
    for (a in which(here & !known)) {
      src <- which(known & here)
      if (length(src)) {
        weight <- f[src, a]
        if (!any(weight > 0)) {
          weight <- st$n[src, c] *
            (if (!is.null(st$disp_matrix)) st$disp_matrix[src, a] else 1)
        }
        b <- src[which.max(weight)]
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
      sub <- d[idx, idx]
      rate <- migration_rate(f[idx, idx, drop = FALSE], st$n[idx, c])
      if (speciation$mode == "ecological") {
        zj <- st$z[idx, j]
        apart <- abs(outer(zj, zj, "-")) >= speciation$threshold &
          rate < speciation$migration_max
        sub[] <- ifelse(apart, sub + 1, 0)
      } else if (speciation$divergence == "cluster") {
        cl <- link_groups(rate >= speciation$migration_min)
        same <- outer(cl, cl, "==")
        tsub <- tg[idx, idx]
        tsub[] <- ifelse(same, tsub + 1, 0)
        sub[] <- ifelse(same, ifelse(tsub >= speciation$reset_generations, 0, sub), sub + 1)
        tg[idx, idx] <- tsub
      } else {
        M <- f[idx, idx, drop = FALSE] + t(f[idx, idx, drop = FALSE])
        sub[] <- pmax(0, sub + 1 - speciation$homogenization * M)
      }
      diag(sub) <- 0
      d[idx, idx] <- sub
    } else if (length(idx) == 1) {
      d[idx, idx] <- 0
    }
    divergence[[j]] <- d
    together[[j]] <- tg
  }
  list(divergence = divergence, together = together)
}

# groups of the populations of species j (column index over all species ever) by single
# linkage on pairs below the limit; a list of sites
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

# groups in the order of budding: first the one that keeps the species (most sites, then most
# individuals), then the others from the most to the least divergent from it
evolving_order <- function(st, j, groups, d, speciation) {
  c <- match(j, st$alive)
  size <- vapply(groups, length, numeric(1))
  total <- vapply(groups, function(s) sum(st$n[s, c]), numeric(1))
  keep <- order(-size, -total)[1]
  away <- function(s) {
    if (speciation$mode == "geographic") {
      mean(d[s, groups[[keep]], drop = FALSE], na.rm = TRUE)
    } else {
      abs(mean(st$z[s, j], na.rm = TRUE) - mean(st$z[groups[[keep]], j], na.rm = TRUE))
    }
  }
  others <- setdiff(seq_along(groups), keep)
  dist <- vapply(groups[others], away, numeric(1))
  c(groups[keep], groups[others[order(-dist)]])
}

# connected components of a symmetric logical adjacency matrix
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

# new species to start this generation: a list with, for each budding mother, `j` (index over
# all species ever), `sites` (a list of the sites of each new species), `amounts` (individuals
# each new species takes in its sites, or NULL for whole populations) and `groups` (for the
# geographic and ecological modes, the groups of populations with the one kept first)
evolving_initiate <- function(st, divergence, speciation, disp_matrix) {
  mode <- speciation$mode
  out <- list()
  if (mode %in% c("geographic", "ecological")) {
    for (j in st$alive) {
      groups <- evolving_groups(st, j, divergence[[j]], speciation)
      if (length(groups) < 2) next
      groups <- evolving_order(st, j, groups, divergence[[j]], speciation)
      out[[length(out) + 1]] <- list(j = j, sites = groups[-1], amounts = NULL,
                                     groups = groups)
    }
  } else if (mode == "point") {
    if (speciation$rate == 0) return(out)
    mut <- matrix(stats::rbinom(st$m * st$k, st$n, speciation$rate), st$m, st$k)
    for (c in which(colSums(mut) > 0)) {
      s <- rep(seq_len(st$m), mut[, c])
      out[[length(out) + 1]] <- list(j = st$alive[c], sites = as.list(s),
                                     amounts = as.list(rep(1L, length(s))))
    }
  } else {
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
  }
  out
}

# two groups of the occupied sites `idx` for fission: the sites nearest (by the kernel) to a
# random one, in a random number, or a random partition
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

# species j (index over all species ever) buds: new species g takes `amounts[[g]]` individuals
# (all, when NULL) in the sites `sites[[g]]`; see evolving_split_many()
evolving_split <- function(st, j, sites, names_new, amounts = NULL, jump_var = 0,
                           mutants = FALSE) {
  evolving_split_many(st, list(list(j = j, sites = sites, amounts = amounts, names = names_new,
                                    ids = st$k_all + seq_along(sites))), jump_var, mutants)
}

# several mothers bud at once: each element of `buds` has the mother `j` (index over all
# species ever), the `sites` and `amounts` of each new species (amounts NULL for whole
# populations), and their `names` and `ids`. Every per-species structure of the state grows once;
# new species take the genetic means and history of the mother in their sites, shifted by one
# Gaussian draw of variance `jump_var` each. Point mutants (`mutants = TRUE`) arrive now, from a
# propagule of their own individuals.
evolving_split_many <- function(st, buds, jump_var = 0, mutants = FALSE) {
  q <- sum(vapply(buds, function(b) length(b$sites), numeric(1)))
  m <- st$m
  mothers <- rep(vapply(buds, function(b) b$j, numeric(1)),
                 vapply(buds, function(b) length(b$sites), numeric(1)))
  names_new <- unlist(lapply(buds, function(b) b$names))
  grow_cols <- function(x, fill) cbind(x, matrix(fill, nrow(x), q))
  grow_vec <- function(x) c(x, x[mothers])
  grow_hist <- function(x, fill) {
    if (is.null(x)) return(NULL)
    d <- dim(x)
    out <- array(fill, c(d[1], d[2] + q, d[3]))
    out[, seq_len(d[2]), ] <- x
    dimnames(out) <- list(dimnames(x)[[1]], c(dimnames(x)[[2]], names_new), dimnames(x)[[3]])
    out
  }
  cols_j <- match(vapply(buds, function(b) b$j, numeric(1)), st$alive)
  for (nm in c("arrival", "propagule")) assign(nm, grow_cols(get(nm, st), NA_integer_), st)
  for (nm in c("immigrants", "load")) assign(nm, grow_cols(get(nm, st), 0), st)
  st$present <- grow_cols(st$present, FALSE)
  st$z <- grow_cols(st$z, NA_real_)
  st$ph <- grow_cols(st$ph, NA_real_)
  st$pool_mean <- grow_cols(st$pool_mean, NA_real_)
  st$history <- grow_hist(st$history, 0L)
  st$z_history <- grow_hist(st$z_history, NA_real_)
  st$dispersal <- grow_vec(st$dispersal)
  st$extirpation <- grow_vec(st$extirpation)
  for (nm in c("g", "b", "ve", "optimum", "sigma", "r_max")) st$ev[[nm]] <- grow_vec(st$ev[[nm]])
  st$fixed_z <- grow_vec(st$fixed_z)
  st$species <- c(st$species, names_new)
  new_n <- matrix(0L, m, q)
  # local copies, modified in place and written back once
  fields <- c("arrival", "propagule", "immigrants", "load", "present", "z", "ph", "pool_mean")
  x <- mget(fields, envir = st)
  optimum <- st$ev$optimum
  n <- st$n
  col <- 0
  for (bi in seq_along(buds)) {
    b <- buds[[bi]]
    j <- b$j
    c_j <- cols_j[bi]
    taken <- integer(m)
    for (g in seq_along(b$sites)) {
      s <- b$sites[[g]]
      id <- b$ids[g]
      col <- col + 1
      a <- if (is.null(b$amounts)) n[s, c_j] else b$amounts[[g]]
      new_n[s, col] <- a
      taken[s] <- taken[s] + a
      x$pool_mean[, id] <- x$pool_mean[, j]
      for (nm in c("arrival", "propagule", "immigrants", "load", "present", "z", "ph")) {
        x[[nm]][s, id] <- x[[nm]][s, j]
      }
      if (mutants) {
        x$arrival[s, id] <- as.integer(st$step)
        x$propagule[s, id] <- as.integer(a)
        x$immigrants[s, id] <- 0
        x$load[s, id] <- 0
      }
      if (jump_var > 0) {
        shift <- stats::rnorm(1, 0, sqrt(jump_var))
        x$z[s, id] <- x$z[s, id] + shift
        x$ph[s, id] <- x$ph[s, id] + shift
        x$pool_mean[, id] <- x$pool_mean[, id] + shift
        optimum[id] <- optimum[id] + shift
      }
    }
    # the mother keeps what is left; where nothing is left, its population is gone
    n[, c_j] <- n[, c_j] - as.integer(taken)
    emptied <- which(taken > 0 & n[, c_j] == 0)
    for (nm in c("arrival", "propagule", "z", "ph")) x[[nm]][emptied, j] <- NA
    for (nm in c("immigrants", "load")) x[[nm]][emptied, j] <- 0
    x$present[emptied, j] <- FALSE
  }
  for (nm in fields) assign(nm, x[[nm]], envir = st)
  st$ev$optimum <- optimum
  st$n <- n
  new_ids <- unlist(lapply(buds, function(b) b$ids))
  st$n <- cbind(st$n, new_n)
  # emigrants of this step stay aligned with the alive species (none for the new ones)
  if (!is.null(st$em)) st$em <- cbind(st$em, matrix(0, m, q))
  st$alive <- c(st$alive, new_ids)
  st$k_all <- st$k_all + q
  st$k <- length(st$alive)
  st$p_em <- rep(st$dispersal[st$alive], each = m)
  st$p_ex <- rep(st$extirpation[st$alive], each = m)
  st$at <- matrix(0, st$k, st$k)
  invisible(st)
}

# the full tree (with extinct species) from the table of species: each species is a lineage
# from its birth to its end, from which the species that budded from it branch off at their
# birth, in their order; species at the start hang from the tips of `tree0` (or from the root,
# alone). Built from the youngest lineages to the oldest, without recursion, since a lineage
# may have many buds.
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

# rows of the series for generation `gen`: for each lineage alive at its start, its occupancy,
# abundance and mean genetic optimum (weighted by abundance) at its end, and the sites it
# colonised and lost since the last generation, out of those it did not and did occupy then.
# Lineages that died out in this generation get a row with their last losses. A site a mother
# passed to its buds in this generation counts neither as a loss nor as exposure; a lineage
# born in this generation has no transitions yet.
evolving_series_rows <- function(st, before, buds, sp_table, gen) {
  dead <- setdiff(which(colSums(before) > 0), st$alive)
  ids <- c(st$alive, dead)
  n <- cbind(st$n, matrix(0, st$m, length(dead)))
  now <- n > 0
  was <- before[, ids, drop = FALSE]
  skip <- matrix(FALSE, st$m, length(ids))
  for (b in buds) {
    c_j <- match(b$j, ids)
    if (is.na(c_j)) next
    given <- unique(unlist(b$sites))
    given <- given[!now[given, c_j]]
    skip[given, c_j] <- TRUE
  }
  newborn <- sp_table$birth[ids] == gen & gen > 0
  zz <- cbind(st$z[, st$alive, drop = FALSE], matrix(0, st$m, length(dead)))
  zz[!now] <- 0
  abundance <- colSums(n)
  keep <- !skip
  out <- data.frame(generation = rep(gen, length(ids)), lineage = st$species[ids],
                    occupancy = colSums(now), abundance = abundance,
                    genetic = ifelse(abundance > 0, colSums(n * zz) / abundance, NA_real_),
                    colonised = colSums(now & !was & keep),
                    extirpated = colSums(!now & was & keep),
                    present_before = colSums(was & keep), absent_before = colSums(!was & keep),
                    stringsAsFactors = FALSE)
  out[newborn, c("colonised", "extirpated", "present_before", "absent_before")] <- NA
  out
}

# extant lineages grouped into good species: a lineage that completed by its end (`status`
# "good") is its own species; an incipient one belongs to the species its parent belonged to at
# its birth, as in the labels of the PBD; `group` maps lineages to species and
# `rep` names, for each species (in the order of the table), the lineage that represents it:
# itself if extant, otherwise its earliest extant member
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

# a record of the extant lineages (sites x lineages) by good species: summed, or taken from the
# most abundant lineage of the species in each site
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

# names for q new species that no species has used
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

# drops from the state the alive species with no individual left, as mc_step() does at the
# start of the next step; returns their indices over all species ever
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
  st$k <- length(st$alive)
  st$p_em <- rep(st$dispersal[st$alive], each = st$m)
  st$p_ex <- rep(st$extirpation[st$alive], each = st$m)
  lost
}
