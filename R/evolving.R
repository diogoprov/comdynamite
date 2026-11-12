#' Speciation settings
#'
#' Settings for speciation in [simulate_evolving_metacommunity()]. Speciation is checked once
#' per generation (`steps_per_generation` of [evolution_control()]) among the populations of
#' each species. When they fall into two or more groups, the species buds: the largest group
#' (most sites, then most individuals) keeps the identity and age of the species, and each other
#' group becomes a new species with its populations, abundances and genetic means. Several new
#' species at once branch off one after another, the most divergent first, at the same time
#' (internal branches of length zero), so that the tree stays binary.
#'
#' Gene flow between two populations \eqn{a}{a} and \eqn{b}{b} is measured by the migration
#' rate, the larger of the expected migrants from one to the other in the generation over the
#' size of the receiving population, so that it does not depend on abundance.
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
#'
#' A new population takes the divergences (or counts) of the population that sent it most
#' migrants, and a population lost at any step of a generation founds anew.
#'
#' @param mode `"geographic"` or `"ecological"`.
#' @param threshold Divergence in generations (`"geographic"`) or difference of genetic means in
#'   units of the environment (`"ecological"`) that separates groups; `Inf` turns speciation
#'   off.
#' @param divergence For `"geographic"`: `"cluster"` or `"gradual"`.
#' @param migration_min For `divergence = "cluster"`: migration rate per generation that links
#'   two populations.
#' @param reset_generations For `divergence = "cluster"`: generations in a row in the same
#'   cluster after which divergence returns to 0.
#' @param homogenization For `divergence = "gradual"`: reduction of divergence per migrant.
#' @param min_generations For `"ecological"`: generations in a row of difference and reduced
#'   gene flow before groups split.
#' @param migration_max For `"ecological"`: migration rate per generation below which gene flow
#'   counts as reduced.
#' @return An `mc_speciation` object.
#' @references Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M., Rangel,
#'   T. F., & Pellissier, L. (2021). gen3sis: A general engine for eco-evolutionary simulations
#'   of the processes that shape Earth's biodiversity. *PLOS Biology*, 19, e3001340.
#'   \doi{10.1371/journal.pbio.3001340}
#' @family simulate metacommunities
#' @examples
#' speciation_control("geographic", threshold = 50)
#' speciation_control("geographic", threshold = 50, divergence = "gradual", homogenization = 0.5)
#' speciation_control("ecological", threshold = 20, min_generations = 30)
#' @export
speciation_control <- function(mode = c("geographic", "ecological"), threshold,
                               divergence = c("cluster", "gradual"), migration_min = 0.001,
                               reset_generations = 10, homogenization = 1,
                               min_generations = 20, migration_max = 0.001) {
  mode <- match.arg(mode)
  divergence <- match.arg(divergence)
  if (missing(threshold)) stop("`threshold` is needed", call. = FALSE)
  one <- function(x) is.numeric(x) && length(x) == 1
  stopifnot(one(threshold), threshold > 0, one(migration_min), migration_min > 0,
            one(reset_generations), reset_generations >= 1, one(homogenization),
            homogenization >= 0, one(min_generations), min_generations >= 1,
            one(migration_max), migration_max >= 0)
  structure(list(mode = mode, threshold = threshold, divergence = divergence,
                 migration_min = migration_min, reset_generations = reset_generations,
                 homogenization = homogenization, min_generations = min_generations,
                 migration_max = migration_max),
            class = "mc_speciation")
}

#' @export
print.mc_speciation <- function(x, ...) {
  cat("<mc_speciation>", x$mode, "| threshold =", format(x$threshold))
  if (x$mode == "geographic") {
    cat(" | divergence:", x$divergence,
        if (x$divergence == "cluster") {
          paste("| migration_min =", format(x$migration_min), "| reset_generations =",
                format(x$reset_generations))
        } else {
          paste("| homogenization =", format(x$homogenization))
        })
  } else {
    cat(" | min_generations =", format(x$min_generations), "| migration_max =",
        format(x$migration_max))
  }
  cat("\n")
  invisible(x)
}

#' Simulate a metacommunity with speciation
#'
#' Runs the model of [simulate_metacommunity()], with the evolution of the niche optimum within
#' populations of [evolution_control()], on species that bud as their populations diverge
#' ([speciation_control()]) and die out when their last population is lost. The tree is not
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
#'   species fixed (`genetic_var = 0`), which allows only geographic speciation.
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
#' @return An `mc_evolving` object: a list with `abundance` (sites x extant species), `tree`
#'   (the reconstructed tree of the extant species, branch lengths in generations, or `NULL`
#'   with fewer than two), `tree_full` (with extinct species), `species` (a table with the
#'   parent, birth, end and fate of every species), `events` (one row per budding, with the new
#'   species), `richness` (extant species at each generation), `genetic` (genetic means of the
#'   populations of the extant species), `traits` (of the extant species, with the mean genetic
#'   optimum over their populations as `optimum`), `alpha` (among the extant species at the end,
#'   without the niche overlap, which varies among sites), the records and colonisation history
#'   of [simulate_metacommunity()] for the extant species, and the `landscape`, `dispersal`,
#'   `extirpation`, `speciation`, `evolution` and `n_steps` used. [as_metacommunity()] turns it
#'   into an `mc_sim` for the other functions of the package.
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
  if (is.null(evolution)) evolution <- evolution_control(0, ne_ratio = Inf)
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
                         end = NA_real_, fate = "extant", order = 0L,
                         stringsAsFactors = FALSE)
  split_time <- split0
  events <- data.frame(time = numeric(0), parent = character(0), new = character(0),
                       stringsAsFactors = FALSE)
  # per species (indexed over all species ever): divergence or count among its populations,
  # generations in a row in one cluster, and expected migrants (from row to column) in the
  # current generation
  empty <- function() matrix(NA_real_, st$m, st$m)
  divergence <- lapply(seq_len(k0), function(j) empty())
  together <- lapply(seq_len(k0), function(j) matrix(0, st$m, st$m))
  flow <- lapply(seq_len(k0), function(j) matrix(0, st$m, st$m))
  richness <- integer(0)
  for (step in seq_len(n_steps)) {
    mc_step(st, step)
    gen <- step / per
    # species lost at this step leave now (mc_step() would drop them at the next one, with no
    # random draws either way), so that the table, richness and output agree
    lost <- evolving_prune(st)
    if (length(lost)) {
      sp_table$end[lost] <- gen
      sp_table$fate[lost] <- "extinct"
    }
    if (st$k > 0) {
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
      next
    }
    upd <- evolving_divergence(st, divergence, together, flow, speciation)
    divergence <- upd$divergence
    together <- upd$together
    for (j in st$alive) {
      groups <- evolving_groups(st, j, divergence[[j]], speciation)
      if (length(groups) < 2) next
      groups <- evolving_order(st, j, groups, divergence[[j]], speciation)
      q <- length(groups) - 1
      if (st$k_all + q > max_species) {
        stop("more than `max_species` (", max_species, ") species; raise it or the ",
             "threshold of speciation", call. = FALSE)
      }
      names_new <- evolving_names(st$species, q)
      new_ids <- st$k_all + seq_len(q)
      evolving_split(st, j, groups, names_new)
      sp_table <- rbind(sp_table, data.frame(species = names_new,
                                             parent = sp_table$species[j], birth = gen,
                                             end = NA_real_, fate = "extant",
                                             order = seq_len(q), stringsAsFactors = FALSE))
      events <- rbind(events, data.frame(time = gen, parent = sp_table$species[j],
                                         new = paste(names_new, collapse = ","),
                                         stringsAsFactors = FALSE))
      split_time <- evolving_split_time(split_time, j, q, gen, names_new)
      d <- divergence[[j]]
      tg <- together[[j]]
      for (g in seq_along(groups)) {
        id <- if (g == 1) j else new_ids[g - 1]
        out <- !(seq_len(st$m) %in% groups[[g]])
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
    flow <- lapply(flow, function(x) x * 0)
    if (competition == "phylo" || nrow(events)) {
      a <- evolving_alpha(split_time[st$alive, st$alive, drop = FALSE], gen, competition,
                          inter_max, inter_min, scale, distance_scale)
      st$at <- t(a)
    }
    richness <- c(richness, length(st$alive))
  }
  out <- mc_finish(st)
  extant <- sp_table$species[sp_table$fate == "extant"]
  sp_table$end[sp_table$fate == "extant"] <- generations
  tree_full <- evolving_tree(sp_table, tree0)
  tree <- if (length(extant) >= 2) ape::keep.tip(tree_full, extant) else NULL
  ext_idx <- match(extant, st$species)
  genetic <- out$evolution$genetic[, extant, drop = FALSE]
  traits_out <- data.frame(optimum = colMeans(genetic, na.rm = TRUE),
                           sigma = st$ev$sigma[ext_idx], r_max = st$ev$r_max[ext_idx],
                           row.names = extant)
  alpha <- evolving_alpha(split_time[extant, extant, drop = FALSE], generations, competition,
                          inter_max, inter_min, scale, distance_scale)
  structure(list(abundance = out$abundance[, extant, drop = FALSE], tree = tree,
                 tree_full = tree_full, species = sp_table, events = events,
                 richness = richness, genetic = genetic, traits = traits_out, alpha = alpha,
                 landscape = landscape, dispersal = st$dispersal[ext_idx],
                 extirpation = st$extirpation[ext_idx], speciation = speciation,
                 evolution = out$evolution$settings, n_steps = n_steps,
                 generations = generations,
                 history = if (!is.null(out$history)) out$history[, extant, , drop = FALSE],
                 arrival = out$arrival[, extant, drop = FALSE],
                 propagule = out$propagule[, extant, drop = FALSE],
                 immigrants = out$immigrants[, extant, drop = FALSE],
                 migrant_load = out$migrant_load[, extant, drop = FALSE]),
            class = "mc_evolving")
}

#' @export
print.mc_evolving <- function(x, ...) {
  fates <- table(factor(x$species$fate, c("extant", "extinct")))
  cat("<mc_evolving>", nrow(x$abundance), "sites,", x$generations, "generations\n")
  cat("  species ever:", nrow(x$species), "| extant:", fates[["extant"]], "| extinct:",
      fates[["extinct"]], "| buddings:", nrow(x$events), "\n")
  cat("  speciation:", x$speciation$mode, "with threshold", format(x$speciation$threshold),
      "\n")
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
                 n_steps = x$n_steps), class = "mc_sim")
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
    width <- 2 * outer(sigma2, sigma2, "+")
    out <- matrix(0, nrow(n), ncol(n))
    for (s in which(rowSums(n) > 0)) {
      overlap <- exp(-outer(p[s, ], p[s, ], "-")^2 / width)
      out[s, ] <- n[s, ] %*% (st$at * overlap)
    }
    out
  }
}

# split times after q species bud from species j at generation `gen`: the new species share
# the history of j up to `gen`
evolving_split_time <- function(split_time, j, q, gen, names_new) {
  k <- nrow(split_time)
  out <- matrix(0, k + q, k + q)
  out[seq_len(k), seq_len(k)] <- split_time
  new <- k + seq_len(q)
  out[new, seq_len(k)] <- matrix(split_time[j, ], q, k, byrow = TRUE)
  out[new, j] <- gen
  out[seq_len(k), new] <- t(out[new, seq_len(k)])
  out[new, new] <- gen
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
      x <- d[s, groups[[keep]], drop = FALSE]
      mean(x, na.rm = TRUE)
    } else {
      abs(mean(st$z[s, j]) - mean(st$z[groups[[keep]], j]))
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

# species j (index over all species ever) buds: it keeps the sites of the first group, and each
# other group becomes a new species, for which every per-species structure of the state grows;
# new species take the abundances, genetic means and history of the mother in their sites
evolving_split <- function(st, j, groups, names_new) {
  q <- length(groups) - 1
  m <- st$m
  grow_cols <- function(x, fill) cbind(x, matrix(fill, nrow(x), q))
  grow_vec <- function(x) c(x, rep(x[j], q))
  grow_hist <- function(x, fill) {
    if (is.null(x)) return(NULL)
    d <- dim(x)
    out <- array(fill, c(d[1], d[2] + q, d[3]))
    out[, seq_len(d[2]), ] <- x
    dimnames(out) <- list(dimnames(x)[[1]], c(dimnames(x)[[2]], names_new), dimnames(x)[[3]])
    out
  }
  c_j <- match(j, st$alive)
  new_ids <- st$k_all + seq_len(q)
  for (nm in c("arrival", "propagule")) assign(nm, grow_cols(get(nm, st), NA_integer_), st)
  for (nm in c("immigrants", "load")) assign(nm, grow_cols(get(nm, st), 0), st)
  st$present <- grow_cols(st$present, FALSE)
  st$z <- grow_cols(st$z, NA_real_)
  st$ph <- grow_cols(st$ph, NA_real_)
  st$pool_mean <- grow_cols(st$pool_mean, NA_real_)
  st$pool_mean[, new_ids] <- st$pool_mean[, j]
  st$history <- grow_hist(st$history, 0L)
  st$z_history <- grow_hist(st$z_history, NA_real_)
  st$dispersal <- grow_vec(st$dispersal)
  st$extirpation <- grow_vec(st$extirpation)
  for (nm in c("g", "b", "ve", "optimum", "sigma", "r_max")) st$ev[[nm]] <- grow_vec(st$ev[[nm]])
  st$fixed_z <- grow_vec(st$fixed_z)
  st$species <- c(st$species, names_new)
  new_n <- matrix(0L, m, q)
  moved <- integer(0)
  for (g in seq_len(q)) {
    s <- groups[[g + 1]]
    id <- new_ids[g]
    moved <- c(moved, s)
    new_n[s, g] <- st$n[s, c_j]
    for (nm in c("arrival", "propagule", "immigrants", "load", "present", "z", "ph")) {
      x <- get(nm, st)
      x[s, id] <- x[s, j]
      assign(nm, x, st)
    }
  }
  # the mother keeps only the sites of the first group
  st$n[moved, c_j] <- 0L
  for (nm in c("arrival", "propagule", "z", "ph")) {
    x <- get(nm, st)
    x[moved, j] <- NA
    assign(nm, x, st)
  }
  for (nm in c("immigrants", "load")) {
    x <- get(nm, st)
    x[moved, j] <- 0
    assign(nm, x, st)
  }
  st$present[moved, j] <- FALSE
  st$n <- cbind(st$n, new_n)
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
# alone)
evolving_tree <- function(sp, tree0) {
  buds <- split(seq_len(nrow(sp)), sp$parent)
  buds <- lapply(buds, function(i) sp$species[i[order(sp$birth[i], sp$order[i])]])
  lineage <- function(s, start, k = 1, extra = 0) {
    i <- match(s, sp$species)
    b <- buds[[s]]
    if (is.null(b) || k > length(b)) {
      return(sprintf("%s:%.10g", s, sp$end[i] - start + extra))
    }
    t <- sp$birth[match(b[k], sp$species)]
    sprintf("(%s,%s):%.10g", lineage(b[k], t), lineage(s, t, k + 1), t - start + extra)
  }
  roots <- sp$species[is.na(sp$parent)]
  if (length(roots) == 1) {
    tr <- ape::read.tree(text = paste0("(", lineage(roots, 0), ");"))
    return(if (length(tr$tip.label) > 1) ape::collapse.singles(tr, root.edge = TRUE) else tr)
  }
  ntip <- length(tree0$tip.label)
  node_newick <- function(v, len) {
    if (v <= ntip) return(lineage(tree0$tip.label[v], 0, extra = len))
    kids <- tree0$edge[tree0$edge[, 1] == v, 2]
    lens <- tree0$edge.length[tree0$edge[, 1] == v]
    inner <- paste(mapply(node_newick, kids, lens), collapse = ",")
    sprintf("(%s):%.10g", inner, len)
  }
  text <- sub(":0$", "", node_newick(ntip + 1L, 0))
  ape::read.tree(text = paste0(text, ";"))
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
  st$k <- length(st$alive)
  st$p_em <- rep(st$dispersal[st$alive], each = st$m)
  st$p_ex <- rep(st$extirpation[st$alive], each = st$m)
  lost
}
