#' Simulations as input for Herodotools, with the known history
#'
#' Herodotools (Nakamura et al. 2024) estimates the biogeographic history of assemblages from a
#' tree, the occurrence of species in sites, the region of each site and the area of each
#' ancestor. `as_herodotools()` gives these inputs from a simulation or an `mc_data` object, with
#' the names of the arguments of `Herodotools::calc_age_arrival()`,
#' `Herodotools::calc_dispersal_from()` and `Herodotools::calc_insitu_diversification()`, so the
#' result can go to them through [do.call()]. When the areas evolved as a discrete trait from
#' [trait_mk()], whose states are the regions, `history` gives the true area of each ancestor;
#' otherwise reconstruct them elsewhere, for example with BioGeoBEARS. [region_pool()] limits
#' each species to the sites of its region in the simulation, and [arrival_ages()] gives the
#' true ages at which lineages arrived in each region, to compare with what Herodotools
#' estimates from the nodes.
#'
#' For an evolving metacommunity from [simulate_evolving_metacommunity()], `history` defaults to
#' the true history of the ranges of its lineages, from [range_history()], whose ranges join
#' regions with `"+"` (`"A+B"`); Herodotools then counts an ancestor in `"A+B"` as present in
#' both regions.
#'
#' Herodotools finds a region among the areas of the ancestors by matching its name as a
#' pattern, so no region name may be contained in another, such as `"A"` and `"AB"`. Sites
#' without species are dropped with a warning, since Herodotools needs at least one species in
#' each site.
#'
#' @param x An `mc_sim`, `mc_evolving` or `mc_data` object, or a community matrix with sites in
#'   rows and species in columns, named by species.
#' @param regions Region of each site, by site name or in the order of the sites.
#' @param history Optional `simmap` tree from [trait_mk()] or [range_history()] whose states are
#'   the regions, or ranges of several regions joined by `"+"`, with the same species as `x`.
#' @param tree The tree, when `x` has none and `history` is not given.
#' @return A list with `W`, the incidence matrix of sites x species; `tree`, a plain `phylo`
#'   object; `biogeo`, a one-column data frame with the region of each site; and, with
#'   `history`, `ancestral.area`, a one-column data frame with the area of each internal node,
#'   in the order of their numbers in `tree$edge` (the root first), as Herodotools reads it.
#' @references Nakamura, G., Rodrigues, A. V., Luza, A. L., Maestri, R., Debastiani, V. J., &
#'   Duarte, L. D. S. (2024). Herodotools: An R package to integrate macroevolution,
#'   biogeography and community ecology. *Journal of Biogeography*, 51, 533--543.
#'   \doi{10.1111/jbi.14610}
#' @seealso [trait_mk()], [region_pool()], [arrival_ages()], [tip_history_metrics()]
#' @family species traits
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(15)
#' areas <- trait_mk(tree, q = 0.5, states = c("A", "B"))
#' regions <- rep(c("A", "B"), each = 10)
#' sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 200,
#'                               pool = region_pool(areas, regions))
#' input <- as_herodotools(sim, regions, history = areas)
#' str(input, max.level = 1)
#' arrival_ages(areas, sim, regions)[c(1, 11), ]
#' # with Herodotools installed: do.call(Herodotools::calc_age_arrival, input)
#' @export
as_herodotools <- function(x, regions, history = NULL, tree = NULL) {
  if (inherits(x, "mc_evolving")) {
    if (is.null(history)) history <- range_history(x, regions)
    x <- as_metacommunity(x)
  }
  d <- if (inherits(x, c("mc_sim", "mc_data"))) as_mc_data(x) else mc_data(x)
  tree <- history %||% tree %||% d$tree
  if (is.null(tree)) stop("give `history` or `tree`: `x` has no tree", call. = FALSE)
  if (!setequal(tree$tip.label, colnames(d$comm))) {
    stop("the tree and `x` must have the same species", call. = FALSE)
  }
  regions <- as.character(site_rows(regions, rownames(d$comm), "`regions`"))
  check_region_names(unique(regions))
  W <- (d$comm[, tree$tip.label, drop = FALSE] > 0) * 1
  empty <- rowSums(W) == 0
  if (any(empty)) {
    warning("sites without species dropped: ", species_list(rownames(W)[empty]), call. = FALSE)
  }
  plain <- tree[intersect(c("edge", "edge.length", "tip.label", "Nnode", "node.label"),
                          names(tree))]
  class(plain) <- "phylo"
  out <- list(W = W[!empty, , drop = FALSE], tree = plain,
              biogeo = data.frame(region = regions[!empty], row.names = rownames(W)[!empty]))
  if (!is.null(history)) {
    areas <- history_states(history)$nodes
    parts <- unique(unlist(strsplit(areas, "+", fixed = TRUE)))
    check_region_names(unique(c(regions, parts)))
    out$ancestral.area <- data.frame(area = areas,
                                     row.names = paste0("N", ape::Ntip(history) + seq_along(areas)))
  }
  out
}

#' Sites where each species can occur, from its region
#'
#' Gives the `pool` of [simulate_metacommunity()] from a history of areas: each species can occur
#' only in the sites of the region (or regions, in a range such as `"A+B"`) it occupies at the
#' tips of `history`.
#'
#' @param history A `simmap` tree from [trait_mk()] whose states are the regions.
#' @param regions Region of each site. Names, when given, name the sites; otherwise they are
#'   `site1`, `site2` and so on, as in [simulate_metacommunity()].
#' @return Logical matrix of sites x species.
#' @inherit as_herodotools examples
#' @family species traits
#' @export
region_pool <- function(history, regions) {
  tips <- history_states(history)$tips
  sites <- names(regions) %||% paste0("site", seq_along(regions))
  pool <- outer(as.character(regions), tips, Vectorize(range_has))
  dimnames(pool) <- list(sites, names(tips))
  pool
}

#' True ages at which lineages arrived in each region
#'
#' For each species in each site, the time before the present at which its lineage last entered
#' the region of the site, from the history of areas: the time spent in that region along the
#' path from the tip towards the root, until the first period in another region, plus the time
#' from the tip to the present when the tree is not ultrametric. A lineage that stayed in the
#' region since the root gets the age of the root. With ranges of several regions (`"A+B"`, from
#' [range_history()]), the age in each region is the time since the range last came to include
#' it. These are the true values of what
#' `Herodotools::calc_age_arrival()` estimates from the areas of the nodes, which can only see
#' arrivals at nodes and so gives younger ages.
#'
#' @inheritParams as_herodotools
#' @param history A `simmap` tree from [trait_mk()] whose states are the regions.
#' @return Matrix of sites x species with the age of arrival, `NA` where the species is absent
#'   or where its tip is in another region.
#' @inherit as_herodotools examples
#' @family species traits
#' @export
arrival_ages <- function(history, x, regions) {
  d <- if (inherits(x, c("mc_sim", "mc_data"))) as_mc_data(x) else mc_data(x)
  regions <- as.character(site_rows(regions, rownames(d$comm), "`regions`"))
  edge <- history$edge
  depth <- ape::node.depth.edgelength(history)
  height <- max(depth)
  edge_to <- match(seq_len(ape::Ntip(history) + history$Nnode), edge[, 2])
  since <- lapply(stats::setNames(seq_len(ape::Ntip(history)), history$tip.label), function(i) {
    edges <- integer(0)
    node <- i
    while (!is.na(edge_to[node])) {
      edges <- c(edges, edge_to[node])
      node <- edge[edge_to[node], 1]
    }
    segments <- unlist(lapply(edges, function(e) rev(history$maps[[e]])))
    # for each region in the range at the tip, the time since the lineage last entered it
    here <- strsplit(names(segments)[1], "+", fixed = TRUE)[[1]]
    vapply(stats::setNames(here, here), function(r) {
      away <- which(!range_has(r, names(segments)))[1]
      age <- if (is.na(away)) sum(segments) else sum(segments[seq_len(away - 1)])
      age + height - depth[i]
    }, numeric(1))
  })
  species <- colnames(d$comm)
  out <- matrix(NA_real_, nrow(d$comm), length(species), dimnames = dimnames(d$comm))
  for (j in species) {
    present <- d$comm[, j] > 0 & regions %in% names(since[[j]])
    out[present, j] <- since[[j]][regions[present]]
  }
  out
}

# does each range (a state such as "A" or "A+B") include the region?
range_has <- function(region, range) {
  vapply(strsplit(range, "+", fixed = TRUE), function(p) region %in% p, logical(1))
}

# states at the internal nodes (root first) and at the tips of a simmap tree
history_states <- function(history) {
  if (!inherits(history, "phylo") || is.null(history$maps)) {
    stop("`history` must be a `simmap` tree, such as from `trait_mk()`", call. = FALSE)
  }
  n_tip <- ape::Ntip(history)
  states <- character(n_tip + history$Nnode)
  states[history$edge[, 1]] <- vapply(history$maps, function(m) names(m)[1], character(1))
  ends <- vapply(history$maps, function(m) names(m)[length(m)], character(1))
  states[history$edge[, 2]] <- ends
  list(nodes = states[n_tip + seq_len(history$Nnode)],
       tips = stats::setNames(states[seq_len(n_tip)], history$tip.label))
}

check_region_names <- function(names) {
  inside <- vapply(names, function(a) sum(grepl(a, names, fixed = TRUE)) > 1, logical(1))
  if (any(inside)) {
    stop("region names may not be contained in one another: ", species_list(names[inside]),
         call. = FALSE)
  }
}
