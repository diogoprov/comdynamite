#' History of the geographic ranges of lineages, as a stochastic map
#'
#' Maps, along the tree of [simulate_evolving_metacommunity()], the range of each lineage at
#' each moment, as the combination of areas in which it occupies at least one site: `"A"`,
#' `"B"`, `"A+B"` and so on, with the areas in the order of their levels. The result is a tree
#' of class `simmap`, as from [trait_mk()] or phytools, with the true history, for
#' [history_rates()], [as_herodotools()], [arrival_ages()] and [region_pool()], and for
#' comparison with ranges reconstructed by models such as DEC (BioGeoBEARS).
#'
#' The range of a lineage changes at the end of each generation in which it colonises the
#' first site, or loses the last site, of an area. A new lineage starts with the sites it took
#' from its mother. Branches of the tree follow lineages: a branch that continues a mother after
#' a budding is the mother, and the stem of a new lineage is the new lineage. In the
#' reconstructed tree a branch may join several lineages of the full tree (an extinct mother and
#' the descendant that survived), and its map joins their ranges in order. Each tip of the
#' reconstructed tree has the range of the lineage that represents its species (with protracted
#' speciation, not the union of the incipient lineages of the species). Time before the start
#' of the simulation (branches of a starting tree, or a stem) takes the range at the start.
#'
#' @param x An `mc_evolving` object from [simulate_evolving_metacommunity()].
#' @param areas Area of each site, in the order of the sites or named by site number; a factor
#'   gives the order of the areas in the names of the ranges. Area names may not contain `"+"`.
#' @param tree `"reconstructed"` (the default) for the tree of extant species, `x$tree`, or
#'   `"full"` for the tree of all lineages, extinct ones included, `x$tree_full`.
#' @return A tree of class `c("simmap", "phylo")` with `maps` (for each edge, the time spent in
#'   each range, from its start) and `mapped.edge` (edges x ranges, the total time), or `NULL`
#'   when the tree is `NULL`.
#' @seealso [history_rates()], [as_herodotools()], [arrival_ages()]
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' ls <- landscape_grid(5, disp_scale = 0.5)
#' areas <- ifelse(ls$coords$x < 2.5, "W", "E")
#' ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
#' sim <- simulate_evolving_metacommunity(ls, ancestor, dispersal = 0.01, extirpation = 0.02,
#'                                        n_steps = 200,
#'                                        speciation = speciation_control("fission", rate = 0.02))
#' h <- range_history(sim, areas)
#' colSums(h$mapped.edge)
#' @export
range_history <- function(x, areas, tree = c("reconstructed", "full")) {
  if (!inherits(x, "mc_evolving")) {
    stop("`x` must come from simulate_evolving_metacommunity()", call. = FALSE)
  }
  if (is.null(x$ranges)) {
    stop("`x` has no `ranges`: rerun it with this version of the package", call. = FALSE)
  }
  tree <- match.arg(tree)
  m <- nrow(x$abundance)
  if (!is.null(names(areas)) && !is.factor(areas)) {
    at <- match(as.character(seq_len(m)), names(areas))
    if (anyNA(at)) at <- match(paste0("site", seq_len(m)), names(areas))
    if (anyNA(at)) stop("`areas` names must be the site numbers", call. = FALSE)
    areas <- areas[at]
  }
  if (length(areas) != m || anyNA(areas)) {
    stop("`areas` needs one area per site (", m, ")", call. = FALSE)
  }
  areas <- if (is.factor(areas)) droplevels(areas) else factor(areas, unique(areas))
  if (any(grepl("+", levels(areas), fixed = TRUE))) {
    stop("area names may not contain \"+\"", call. = FALSE)
  }
  full <- x$tree_full
  timeline <- range_timelines(x$ranges, areas)
  full_maps <- range_full_maps(full, x$species, timeline)
  if (tree == "full") return(ranges_simmap(full, full_maps))
  if (is.null(x$tree)) return(NULL)
  rec <- x$tree
  rep <- range_representatives(x)
  rec_maps <- range_project_maps(full, full_maps, rec, rep[rec$tip.label])
  ranges_simmap(rec, rec_maps)
}

# for each lineage, the times (generations) at which its range changed and the range from then
range_timelines <- function(ranges, areas) {
  levels <- levels(areas)
  ranges <- ranges[order(ranges$generation), , drop = FALSE]
  by_lineage <- split(ranges, factor(ranges$lineage, unique(ranges$lineage)))
  lapply(by_lineage, function(r) {
    count <- stats::setNames(integer(length(levels)), levels)
    gens <- unique(r$generation)
    state <- character(length(gens))
    for (i in seq_along(gens)) {
      ri <- r[r$generation == gens[i], , drop = FALSE]
      delta <- tapply(ifelse(ri$present, 1L, -1L), factor(areas[ri$site], levels), sum)
      delta[is.na(delta)] <- 0L
      count <- count + as.integer(delta)
      state[i] <- paste(levels[count > 0], collapse = "+")
    }
    keep <- c(TRUE, state[-1] != state[-length(state)])
    data.frame(time = gens[keep], state = state[keep], stringsAsFactors = FALSE)
  })
}

# the ranges of a lineage over [from, to] as a named vector of durations; time before its first
# record takes its first range
range_segment <- function(tl, from, to) {
  # changes at the end of the segment count, with no duration, so that its end has the range
  # reached then; an empty range (at extinction) is left out
  change <- tl$time[tl$time > from & tl$time <= to]
  starts <- c(from, change)
  ends <- c(change, to)
  state <- vapply(starts, function(t) {
    i <- which(tl$time <= t)
    tl$state[if (length(i)) max(i) else 1]
  }, character(1))
  d <- stats::setNames(ends - starts, state)
  d <- d[names(d) != "" & (d > 0 | seq_along(d) == length(d))]
  if (!length(d)) d <- stats::setNames(to - from, tl$state[tl$state != ""][1])
  # merge consecutive pieces in the same range
  run <- cumsum(c(TRUE, names(d)[-1] != names(d)[-length(d)]))
  stats::setNames(as.vector(tapply(d, run, sum)), names(d)[!duplicated(run)])
}

# maps of the edges of the full tree: each edge follows the lineage that has, among the tips it
# leads to, the earliest birth (ties by the order of budding, then of the table)
range_full_maps <- function(full, sp, timeline) {
  n_tip <- ape::Ntip(full)
  depth <- ape::node.depth.edgelength(full)
  tip_sp <- match(full$tip.label, sp$species)
  offset <- stats::median(depth[seq_len(n_tip)] - sp$end[tip_sp])
  key <- order(order(sp$birth[tip_sp], sp$order[tip_sp], tip_sp))
  # earliest tip below each node
  first <- integer(n_tip + full$Nnode)
  first[seq_len(n_tip)] <- seq_len(n_tip)
  first[-seq_len(n_tip)] <- NA_integer_
  for (e in ape::reorder.phylo(full, "postorder", index.only = TRUE)) {
    p <- full$edge[e, 1]
    c <- full$edge[e, 2]
    if (is.na(first[p]) || key[first[c]] < key[first[p]]) first[p] <- first[c]
  }
  lapply(seq_len(nrow(full$edge)), function(e) {
    lineage <- full$tip.label[first[full$edge[e, 2]]]
    from <- depth[full$edge[e, 1]] - offset
    to <- depth[full$edge[e, 2]] - offset
    # a lineage never present at the end of a generation (born and lost in the same one)
    # takes the range of its mother
    tl <- timeline[[lineage]] %||% timeline[[sp$parent[match(lineage, sp$species)]]]
    range_segment(tl, from, to)
  })
}

# the lineage of the full tree that represents each species of the reconstructed tree
range_representatives <- function(x) {
  obs <- evolving_observed(x$species)
  obs$rep
}

# maps of the reconstructed tree from those of the full tree: each of its edges joins the
# edges of the full tree between the nodes that are the ancestors of the same tips
range_project_maps <- function(full, full_maps, rec, tips) {
  n_rec <- ape::Ntip(rec)
  tip_full <- match(tips, full$tip.label)
  # tips of the full tree below each node of the reconstructed one, and their ancestor there
  below <- lapply(seq_len(n_rec + rec$Nnode), function(v) {
    if (v <= n_rec) return(v)
    range_tips_below(rec, v)
  })
  node_full <- vapply(below, function(t) {
    if (length(t) == 1) return(tip_full[t])
    as.integer(ape::getMRCA(full, tip_full[t]))
  }, integer(1))
  edge_to <- match(seq_len(ape::Ntip(full) + full$Nnode), full$edge[, 2])
  lapply(seq_len(nrow(rec$edge)), function(e) {
    top <- node_full[rec$edge[e, 1]]
    v <- node_full[rec$edge[e, 2]]
    path <- integer(0)
    while (v != top) {
      path <- c(edge_to[v], path)
      v <- full$edge[edge_to[v], 1]
    }
    d <- unlist(unname(full_maps[path]))
    run <- cumsum(c(TRUE, names(d)[-1] != names(d)[-length(d)]))
    stats::setNames(as.vector(tapply(d, run, sum)), names(d)[!duplicated(run)])
  })
}

# tips below an internal node
range_tips_below <- function(tree, node) {
  n_tip <- ape::Ntip(tree)
  out <- integer(0)
  todo <- node
  while (length(todo)) {
    kids <- tree$edge[tree$edge[, 1] %in% todo, 2]
    out <- c(out, kids[kids <= n_tip])
    todo <- kids[kids > n_tip]
  }
  out
}

ranges_simmap <- function(tree, maps) {
  states <- unique(unlist(lapply(maps, names), use.names = FALSE))
  mapped <- matrix(0, length(maps), length(states), dimnames = list(NULL, states))
  for (e in seq_along(maps)) {
    total <- tapply(maps[[e]], names(maps[[e]]), sum)
    mapped[e, names(total)] <- total
  }
  rownames(mapped) <- paste(tree$edge[, 1], tree$edge[, 2], sep = ",")
  tree$maps <- maps
  tree$mapped.edge <- mapped
  class(tree) <- c("simmap", "phylo")
  tree
}
