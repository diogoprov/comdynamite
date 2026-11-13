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
#' speciation, not the union of the incipient lineages of the species).
#'
#' A node has the range of the mother just before the split (the union of the ranges that the
#' mother and the new lineage start with, apart from changes in the same generation), and the
#' tree keeps the range at the root, which no edge shows, as `root_range`. The stem from the
#' start of the simulation to the first split (`root.edge`) is not mapped. Nodes of a starting
#' tree, before the simulation, take the union of the ranges at the start of the species below
#' them, and the parts of edges before the start take the range at the start.
#'
#' @param x An `mc_evolving` object from [simulate_evolving_metacommunity()].
#' @param areas Area of each site, in the order of the sites or named by site number; a factor
#'   gives the order of the areas in the names of the ranges. Area names may not contain `"+"`.
#' @param tree `"reconstructed"` (the default) for the tree of extant species, `x$tree`, or
#'   `"full"` for the tree of all lineages, extinct ones included, `x$tree_full`.
#' @return A tree of class `c("simmap", "phylo")` with `maps` (for each edge, the time spent in
#'   each range, from its start; a change in the last generation of a tip is kept as a piece
#'   without duration, so that the tip has the final range), `mapped.edge` (edges x ranges, the
#'   total time) and `root_range`, or `NULL` when the tree is `NULL`.
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
  if (tree == "reconstructed" && is.null(x$tree)) return(NULL)
  m <- nrow(x$abundance)
  if (!is.null(names(areas))) {
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
  mapped <- range_full_maps(full, x$species, timeline, levels(areas))
  if (tree == "full") return(ranges_simmap(full, mapped$maps, mapped$nodes[ape::Ntip(full) + 1]))
  rec <- x$tree
  rep <- range_representatives(x)
  projected <- range_project_maps(full, mapped$maps, rec, rep[rec$tip.label])
  ranges_simmap(rec, projected$maps, mapped$nodes[projected$root_full])
}

# for each lineage, the times (generations) at which its range changed and the range from then
range_timelines <- function(ranges, areas) {
  levels <- levels(areas)
  if (!nrow(ranges)) return(list())
  # one row per lineage and generation with changes, in order; columns are areas
  event <- paste(ranges$lineage, ranges$generation, sep = "\r")
  ord <- order(match(ranges$lineage, unique(ranges$lineage)), ranges$generation)
  key <- unique(event[ord])
  delta <- rowsum(outer(ifelse(ranges$present, 1, -1), rep(1, length(levels))) *
                    outer(as.integer(areas[ranges$site]), seq_along(levels), "=="),
                  factor(event, key), reorder = TRUE)
  lineage <- sub("\r.*$", "", key)
  time <- as.numeric(sub("^.*\r", "", key))
  # occupied sites of each area after each event, accumulated within each lineage
  count <- apply(delta, 2, function(d) stats::ave(d, lineage, FUN = cumsum))
  count <- matrix(count, nrow(delta))
  code <- as.vector((count > 0) %*% 2^(seq_along(levels) - 1))
  labels <- vapply(sort(unique(code)), function(cd) {
    paste(levels[bitwAnd(cd, 2^(seq_along(levels) - 1)) > 0], collapse = "+")
  }, character(1))
  state <- labels[match(code, sort(unique(code)))]
  keep <- !duplicated(lineage) | c(TRUE, state[-1] != state[-length(state)])
  split(data.frame(time = time[keep], state = state[keep], stringsAsFactors = FALSE),
        factor(lineage[keep], unique(lineage)))
}

range_eps <- 1e-6

# the ranges of a lineage over [from, to] as a named vector of durations, each range taken
# after the changes at its start; changes at `to` count, with no duration, only for terminal
# edges (the range at a tip), and an empty range (at extinction) is left out. Time before the
# first record takes the first range.
range_segment <- function(tl, from, to, terminal = FALSE) {
  eps <- range_eps
  change <- tl$time[tl$time > from + eps &
                      (tl$time < to - eps | (terminal & abs(tl$time - to) <= eps))]
  change <- pmin(change, to)
  starts <- c(from, change)
  ends <- c(change, to)
  state <- vapply(starts, function(t) {
    i <- which(tl$time <= t + eps)
    tl$state[if (length(i)) max(i) else 1]
  }, character(1))
  range_pieces(stats::setNames(ends - starts, state), tl$state[tl$state != ""][1])
}

# drops empty ranges and pieces without duration other than the last, and merges runs
range_pieces <- function(d, fallback = NULL) {
  d <- d[names(d) != "" & (d > 0 | seq_along(d) == length(d))]
  if (!length(d)) return(stats::setNames(0, fallback))
  run <- cumsum(c(TRUE, names(d)[-1] != names(d)[-length(d)]))
  stats::setNames(as.vector(tapply(d, run, sum)), names(d)[!duplicated(run)])
}

# the range of a lineage just before the changes at time t (at a node, before the split)
range_before <- function(tl, t) {
  i <- which(tl$time < t - range_eps)
  tl$state[if (length(i)) max(i) else 1]
}

range_union <- function(states, levels) {
  parts <- unlist(strsplit(states, "+", fixed = TRUE))
  paste(levels[levels %in% parts], collapse = "+")
}

# maps of the edges of the full tree, and the range at each node. Each edge follows the lineage
# that has, among the tips it leads to, the earliest birth (ties by the order of budding, then
# of the table); a node has the range of that lineage before the split. Nodes before the start
# (of a starting tree) take the union of the ranges at the start of the species below them.
range_full_maps <- function(full, sp, timeline, levels) {
  n_tip <- ape::Ntip(full)
  n_node <- n_tip + full$Nnode
  depth <- ape::node.depth.edgelength(full)
  tip_sp <- match(full$tip.label, sp$species)
  offset <- stats::median(depth[seq_len(n_tip)] - sp$end[tip_sp])
  # times in generations; they carry the rounding of the Newick text, so comparisons with the
  # generations of the log use a tolerance (`range_eps`)
  t <- depth - offset
  tl_of <- function(lineage) {
    timeline[[lineage]] %||% data.frame(time = 0, state = "", stringsAsFactors = FALSE)
  }
  key <- order(order(sp$birth[tip_sp], sp$order[tip_sp], tip_sp))
  first <- c(seq_len(n_tip), rep(NA_integer_, full$Nnode))
  post <- if (nrow(full$edge) > 1) {
    ape::reorder.phylo(full, "postorder", index.only = TRUE)
  } else {
    seq_len(nrow(full$edge))
  }
  # the start ranges of the species below each node, for nodes before the start
  below_start <- vector("list", n_node)
  for (i in which(sp$birth[tip_sp] == 0)) {
    below_start[[i]] <- range_before(tl_of(full$tip.label[i]), 2 * range_eps)
  }
  for (e in post) {
    p <- full$edge[e, 1]
    c <- full$edge[e, 2]
    if (is.na(first[p]) || key[first[c]] < key[first[p]]) first[p] <- first[c]
    below_start[p] <- list(c(below_start[[p]], below_start[[c]]))
  }
  node_state <- character(n_node)
  for (v in seq_len(n_node)[-seq_len(n_tip)]) {
    node_state[v] <- if (t[v] < -range_eps) range_union(below_start[[v]], levels) else
      range_before(tl_of(full$tip.label[first[v]]), t[v])
  }
  maps <- lapply(seq_len(nrow(full$edge)), function(e) {
    p <- full$edge[e, 1]
    c <- full$edge[e, 2]
    if (c > n_tip && t[c] < -range_eps) return(stats::setNames(t[c] - t[p], node_state[c]))
    range_segment(tl_of(full$tip.label[first[c]]), t[p], t[c], terminal = c <= n_tip)
  })
  list(maps = maps, nodes = node_state)
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
  node_full <- vapply(seq_len(n_rec + rec$Nnode), function(v) {
    if (v <= n_rec) return(tip_full[v])
    as.integer(ape::getMRCA(full, tip_full[range_tips_below(rec, v)]))
  }, integer(1))
  edge_to <- match(seq_len(ape::Ntip(full) + full$Nnode), full$edge[, 2])
  maps <- lapply(seq_len(nrow(rec$edge)), function(e) {
    top <- node_full[rec$edge[e, 1]]
    v <- node_full[rec$edge[e, 2]]
    path <- integer(0)
    while (v != top) {
      path <- c(edge_to[v], path)
      v <- full$edge[edge_to[v], 1]
    }
    range_pieces(unlist(unname(full_maps[path])))
  })
  list(maps = maps, root_full = node_full[n_rec + 1])
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

ranges_simmap <- function(tree, maps, root_range = NULL) {
  states <- unique(unlist(lapply(maps, names), use.names = FALSE))
  mapped <- matrix(0, length(maps), length(states), dimnames = list(NULL, states))
  for (e in seq_along(maps)) {
    total <- tapply(maps[[e]], names(maps[[e]]), sum)
    mapped[e, names(total)] <- total
  }
  rownames(mapped) <- paste(tree$edge[, 1], tree$edge[, 2], sep = ",")
  tree$maps <- maps
  tree$mapped.edge <- mapped
  # the range at the root, before its first split, which its edges cannot show
  tree$root_range <- root_range
  class(tree) <- c("simmap", "phylo")
  tree
}
