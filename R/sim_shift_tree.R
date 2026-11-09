# Internal: trees with heritable shifts in the rate of speciation, for the simulations of the
# DAMOCLES paper (data-raw/s1_type1_shifts.R and later blocks; dev/plano-artigo-damocles.md).
#
# A forward birth-death process with one rate per lineage. Speciation starts at `lambda0` and
# extinction is `epsilon` times speciation in every lineage. Each of `n_shifts` shifts
# multiplies the speciation (and so the extinction) rate of one living lineage and all its
# descendants by `ratio`; shift k happens when the number of living lineages first reaches a
# number drawn between 3 and `n_tips / 4`, on a lineage drawn at random, so shifts can be
# nested. The tree is conditioned on `n_tips` living lineages with the general sampling
# approach (GSA) of Hartmann, Wong & Stadler (2010, Syst. Biol. 59: 465-476,
# doi:10.1093/sysbio/syq026): each run goes on until `gsa_factor * n_tips` lineages are alive or
# all are extinct; a run contributes a tree with probability proportional to the total time Phi
# it spent with exactly `n_tips` lineages (their sampling rate r, here min(1, Phi / Phi_ref) with
# Phi_ref = gsa_scale / (n_tips * lambda0 * (1 + epsilon)), gsa_scale times the expected length
# of one such interval under the background rates); the tree is cut in one of those intervals,
# drawn with weight equal to its length, at a time uniform within it. Extinct lineages are
# pruned and the tree rescaled to height one. Trees in which a regime (the background included)
# holds fewer than `min_fraction` of the tips are drawn again (a further condition, not part of
# the GSA). The attribute "capped" of the result is TRUE when Phi exceeded Phi_ref, which
# truncates the weight; data-raw/validation_shift_tree.R reports how often.
#
# Returns a list with `tree` (a binary ultrametric phylo of height one), `regimes` (a factor,
# the regime of each edge in the order of tree$edge: "0" for the background and k for the
# lineages below shift k, the regime of the lower end of the edge), `tip_regime` (integer, per
# tip) and `rates` (speciation rate of each regime before rescaling).
sim_shift_tree <- function(n_tips, n_shifts = 1, ratio = 2, lambda0 = 1, epsilon = 0.5,
                           min_fraction = 0.05, gsa_factor = 2, gsa_scale = 15,
                           max_tries = 10000) {
  stopifnot(n_tips >= 8, n_shifts >= 0, ratio > 0, lambda0 > 0, epsilon >= 0, epsilon < 1,
            gsa_factor > 1)
  for (try in seq_len(max_tries)) {
    sim <- shift_tree_once(n_tips, n_shifts, ratio, lambda0, epsilon, gsa_factor, gsa_scale)
    if (is.null(sim)) next
    tips <- table(factor(sim$tip_regime, levels = seq_len(n_shifts + 1) - 1))
    if (all(tips >= min_fraction * n_tips)) return(sim)
  }
  stop("no tree met the conditions in ", max_tries, " tries", call. = FALSE)
}

shift_tree_once <- function(n_tips, n_shifts, ratio, lambda0, epsilon, gsa_factor, gsa_scale) {
  thresholds <- sort(sample(3:max(4, floor(n_tips / 4)), n_shifts, replace = TRUE))
  m <- ceiling(gsa_factor * n_tips)
  # Lineages alive: node they start from, start time, speciation rate, regime.
  from <- 0L
  start <- 0
  rate <- lambda0
  regime <- 0L
  alive <- 1L
  # Segments of lineages that ended: from, to (node ids), start, end, regime, kind of end
  # (1 speciation or shift, 2 extinction). Node ids are the row numbers.
  segments <- matrix(0, 64, 6)
  n_seg <- 0L
  close <- function(i, time, kind) {
    n_seg <<- n_seg + 1L
    if (n_seg > nrow(segments)) segments <<- rbind(segments, matrix(0, nrow(segments), 6))
    segments[n_seg, ] <<- c(from[i], n_seg, start[i], time, regime[i], kind)
    n_seg
  }
  rates <- lambda0
  t <- 0
  shifts <- 0L
  # Intervals with exactly n_tips lineages alive.
  n_start <- numeric(0)
  n_end <- numeric(0)
  while (length(alive) < m) {
    while (shifts < n_shifts && length(alive) >= thresholds[shifts + 1]) {
      shifts <- shifts + 1L
      i <- alive[sample.int(length(alive), 1)]
      from[i] <- close(i, t, 1)
      start[i] <- t
      rate[i] <- rate[i] * ratio
      regime[i] <- shifts
      rates[shifts + 1] <- rate[i]
    }
    dt <- stats::rexp(1, sum(rate[alive]) * (1 + epsilon))
    if (length(alive) == n_tips) {
      n_start <- c(n_start, t)
      n_end <- c(n_end, t + dt)
    }
    t <- t + dt
    i <- alive[if (length(alive) == 1) 1 else sample.int(length(alive), 1, prob = rate[alive])]
    if (stats::runif(1) < 1 / (1 + epsilon)) {
      node <- close(i, t, 1)
      new <- length(from) + 1:2
      from[new] <- node
      start[new] <- t
      rate[new] <- rate[i]
      regime[new] <- regime[i]
      alive <- c(alive[alive != i], new)
    } else {
      close(i, t, 2)
      alive <- alive[alive != i]
      # Extinct: intervals with n_tips lineages before extinction still count.
      if (!length(alive)) break
    }
  }
  if (shifts < n_shifts || !length(n_start)) return(NULL)
  phi <- sum(n_end - n_start)
  phi_ref <- gsa_scale / (n_tips * lambda0 * (1 + epsilon))
  if (stats::runif(1) > phi / phi_ref) return(NULL)
  # Cut time: uniform over the intervals with n_tips lineages.
  k <- if (length(n_start) == 1) 1 else sample.int(length(n_start), 1, prob = n_end - n_start)
  cut <- stats::runif(1, n_start[k], n_end[k])
  for (i in alive) close(i, t, 1)
  out <- shift_tree_cut(segments[seq_len(n_seg), , drop = FALSE], cut, rates)
  attr(out, "capped") <- phi > phi_ref
  out
}

# Keep the segments that start before `cut`; segments that span it end there, as living tips.
shift_tree_cut <- function(segments, cut, rates) {
  colnames(segments) <- c("from", "to", "start", "end", "regime", "kind")
  keep <- segments[, "start"] < cut
  seg <- segments[keep, , drop = FALSE]
  spans <- seg[, "end"] > cut
  seg[spans, "end"] <- cut
  seg[spans, "kind"] <- 0
  # Living tips end at the cut; extinct tips (kind 2) are dropped below.
  is_parent <- seg[, "to"] %in% seg[, "from"]
  tip <- !is_parent
  n_tip <- sum(tip)
  ids <- c(0, seg[, "to"])
  internal <- setdiff(ids, seg[tip, "to"])
  number <- integer(max(ids) + 1)
  number[seg[tip, "to"] + 1] <- seq_len(n_tip)
  number[internal + 1] <- n_tip + seq_along(internal)
  tree <- list(edge = cbind(number[seg[, "from"] + 1], number[seg[, "to"] + 1]),
               edge.length = seg[, "end"] - seg[, "start"], Nnode = length(internal),
               tip.label = paste0("t", seq_len(n_tip)))
  class(tree) <- "phylo"
  # Regime of each node: that of the segment ending in it; kept through pruning as labels.
  node_regime <- integer(n_tip + length(internal))
  node_regime[tree$edge[, 2]] <- seg[, "regime"]
  extinct <- tree$edge[seg[, "kind"] == 2 & tip, 2]
  tree$tip.label <- paste0("t", seq_len(n_tip), "_", node_regime[seq_len(n_tip)])
  tree$node.label <- as.character(node_regime[n_tip + seq_along(internal)])
  tree <- ape::drop.tip(tree, extinct, collapse.singles = FALSE)
  tree <- ape::collapse.singles(tree)
  n_tip <- ape::Ntip(tree)
  tip_regime <- as.integer(sub(".*_", "", tree$tip.label))
  lower <- tree$edge[, 2]
  regimes <- ifelse(lower <= n_tip, tip_regime[pmin(lower, n_tip)],
                    as.integer(tree$node.label[pmax(lower - n_tip, 1)]))
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  tree$tip.label <- paste0("t", seq_len(n_tip))
  tree$node.label <- NULL
  list(tree = tree, regimes = factor(regimes), tip_regime = tip_regime, rates = rates)
}
