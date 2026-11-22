#' Internal: trees with heritable shifts in the rate of speciation, for the simulations of the
#' DAMOCLES paper (data-raw/s1_type1_shifts.R and later blocks; dev/plano-artigo-damocles.md).
#'
#' A forward birth-death process with one rate per lineage. Speciation starts at `lambda0` and
#' extinction is `epsilon` times speciation in every lineage. Each of `n_shifts` shifts
#' multiplies the speciation (and so the extinction) rate of one living lineage and all its
#' descendants by `ratio`; shift k happens when the number of living lineages first reaches a
#' number drawn between 3 and `n_tips / 4`, on a lineage drawn at random, so shifts can be
#' nested. The tree is conditioned on `n_tips` living lineages with the general sampling
#' approach (GSA) of Hartmann, Wong & Stadler (2010, Syst. Biol. 59: 465-476,
#' doi:10.1093/sysbio/syq026): each run goes on until `gsa_factor * n_tips` lineages are alive or
#' all are extinct; a run contributes a tree with probability proportional to the total time Phi
#' it spent with exactly `n_tips` lineages (their sampling rate r, here min(1, Phi / Phi_ref) with
#' Phi_ref = gsa_scale / (n_tips * lambda0 * (1 + epsilon)), gsa_scale times the expected length
#' of one such interval under the background rates); the tree is cut in one of those intervals,
#' drawn with weight equal to its length, at a time uniform within it. Extinct lineages are
#' pruned and the tree rescaled to height one. Trees in which a regime (the background included)
#' holds fewer than `min_fraction` of the tips are drawn again (a further condition, not part of
#' the GSA). The attribute "capped" of the result is TRUE when Phi exceeded Phi_ref, which
#' truncates the weight; data-raw/validation_shift_tree.R reports how often.
#'
#' Returns a list with `tree` (a binary ultrametric phylo of height one), `regimes` (a factor,
#' the regime of each edge in the order of tree$edge: "0" for the background and k for the
#' lineages below shift k, the regime of the lower end of the edge), `tip_regime` (integer, per
#' tip) and `rates` (speciation rate of each regime before rescaling).
#' @noRd
sim_shift_tree <- function(n_tips, n_shifts = 1, ratio = 2, lambda0 = 1, epsilon = 0.5,
                           min_fraction = 0.05, gsa_factor = 2, gsa_scale = 15,
                           max_tries = 10000) {
  check_count(n_tips, "n_tips", min = 8)
  check_count(n_shifts, "n_shifts", min = 0)
  check_number(ratio, "ratio", lower = 0, open = c(TRUE, FALSE))
  check_number(lambda0, "lambda0", lower = 0, open = c(TRUE, FALSE))
  check_number(epsilon, "epsilon", 0, 1, open = c(FALSE, TRUE))
  check_number(min_fraction, "min_fraction", 0, 1)
  check_number(gsa_factor, "gsa_factor", lower = 1, open = c(TRUE, FALSE))
  check_number(gsa_scale, "gsa_scale", lower = 0, open = c(TRUE, FALSE))
  check_count(max_tries, "max_tries")
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
  s <- shift_tree_state(lambda0)
  # Intervals with exactly n_tips lineages alive.
  n_start <- numeric(0)
  n_end <- numeric(0)
  while (length(s$alive) < m) {
    while (s$shifts < n_shifts && length(s$alive) >= thresholds[s$shifts + 1]) {
      shift_tree_shift(s, ratio)
    }
    dt <- stats::rexp(1, sum(s$rate[s$alive]) * (1 + epsilon))
    if (length(s$alive) == n_tips) {
      n_start <- c(n_start, s$t)
      n_end <- c(n_end, s$t + dt)
    }
    s$t <- s$t + dt
    shift_tree_event(s, epsilon)
    # Extinct: intervals with n_tips lineages before extinction still count.
    if (!length(s$alive)) break
  }
  if (s$shifts < n_shifts) return(NULL)
  shift_tree_sample(s, n_start, n_end, gsa_scale / (n_tips * lambda0 * (1 + epsilon)))
}

#' the tree of a run `s` cut at a time drawn within the intervals with the number of tips
#' wanted (from `n_start` to `n_end`), or NULL when the run had none or is not sampled: it is
#' kept with probability equal to the total length of those intervals over `phi_ref`
#' @noRd
shift_tree_sample <- function(s, n_start, n_end, phi_ref) {
  if (!length(n_start)) return(NULL)
  phi <- sum(n_end - n_start)
  if (stats::runif(1) > phi / phi_ref) return(NULL)
  # Cut time: uniform over the intervals with n_tips lineages.
  k <- if (length(n_start) == 1) 1 else sample.int(length(n_start), 1, prob = n_end - n_start)
  cut <- stats::runif(1, n_start[k], n_end[k])
  for (i in s$alive) shift_tree_close(s, i, s$t, 1)
  out <- shift_tree_cut(s$segments[seq_len(s$n), , drop = FALSE], cut, s$rates)
  attr(out, "capped") <- phi > phi_ref
  out
}

#' the state of a run of shift_tree_once() at time zero, as an environment that the steps of
#' the run change in place. Lineages: the node they start `from`, their `start` time,
#' speciation `rate` and `regime`, and which are `alive`. `segments` holds the `n` lineages
#' that ended: from, to (node ids, the row numbers), start, end, regime and kind of end (1
#' speciation or shift, 2 extinction). `rates` is the speciation rate of each regime, `t` the
#' time and `shifts` the number of shifts so far
#' @noRd
shift_tree_state <- function(lambda0) {
  s <- new.env(parent = emptyenv())
  s$from <- 0L
  s$start <- 0
  s$rate <- lambda0
  s$regime <- 0L
  s$alive <- 1L
  s$segments <- matrix(0, 64, 6)
  s$n <- 0L
  s$rates <- lambda0
  s$t <- 0
  s$shifts <- 0L
  s
}

#' ends lineage `i` of the run `s` at `time`, with the `kind` of end; returns the id of the
#' node where it ends
#' @noRd
shift_tree_close <- function(s, i, time, kind) {
  n <- s$n + 1L
  if (n > nrow(s$segments)) s$segments <- rbind(s$segments, matrix(0, nrow(s$segments), 6))
  s$segments[n, ] <- c(s$from[i], n, s$start[i], time, s$regime[i], kind)
  s$n <- n
  n
}

#' a shift in the run `s`: the speciation rate of one living lineage, drawn at random, is
#' multiplied by `ratio` and the lineage starts a new regime
#' @noRd
shift_tree_shift <- function(s, ratio) {
  s$shifts <- s$shifts + 1L
  i <- s$alive[sample.int(length(s$alive), 1)]
  node <- shift_tree_close(s, i, s$t, 1)
  s$from[i] <- node
  s$start[i] <- s$t
  s$rate[i] <- s$rate[i] * ratio
  s$regime[i] <- s$shifts
  s$rates[s$shifts + 1] <- s$rate[i]
  invisible(s)
}

#' an event of the run `s` at its current time: one living lineage, drawn with weight equal to
#' its rate, splits in two or goes extinct
#' @noRd
shift_tree_event <- function(s, epsilon) {
  alive <- s$alive
  i <- alive[if (length(alive) == 1) 1 else sample.int(length(alive), 1, prob = s$rate[alive])]
  if (stats::runif(1) < 1 / (1 + epsilon)) {
    node <- shift_tree_close(s, i, s$t, 1)
    new <- length(s$from) + 1:2
    s$from[new] <- node
    s$start[new] <- s$t
    s$rate[new] <- s$rate[i]
    s$regime[new] <- s$regime[i]
    s$alive <- c(alive[alive != i], new)
  } else {
    shift_tree_close(s, i, s$t, 2)
    s$alive <- alive[alive != i]
  }
  invisible(s)
}

#' Keep the segments that start before `cut`; segments that span it end there, as living tips.
#' @noRd
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
