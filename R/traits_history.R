#' A discrete trait evolving along a tree, with its whole history
#'
#' Simulates a discrete character under the Mk model (Lewis 2001): along each branch the state
#' changes after exponential waiting times, with the rates of `q`. The result keeps every change,
#' in the format of a stochastic character map of phytools (class `simmap`), so it is the known
#' history against which reconstructions can be compared, for example by [tip_history_metrics()].
#'
#' Rates can vary among branches in two ways. With `alpha` (or `edge_rates`), all rates of a
#' branch are multiplied by one factor, drawn for each branch from a gamma distribution with
#' shape and rate `alpha` (mean 1), as in the model that `phytools::fitgammaMk()` fits (Revell
#' 2024). With `rate_classes`, the character has hidden rate classes, as in the hidden-rates
#' models of corHMM (Beaulieu et al. 2013): in class k, changes between states happen at the
#' rates of `q` times `rate_classes[k]`, and the lineage moves between classes at rate
#' `class_switch` without changing its state. The class at the root is drawn with equal
#' probabilities.
#'
#' @param tree A rooted `phylo` tree with branch lengths.
#' @param root Deprecated; use `root_state`.
#' @param q One rate for every change between states (the equal-rates model), or a square matrix
#'   of rates from the state of each row to the state of each column, with the states as names;
#'   its diagonal is ignored.
#' @param states Names of the states, when `q` is one rate.
#' @param root_state State at the root; `NULL` draws one with equal probabilities.
#' @param alpha Shape (and rate) of the gamma distribution of the rate multiplier of each branch;
#'   `NULL` for none. Small values make rates more uneven among branches.
#' @param edge_rates Optional multipliers of the rates of each branch, in the order of the rows of
#'   `tree$edge`, instead of drawing them with `alpha`.
#' @param rate_classes Optional multipliers of the rates in each hidden class, such as
#'   `c(slow = 0.1, fast = 1)`; classes without names are called `R1`, `R2` and so on.
#' @param class_switch Rate of moving from one hidden class to each other one.
#' @return The tree, reordered "cladewise", with class `c("simmap", "phylo")` and, as in
#'   `phytools::sim.history()`: `maps`, a list with one element per edge giving the time spent in
#'   each state in the order they occurred; `mapped.edge`, the total time in each state per edge;
#'   `node.states`, the states at the start and end of each edge; and `states`, the states of
#'   the tips. `Q` holds the rate matrix used. With rates varying among branches, `edge_rates`
#'   holds the multiplier of each edge, in the new order of the edges. With `rate_classes`,
#'   `hidden` holds the same history with the class of each segment, as states named
#'   `"state|class"`, with its own `maps`, `node.states`, `states` and `Q`, and `classes` holds
#'   the class of each tip.
#' @references Beaulieu, J. M., O'Meara, B. C., & Donoghue, M. J. (2013). Identifying hidden rate
#'   changes in the evolution of a binary morphological character: the evolution of plant habit
#'   in campanulid angiosperms. *Systematic Biology*, 62, 725--737. \doi{10.1093/sysbio/syt034}
#'
#'   Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from discrete
#'   morphological character data. *Systematic Biology*, 50, 913--925.
#'   \doi{10.1080/106351501753462876}
#'
#'   Revell, L. J. (2024). Testing for heterogeneous rates of discrete character evolution on
#'   phylogenies. *Journal of Evolutionary Biology*, 37, 1591--1602. \doi{10.1093/jeb/voae119}
#' @seealso [tip_history_metrics()], [trait_rates()] for continuous traits, and [trait_ou()],
#'   which takes the mapped history as `regimes` (optima and rates that change with the state).
#' @family species traits
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(10)
#' h <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
#' h$states
#' tip_history_metrics(h)
#' # rates that vary among branches
#' g <- trait_mk(tr, q = 0.5, alpha = 0.3)
#' round(g$edge_rates, 2)
#' hr <- trait_mk(tr, q = 0.5, rate_classes = c(slow = 0.1, fast = 2), class_switch = 0.2)
#' hr$classes
#' @export
trait_mk <- function(tree, q = 1, states = c("0", "1"), root_state = NULL, alpha = NULL,
                     edge_rates = NULL, rate_classes = NULL, class_switch = 0.1, root = NULL) {
  check_tree(tree)
  root <- renamed_arg(root_state, root, "root_state", "root")
  if (is.matrix(q)) {
    states <- rownames(q)
    if (is.null(states) || !identical(states, colnames(q))) {
      stop("`q` needs the same state names in its rows and columns", call. = FALSE)
    }
  } else {
    stopifnot(length(q) == 1, length(states) >= 2, !anyDuplicated(states))
    q <- matrix(q, length(states), length(states), dimnames = list(states, states))
  }
  diag(q) <- 0
  if (any(!is.finite(q)) || any(q < 0)) stop("rates in `q` must be zero or more", call. = FALSE)
  if (!is.null(alpha) && !is.null(edge_rates)) {
    stop("give `alpha` or `edge_rates`, not both", call. = FALSE)
  }
  n_edge <- nrow(tree$edge)
  if (!is.null(alpha)) {
    if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("`alpha` must be one positive number", call. = FALSE)
    }
    edge_rates <- stats::rgamma(n_edge, shape = alpha, rate = alpha)
  }
  if (!is.null(edge_rates) &&
      (length(edge_rates) != n_edge || any(!is.finite(edge_rates)) || any(edge_rates < 0))) {
    stop("`edge_rates` needs one number, zero or more, per row of `tree$edge`", call. = FALSE)
  }
  observed_q <- q
  diag(observed_q) <- -rowSums(q)
  root <- root %||% sample(states, 1)
  if (!root %in% states) stop("`root_state` must be one of the states", call. = FALSE)

  if (!is.null(rate_classes)) {
    if (!is.numeric(rate_classes) || length(rate_classes) < 2 || any(rate_classes < 0)) {
      stop("`rate_classes` needs two or more multipliers, zero or more", call. = FALSE)
    }
    if (length(class_switch) != 1 || !is.finite(class_switch) || class_switch < 0) {
      stop("`class_switch` must be one number, zero or more", call. = FALSE)
    }
    classes <- names(rate_classes) %||% paste0("R", seq_along(rate_classes))
    if (anyDuplicated(classes) || any(classes == "")) {
      stop("`rate_classes` needs different names", call. = FALSE)
    }
    names(rate_classes) <- classes
    full <- expand.grid(state = states, class = classes, stringsAsFactors = FALSE)
    labels <- paste(full$state, full$class, sep = "|")
    hidden_q <- matrix(0, nrow(full), nrow(full), dimnames = list(labels, labels))
    for (i in seq_len(nrow(full))) for (j in seq_len(nrow(full))) {
      if (full$class[i] == full$class[j]) {
        hidden_q[i, j] <- q[full$state[i], full$state[j]] * rate_classes[[full$class[i]]]
      } else if (full$state[i] == full$state[j]) {
        hidden_q[i, j] <- class_switch
      }
    }
    diag(hidden_q) <- 0
    start <- paste(root, sample(classes, 1), sep = "|")
  } else {
    hidden_q <- q
    start <- root
  }

  order <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  tree <- ape::reorder.phylo(tree, "cladewise")
  if (!is.null(edge_rates)) edge_rates <- edge_rates[order]
  run <- mk_history(tree, hidden_q, start, edge_rates %||% rep(1, n_edge))
  if (is.null(rate_classes)) {
    out <- as_simmap(tree, run$maps, run$at_node, states)
  } else {
    observed <- stats::setNames(full$state, labels)
    collapsed <- lapply(run$maps, function(m) {
      runs <- rle(unname(observed[names(m)]))
      stats::setNames(vapply(split(unname(m), rep(seq_along(runs$lengths), runs$lengths)), sum,
                             numeric(1)), runs$values)
    })
    out <- as_simmap(tree, collapsed, unname(observed[run$at_node]), states)
    hidden <- as_simmap(tree, run$maps, run$at_node, labels)
    diag(hidden_q) <- -rowSums(hidden_q)
    hidden$Q <- hidden_q
    out$hidden <- hidden
    out$classes <- stats::setNames(full$class[match(run$at_node[seq_len(ape::Ntip(tree))],
                                                    labels)], tree$tip.label)
  }
  out$Q <- observed_q
  if (!is.null(edge_rates)) out$edge_rates <- edge_rates
  out
}

# the history of a continuous-time Markov chain with rates `q` (zero diagonal) along a tree in
# cladewise order, starting from `root`; the rates of edge e are multiplied by `multiplier[e]`
mk_history <- function(tree, q, root, multiplier) {
  states <- rownames(q)
  n_tip <- ape::Ntip(tree)
  at_node <- character(n_tip + tree$Nnode)
  at_node[n_tip + 1] <- root
  leaving <- rowSums(q)
  maps <- vector("list", nrow(tree$edge))
  for (e in seq_len(nrow(tree$edge))) {  # cladewise: parents before children
    state <- at_node[tree$edge[e, 1]]
    left <- tree$edge.length[e]
    times <- numeric(0)
    visited <- character(0)
    repeat {
      rate <- leaving[[state]] * multiplier[e]
      wait <- if (rate > 0) stats::rexp(1, rate) else Inf
      if (wait >= left) {
        times <- c(times, left)
        visited <- c(visited, state)
        break
      }
      times <- c(times, wait)
      visited <- c(visited, state)
      left <- left - wait
      others <- setdiff(states, state)
      state <- if (length(others) == 1) others else sample(others, 1, prob = q[state, others])
    }
    maps[[e]] <- stats::setNames(times, visited)
    at_node[tree$edge[e, 2]] <- state
  }
  list(maps = maps, at_node = at_node)
}

# a cladewise tree with the history of a discrete character, in the format of phytools
as_simmap <- function(tree, maps, at_node, states) {
  n_tip <- ape::Ntip(tree)
  labels <- paste(tree$edge[, 1], tree$edge[, 2], sep = ",")
  mapped <- t(vapply(maps, function(m) vapply(states, function(s) sum(m[names(m) == s]),
                                                numeric(1)), numeric(length(states))))
  dimnames(mapped) <- list(labels, state = states)
  tree$maps <- maps
  tree$node.states <- cbind(at_node[tree$edge[, 1]], at_node[tree$edge[, 2]])
  tree$states <- stats::setNames(at_node[seq_len(n_tip)], tree$tip.label)
  tree$mapped.edge <- mapped
  class(tree) <- c("simmap", "phylo")
  attr(tree, "order") <- "cladewise"
  tree
}

#' Tip-based metrics of a known trait history
#'
#' Computes, for each tip, the metrics of trait evolution of Luza et al. (2021) from a history of
#' a discrete trait, as `Herodotools::calc_tip_based_trait_evo()` computes them from a stochastic
#' character map. Applied to a history from [trait_mk()], they are the true values that a
#' reconstruction from the tips should recover. Two metrics in continuous time are added, which
#' count the changes inside branches too.
#'
#' Following the code of `calc_tip_based_trait_evo()`, the state of a node is the state at the
#' start of the edges that leave it, and the path of a tip runs from its parent node to the root.
#' `transitions` counts the changes of state between consecutive nodes of the path, starting
#' from the state of the tip; `prop_transitions` divides it by the number of nodes in the path,
#' `total_nodes`. `last_transition_time` sums the length of the terminal edge and of the edges
#' leading to the nodes of the path, from the tip towards the root, while their state is the
#' state of the tip; the root, which has no edge leading to it, adds nothing. `stasis_time` is
#' the longest time spent in one state within a single edge, over the terminal edge and the
#' edges leading to the nodes of the path. Only nodes enter `transitions` and
#' `last_transition_time`, so changes that revert within a branch are missed; `changes` counts
#' every change along the path from the root to the tip and `time_since_change` is the time
#' since the last of them, or the age of the tip from the root when there was none.
#'
#' @param history A `simmap` tree, from [trait_mk()] or phytools.
#' @return Data frame with one row per tip and the columns `transitions`, `total_nodes`,
#'   `prop_transitions`, `last_transition_time`, `stasis_time`, `changes` and
#'   `time_since_change`.
#' @references Luza, A. L., Maestri, R., Debastiani, V. J., Patterson, B. D., Hartz, S. M., &
#'   Duarte, L. D. S. (2021). Is evolution faster at ecotones? A test using rates and tempo of
#'   diet transitions in Neotropical Sigmodontinae (Rodentia, Cricetidae). *Ecology and
#'   Evolution*, 11, 18676--18690. \doi{10.1002/ece3.8476}
#' @seealso [trait_mk()]
#' @inherit trait_mk examples
#' @family species traits
#' @export
tip_history_metrics <- function(history) {
  states <- history_states(history)
  edge <- history$edge
  maps <- history$maps
  n_tip <- ape::Ntip(history)
  node_state <- c(unname(states$tips), states$nodes)
  tip_state <- unname(states$tips)
  edge_to <- match(seq_len(n_tip + history$Nnode), edge[, 2])  # NA at the root
  longest <- vapply(maps, max, numeric(1))
  rows <- lapply(seq_len(n_tip), function(i) {
    edges <- integer(0)
    node <- i
    while (!is.na(edge_to[node])) {
      edges <- c(edges, edge_to[node])
      node <- edge[edge_to[node], 1]
    }
    path <- edge[edges, 1]  # nodes from the parent of the tip to the root
    states <- node_state[path]
    sequence <- c(tip_state[i], states)
    transitions <- sum(sequence[-1] != sequence[-length(sequence)])
    same <- cumprod(states == tip_state[i]) == 1
    lengths <- c(history$edge.length[edges[-1]], 0)  # edges leading to the nodes of the path
    last <- history$edge.length[edges[1]] + sum(lengths[same])
    segments <- unlist(lapply(edges, function(e) rev(maps[[e]])), use.names = TRUE)
    first_change <- which(names(segments) != tip_state[i])[1]
    since <- if (is.na(first_change)) sum(segments) else sum(segments[seq_len(first_change - 1)])
    data.frame(transitions = transitions, total_nodes = length(path),
               prop_transitions = transitions / length(path), last_transition_time = last,
               stasis_time = max(longest[edges]),
               changes = sum(vapply(maps[edges], length, integer(1)) - 1L),
               time_since_change = since)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- history$tip.label
  out
}

# values at the internal nodes (root included) of a Brownian motion with Pagel's lambda, drawn
# given the values at the tips, so that tips and nodes have their joint distribution
bm_nodes <- function(tree, lambda, tips) {
  n_tip <- ape::Ntip(tree)
  scaled <- tree
  terminal <- tree$edge[, 2] <= n_tip
  depth <- ape::node.depth.edgelength(tree)
  scaled$edge.length[!terminal] <- lambda * tree$edge.length[!terminal]
  scaled$edge.length[terminal] <- tree$edge.length[terminal] +
    (1 - lambda) * depth[tree$edge[terminal, 1]]
  covariance <- matrix(ape::node.depth.edgelength(scaled)[ape::mrca(scaled, full = TRUE)],
                       n_tip + tree$Nnode)
  t <- seq_len(n_tip)
  n <- n_tip + seq_len(tree$Nnode)
  weights <- covariance[n, t, drop = FALSE] %*% solve(covariance[t, t])
  centre <- drop(weights %*% tips[tree$tip.label])
  spread <- covariance[n, n, drop = FALSE] - weights %*% covariance[t, n, drop = FALSE]
  decomposition <- eigen((spread + t(spread)) / 2, symmetric = TRUE)
  noise <- decomposition$vectors %*% (sqrt(pmax(decomposition$values, 0)) *
                                        stats::rnorm(tree$Nnode))
  stats::setNames(centre + drop(noise), n)
}
