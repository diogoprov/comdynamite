#' A continuous trait whose rate of evolution varies among branches
#'
#' Simulates Brownian motion in which each branch has its own rate (variance per unit time),
#' and keeps those rates, so that methods which estimate rates per branch can be compared with
#' the rates that generated the data. Three models of how rates vary are available:
#'
#' * `"shift"`: the rate is `sigma2` at the root and changes at a few points of the tree. With
#'   `shift_branch = FALSE` the new rate is inherited by every branch descending from the
#'   shift, as in the models of Eastman et al. (2011), traitMedusa in motmot (Thomas &
#'   Freckleton 2012), BAMM (Rabosky 2014) or bifrost; with `shift_branch = TRUE` it lasts only
#'   to the end of the branch. By default `shifts` shifts happen at the start of branches drawn
#'   with probabilities proportional to their length, or of the branches in `shift_edges`. With
#'   `shift_rate`, shifts instead happen as a Poisson process along the branches, as under the
#'   prior of BAMM: their number is Poisson with mean `shift_rate` times the total length of the
#'   tree, and they fall at uniform positions over it, so also within branches. At each shift
#'   the rate is multiplied by `shift_size` (`relative = TRUE`), or becomes `sigma2` times
#'   `shift_size` whatever it was before (`relative = FALSE`), as each new regime of BAMM or
#'   bifrost has a rate of its own. `shift_size` can be a function of the number of shifts that
#'   draws them, such as `function(n) stats::rlnorm(n, 0, 1)`. The rate of a branch is its mean
#'   over the branch, weighted by time.
#' * `"gbm"`: the logarithm of the rate evolves itself by Brownian motion with standard
#'   deviation `rate_sd` per square root of time, starting from `log(sigma2)` at the root, as in
#'   the variable-rate models of Revell (2021) and Martin et al. (2023). The rate of a branch is
#'   the geometric mean of the rates at its two ends.
#' * `"iid"`: the rate of each branch is `sigma2` times an independent lognormal multiplier with
#'   median 1 and `rate_sd` as the standard deviation of its logarithm.
#'
#' The trait is then Brownian motion on the tree whose branch lengths are multiplied by their
#' rates, drawn at the tips and then at the internal nodes given the tips, as in [trait_bm()].
#'
#' @param tree A rooted `phylo` tree with branch lengths.
#' @param model How rates vary among branches: `"shift"`, `"gbm"` or `"iid"`.
#' @param sigma2 Rate at the root (`"shift"`, `"gbm"`) or median rate (`"iid"`).
#' @param shifts Number of shifts, for `"shift"` when neither `shift_edges` nor `shift_rate` is
#'   given.
#' @param shift_size Multipliers of the rate at each shift, recycled to the number of shifts, or
#'   a function of the number of shifts that gives them. Values below 1 slow evolution down.
#' @param shift_edges Optional row numbers of `tree$edge` at whose start the shifts happen.
#' @param shift_rate Optional rate of shifts per unit of branch length, for shifts as a Poisson
#'   process along the branches.
#' @param shift_branch Should a shift last only to the end of its branch, rather than pass to
#'   all branches descending from it?
#' @param relative Multiply the rate before each shift by `shift_size` (`TRUE`), or set the rate
#'   after it to `sigma2` times `shift_size` (`FALSE`)?
#' @param rate_sd Standard deviation of the logarithm of the rate: per square root of time for
#'   `"gbm"` and among branches for `"iid"`.
#' @param root_value Trait value at the root.
#' @param edge_rates Optional rate of each branch, in the order of the rows of `tree$edge`, such as
#'   those of [history_rates()]. When given, `model` and its arguments are ignored.
#' @return A list with
#'   * `tips`: the trait at the tips, named by species;
#'   * `nodes`: the trait at the internal nodes, named by their numbers in `tree$edge`, the
#'     root first;
#'   * `edge_rate`: the rate of each branch, in the order of the rows of `tree$edge`;
#'   * `node_rate`: for `"gbm"`, the rate at every node, tips first, named by node number;
#'     `NULL` otherwise;
#'   * `shifts`: for `"shift"`, a data frame with the `edge` (row of `tree$edge`), the `node`
#'     at its end, the `position` of the shift as time from the start of the edge, the
#'     `multiplier` and the `rate` after the shift; `NULL` otherwise;
#'   * `tree`: the tree given, and `scaled`: the tree with each branch length multiplied by its
#'     rate, under which the trait is ordinary Brownian motion.
#' @references Rabosky, D. L. (2014). Automatic detection of key innovations, rate shifts, and
#'   diversity-dependence on phylogenetic trees. *PLoS ONE*, 9, e89543.
#'   \doi{10.1371/journal.pone.0089543}
#'
#'   Eastman, J. M., Alfaro, M. E., Joyce, P., Hipp, A. L., & Harmon, L. J. (2011). A
#'   novel comparative method for identifying shifts in the rate of character evolution on
#'   trees. *Evolution*, 65, 3578--3589. \doi{10.1111/j.1558-5646.2011.01401.x}
#'
#'   Martin, B. S., Bradburd, G. S., Harmon, L. J., & Weber, M. G. (2023). Modeling the evolution
#'   of rates of continuous trait evolution. *Systematic Biology*, 72, 590--605.
#'   \doi{10.1093/sysbio/syac068}
#'
#'   Revell, L. J. (2021). A variable-rate quantitative trait evolution model using
#'   penalized-likelihood. *PeerJ*, 9, e11997. \doi{10.7717/peerj.11997}
#'
#'   Thomas, G. H., & Freckleton, R. P. (2012). MOTMOT: models of trait macroevolution on trees.
#'   *Methods in Ecology and Evolution*, 3, 145--151. \doi{10.1111/j.2041-210X.2011.00132.x}
#' @seealso [trait_bm()] for a single rate, [trait_mk()] for discrete traits with rates that
#'   vary among branches.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(20)
#' x <- trait_rates(tr, "shift", shifts = 1, shift_size = 8)
#' x$shifts
#' head(x$tips)
#' # shifts as a Poisson process, each new regime with a rate of its own, as in BAMM
#' b <- trait_rates(tr, "shift", shift_rate = 0.1, relative = FALSE,
#'                  shift_size = function(n) stats::rexp(n))
#' b$shifts
#' g <- trait_rates(tr, "gbm", rate_sd = 0.5)
#' summary(g$edge_rate)
#' # as the environmental optima of a metacommunity
#' traits <- species_traits(tr, optima = g$tips)
#' @export
trait_rates <- function(tree, model = c("shift", "gbm", "iid"), sigma2 = 1, shifts = 1,
                        shift_size = 4, shift_edges = NULL, shift_rate = NULL,
                        shift_branch = FALSE, relative = TRUE, rate_sd = 0.5, root_value = 0,
                        edge_rates = NULL) {
  stopifnot(inherits(tree, "phylo"), !is.null(tree$edge.length))
  model <- if (is.null(edge_rates)) match.arg(model) else "given"
  if (length(sigma2) != 1 || !is.finite(sigma2) || sigma2 <= 0) {
    stop("`sigma2` must be one positive number", call. = FALSE)
  }
  if (length(rate_sd) != 1 || !is.finite(rate_sd) || rate_sd < 0) {
    stop("`rate_sd` must be one number, zero or more", call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  n_edge <- nrow(tree$edge)
  cladewise <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)  # parents first
  root <- n_tip + 1
  node_rate <- NULL
  shift_table <- NULL

  if (model == "given") {
    if (!is.numeric(edge_rates) || length(edge_rates) != n_edge || anyNA(edge_rates) ||
        any(edge_rates <= 0)) {
      stop("`edge_rates` needs one positive number per row of `tree$edge`", call. = FALSE)
    }
    edge_rate <- unname(edge_rates)
  } else if (model == "shift") {
    if (!is.null(shift_rate)) {
      if (length(shift_rate) != 1 || !is.finite(shift_rate) || shift_rate < 0) {
        stop("`shift_rate` must be one number, zero or more", call. = FALSE)
      }
      n_shifts <- stats::rpois(1, shift_rate * sum(tree$edge.length))
      shift_edges <- sample.int(n_edge, n_shifts, replace = TRUE, prob = tree$edge.length)
      position <- stats::runif(n_shifts) * tree$edge.length[shift_edges]
    } else {
      if (is.null(shift_edges)) {
        stopifnot(length(shifts) == 1, shifts >= 0, shifts <= n_edge)
        shift_edges <- if (shifts == 0) integer(0) else {
          sample.int(n_edge, shifts, prob = tree$edge.length)
        }
      }
      shift_edges <- as.integer(shift_edges)
      if (anyNA(shift_edges) || any(shift_edges < 1 | shift_edges > n_edge) ||
          anyDuplicated(shift_edges)) {
        stop("`shift_edges` must be different row numbers of `tree$edge`", call. = FALSE)
      }
      position <- numeric(length(shift_edges))
    }
    n_shifts <- length(shift_edges)
    size <- if (is.function(shift_size)) shift_size(n_shifts) else {
      rep_len(shift_size, n_shifts)
    }
    if (n_shifts > 0 && (!is.numeric(size) || length(size) != n_shifts || anyNA(size) ||
                         any(size <= 0))) {
      stop("`shift_size` must give one positive multiplier per shift", call. = FALSE)
    }
    at_node <- rep(sigma2, n_tip + tree$Nnode)
    edge_rate <- numeric(n_edge)
    after <- numeric(n_shifts)
    for (e in cladewise) {
      start <- at_node[tree$edge[e, 1]]
      here <- which(shift_edges == e)
      here <- here[order(position[here])]
      rate <- start
      integral <- 0
      from <- 0
      for (k in here) {
        integral <- integral + rate * (position[k] - from)
        rate <- if (relative) rate * size[k] else sigma2 * size[k]
        after[k] <- rate
        from <- position[k]
      }
      integral <- integral + rate * (tree$edge.length[e] - from)
      edge_rate[e] <- if (tree$edge.length[e] > 0) integral / tree$edge.length[e] else rate
      at_node[tree$edge[e, 2]] <- if (shift_branch) start else rate
    }
    shift_table <- data.frame(edge = shift_edges, node = tree$edge[shift_edges, 2],
                              position = position, multiplier = size, rate = after)
    shift_table <- shift_table[order(shift_table$edge, shift_table$position), , drop = FALSE]
    rownames(shift_table) <- NULL
  } else if (model == "gbm") {
    log_rate <- numeric(n_tip + tree$Nnode)
    log_rate[root] <- log(sigma2)
    steps <- stats::rnorm(n_edge, 0, rate_sd * sqrt(tree$edge.length))
    for (e in cladewise) log_rate[tree$edge[e, 2]] <- log_rate[tree$edge[e, 1]] + steps[e]
    edge_rate <- exp((log_rate[tree$edge[, 1]] + log_rate[tree$edge[, 2]]) / 2)
    node_rate <- stats::setNames(exp(log_rate), seq_along(log_rate))
  } else {
    edge_rate <- sigma2 * exp(stats::rnorm(n_edge, 0, rate_sd))
  }

  scaled <- tree
  scaled$edge.length <- tree$edge.length * edge_rate
  tips <- trait_bm(scaled, 1, nodes = TRUE)
  nodes <- attr(tips, "nodes") + root_value
  attr(tips, "nodes") <- NULL
  list(tips = tips + root_value, nodes = nodes, edge_rate = edge_rate, node_rate = node_rate,
       shifts = shift_table, tree = tree, scaled = scaled)
}

#' Rates of trait evolution from a history along the tree
#'
#' Turns the history of a discrete character along a tree, such as the habitat or region of
#' each lineage from [trait_mk()] or a stochastic map of phytools, into the rate of evolution of
#' a continuous trait on each branch, for `trait_rates(edge_rates = )`. The rate is a function
#' of the history at each moment, and the rate of a branch is its mean over the branch, weighted
#' by time. Three models link the history to the rate:
#'
#' * `"state"`: each state has its own rate, `state_rates`, as in Brownian motion with
#'   several regimes (O'Meara et al. 2006).
#' * `"burst"`: each change of state, such as the colonisation of a new habitat or region,
#'   multiplies the rate by `1 + burst`, and the excess decays exponentially with time scale
#'   `decay`, also along the descending branches, until the next change. The rate at time `u`
#'   after the last change is `sigma2 * (1 + burst * exp(-u / decay))`; before any change it is
#'   `sigma2`. This follows the expectation of McPeek (2017) that selection, and so evolution,
#'   is strongest when a lineage first adapts to a new community.
#' * `"sympatry"`: the rate depends on the number of lineages that share the state at each
#'   moment, counting the lineage itself, as `sigma2 * exp(r * n)`; a negative `r` slows
#'   evolution down as lineages accumulate in the same region, as in the diversity-dependent model with
#'   biogeography of Drury et al. (2016).
#'
#' @param history A tree with the history of a discrete character, of class `simmap`, from
#'   [trait_mk()] or phytools.
#' @param model `"state"`, `"burst"` or `"sympatry"`.
#' @param sigma2 Rate before any change of state (`"burst"`) or scale of the rate
#'   (`"sympatry"`).
#' @param state_rates For `"state"`, the rate in each state, named by state.
#' @param burst,decay For `"burst"`, the relative increase of the rate at each change of state
#'   and the time scale of its decay.
#' @param r For `"sympatry"`, the change in the logarithm of the rate per lineage in the same
#'   state.
#' @return Numeric vector with the rate of each branch, in the order of the rows of
#'   `history$edge`.
#' @references Drury, J., Clavel, J., Manceau, M., & Morlon, H. (2016). Estimating the effect of
#'   competition on trait evolution using maximum likelihood inference. *Systematic Biology*, 65,
#'   700--710. \doi{10.1093/sysbio/syw020}
#'
#'   McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton University Press.
#'   \doi{10.23943/princeton/9780691088778.001.0001}
#'
#'   O'Meara, B. C., Ané, C., Sanderson, M. J., & Wainwright, P. C. (2006). Testing for different
#'   rates of continuous trait evolution using likelihood. *Evolution*, 60, 922--933.
#'   \doi{10.1111/j.0014-3820.2006.tb01171.x}
#' @seealso [trait_rates()], [trait_mk()], [as_herodotools()]
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(30)
#' habitat <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
#' rates <- history_rates(habitat, "state", state_rates = c(forest = 1, open = 5))
#' x <- trait_rates(habitat, edge_rates = rates)
#' # a burst of evolution after each change of habitat
#' b <- history_rates(habitat, "burst", burst = 9, decay = 0.1)
#' summary(b)
#' # slower evolution as lineages accumulate in the same habitat
#' s <- history_rates(habitat, "sympatry", r = -0.1)
#' @export
history_rates <- function(history, model = c("state", "burst", "sympatry"), sigma2 = 1,
                          state_rates = NULL, burst = 4, decay = 0.1, r = -0.1) {
  if (!inherits(history, "phylo") || is.null(history$maps)) {
    stop("`history` must be a tree with a history of class `simmap`", call. = FALSE)
  }
  model <- match.arg(model)
  if (length(sigma2) != 1 || !is.finite(sigma2) || sigma2 <= 0) {
    stop("`sigma2` must be one positive number", call. = FALSE)
  }
  maps <- history$maps
  lengths <- history$edge.length
  states <- unique(unlist(lapply(maps, names), use.names = FALSE))
  integral <- switch(model,
    state = {
      if (is.null(state_rates) || is.null(names(state_rates)) ||
          !all(states %in% names(state_rates))) {
        stop("`state_rates` needs a rate named by each state: ",
             paste(states, collapse = ", "), call. = FALSE)
      }
      if (any(!is.finite(state_rates[states])) || any(state_rates[states] <= 0)) {
        stop("`state_rates` must be positive", call. = FALSE)
      }
      vapply(maps, function(m) sum(m * state_rates[names(m)]), numeric(1))
    },
    burst = {
      if (length(burst) != 1 || !is.finite(burst) || burst < 0 ||
          length(decay) != 1 || !is.finite(decay) || decay <= 0) {
        stop("`burst` must be zero or more and `decay` positive", call. = FALSE)
      }
      burst_integral(history, burst, decay) * sigma2
    },
    sympatry = {
      if (length(r) != 1 || !is.finite(r)) stop("`r` must be one number", call. = FALSE)
      sympatry_integral(history, r) * sigma2
    })
  rate <- integral / lengths
  zero <- lengths <= 0  # without length the rate does not matter; keep it finite and positive
  if (any(zero)) {
    ends <- vapply(maps[zero], function(m) names(m)[length(m)], character(1))
    rate[zero] <- if (model == "state") state_rates[ends] else sigma2
  }
  unname(rate)
}

# integral of 1 + burst * exp(-u / decay) along each edge, where u is the time since the last
# change of state on the path from the root
burst_integral <- function(history, burst, decay) {
  edge <- history$edge
  since <- rep(Inf, ape::Ntip(history) + history$Nnode)
  integral <- numeric(nrow(edge))
  for (e in ape::reorder.phylo(history, "cladewise", index.only = TRUE)) {
    u <- since[edge[e, 1]]
    m <- history$maps[[e]]
    total <- 0
    for (k in seq_along(m)) {
      if (k > 1) u <- 0  # a change of state at the start of this segment
      span <- m[[k]]
      total <- total + span + burst * decay * (exp(-u / decay) - exp(-(u + span) / decay))
      u <- u + span
    }
    integral[e] <- total
    since[edge[e, 2]] <- u
  }
  integral
}

# integral of exp(r * n) along each edge, where n is the number of lineages in the same state
# at each moment, the lineage itself included
sympatry_integral <- function(history, r) {
  edge <- history$edge
  depth <- ape::node.depth.edgelength(history)
  pieces <- do.call(rbind, lapply(seq_len(nrow(edge)), function(e) {
    m <- history$maps[[e]]
    end <- depth[edge[e, 1]] + cumsum(unname(m))
    data.frame(edge = e, state = names(m), from = c(depth[edge[e, 1]], end[-length(end)]),
               to = end, stringsAsFactors = FALSE)
  }))
  pieces <- pieces[pieces$to > pieces$from, , drop = FALSE]
  breaks <- sort(unique(c(pieces$from, pieces$to)))
  integral <- numeric(nrow(edge))
  for (i in seq_len(length(breaks) - 1)) {
    alive <- which(pieces$from <= breaks[i] & pieces$to >= breaks[i + 1])
    n <- table(pieces$state[alive])[pieces$state[alive]]
    add <- exp(r * as.numeric(n)) * (breaks[i + 1] - breaks[i])
    by_edge <- rowsum(add, pieces$edge[alive])
    at <- as.integer(rownames(by_edge))
    integral[at] <- integral[at] + by_edge[, 1]
  }
  integral
}
