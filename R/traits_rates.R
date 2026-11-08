#' A continuous trait whose rate of evolution varies among branches
#'
#' Simulates Brownian motion in which each branch has its own rate (variance per unit time),
#' and keeps those rates, so that methods which estimate rates per branch can be compared with
#' the rates that generated the data. Three models of how rates vary are available:
#'
#' * `"shift"`: the rate is `sigma2` at the root and is multiplied by `shift_size` at a few
#'   branches. With `shift_branch = FALSE` the new rate is inherited by every branch descending
#'   from the shifted one, as in the models of Eastman et al. (2011), traitMedusa in motmot
#'   (Thomas & Freckleton 2012) or bifrost; with `shift_branch = TRUE` only the shifted branch
#'   changes. The branches are drawn with probabilities proportional to their length, or given
#'   in `shift_edges`.
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
#' @param shifts Number of shifts, for `"shift"` when `shift_edges` is not given.
#' @param shift_size Multipliers of the rate at each shift, recycled to the number of shifts.
#'   Values below 1 slow evolution down.
#' @param shift_edges Optional row numbers of `tree$edge` where the shifts happen.
#' @param shift_branch Should a shift change only its own branch, rather than all branches
#'   descending from it?
#' @param rate_sd Standard deviation of the logarithm of the rate: per square root of time for
#'   `"gbm"` and among branches for `"iid"`.
#' @param root_value Trait value at the root.
#' @return A list with
#'   * `tips`: the trait at the tips, named by species;
#'   * `nodes`: the trait at the internal nodes, named by their numbers in `tree$edge`, the
#'     root first;
#'   * `edge_rate`: the rate of each branch, in the order of the rows of `tree$edge`;
#'   * `node_rate`: for `"gbm"`, the rate at every node, tips first, named by node number;
#'     `NULL` otherwise;
#'   * `shifts`: for `"shift"`, a data frame with the `edge` (row of `tree$edge`), the `node`
#'     at its end and the `multiplier`; `NULL` otherwise;
#'   * `tree`: the tree given, and `scaled`: the tree with each branch length multiplied by its
#'     rate, under which the trait is ordinary Brownian motion.
#' @references Eastman, J. M., Alfaro, M. E., Joyce, P., Hipp, A. L., & Harmon, L. J. (2011). A
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
#' g <- trait_rates(tr, "gbm", rate_sd = 0.5)
#' summary(g$edge_rate)
#' # as the environmental optima of a metacommunity
#' traits <- species_traits(tr, optima = g$tips)
#' @export
trait_rates <- function(tree, model = c("shift", "gbm", "iid"), sigma2 = 1, shifts = 1,
                        shift_size = 4, shift_edges = NULL, shift_branch = FALSE,
                        rate_sd = 0.5, root_value = 0) {
  stopifnot(inherits(tree, "phylo"), !is.null(tree$edge.length))
  model <- match.arg(model)
  if (length(sigma2) != 1 || !is.finite(sigma2) || sigma2 <= 0) {
    stop("`sigma2` must be one positive number", call. = FALSE)
  }
  if (length(rate_sd) != 1 || !is.finite(rate_sd) || rate_sd < 0) {
    stop("`rate_sd` must be one number, zero or more", call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  n_edge <- nrow(tree$edge)
  order <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)  # parents before children
  root <- n_tip + 1
  node_rate <- NULL
  shift_table <- NULL

  if (model == "shift") {
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
    if (length(shift_edges) > 0 &&
        (!is.numeric(shift_size) || length(shift_size) == 0 || any(shift_size <= 0))) {
      stop("`shift_size` must be positive", call. = FALSE)
    }
    multiplier <- rep(1, n_edge)
    multiplier[shift_edges] <- rep_len(shift_size, length(shift_edges))
    at_node <- rep(1, n_tip + tree$Nnode)
    relative <- numeric(n_edge)
    for (e in order) {
      parent <- at_node[tree$edge[e, 1]]
      relative[e] <- parent * multiplier[e]
      at_node[tree$edge[e, 2]] <- if (shift_branch) parent else relative[e]
    }
    edge_rate <- sigma2 * relative
    shift_table <- data.frame(edge = shift_edges, node = tree$edge[shift_edges, 2],
                              multiplier = multiplier[shift_edges])
  } else if (model == "gbm") {
    log_rate <- numeric(n_tip + tree$Nnode)
    log_rate[root] <- log(sigma2)
    steps <- stats::rnorm(n_edge, 0, rate_sd * sqrt(tree$edge.length))
    for (e in order) log_rate[tree$edge[e, 2]] <- log_rate[tree$edge[e, 1]] + steps[e]
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
