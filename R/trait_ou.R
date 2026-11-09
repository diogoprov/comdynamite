#' Traits evolving by an Ornstein-Uhlenbeck process with one or several optima
#'
#' Simulates traits along a tree under the Ornstein-Uhlenbeck (OU) model of
#' Hansen (1997), with the multi-optimum ("multi-peak") extension of Butler &
#' King (2004). Along each branch the trait is pulled towards the optimum
#' `theta` of the selective regime of that branch with strength `alpha`, while
#' drifting with rate `sigma2`. Values are drawn exactly from the transition
#' distribution of the process: over a branch of length \eqn{t}, a trait with value
#' \eqn{x} at the start has, at the end, a normal distribution with mean
#' \eqn{\theta + (x - \theta) e^{-\alpha t}}{theta + (x - theta) * exp(-alpha * t)} and variance
#' \eqn{\frac{\sigma^2}{2 \alpha} (1 - e^{-2 \alpha t})}{sigma2 / (2 * alpha) * (1 - exp(-2 *
#' alpha * t))}, with \eqn{\sigma^2}{sigma2} = `sigma2`, \eqn{\alpha}{alpha} = `alpha` and
#' \eqn{\theta}{theta} = `theta`.
#'
#' With `alpha = 0` the process is Brownian motion with rate `sigma2`, and the
#' optima play no role. With regimes given as a vector, each branch belongs to a single
#' regime, so a shift in optimum happens at the start of the branch where the regime changes.
#' With regimes given as a mapped history (class `simmap`, from [trait_mk()] or phytools),
#' regimes can also change within branches, at the times of the map.
#'
#' Several traits evolve independently with `sigma2` as one number, or correlated with an
#' evolutionary rate matrix, whose diagonal holds the rate of each trait
#' and whose off-diagonal elements hold the evolutionary covariances. A named list of matrices
#' gives each regime its own matrix, as in the multi-rate models of mvMORPH (Clavel et al.
#' 2015); [rate_matrix()] builds them with the structures that mvMORPH fits. The strength of
#' selection `alpha` is the same for all traits; the optimum `theta` can differ between traits
#' (a matrix with one column per trait), so that regimes shift each trait by a different
#' amount. Under Brownian motion
#' (`alpha = 0`) the traits of the tips have covariance \eqn{\sum_k R_k \otimes C_k}{sum_k
#' kronecker(R_k, C_k)}, where \eqn{R_k} is the rate matrix of regime \eqn{k} and \eqn{C_k} the
#' phylogenetic covariance matrix built from the time spent in regime \eqn{k} (Clavel et al.
#' 2015).
#'
#' @param tree A rooted `phylo` tree with branch lengths.
#' @param n_traits Number of traits. Taken from the rate matrices when `sigma2` gives them, or
#'   else from the columns of `theta` when it is a matrix.
#' @param alpha Strength of selection towards the optimum (0 = Brownian
#'   motion).
#' @param sigma2 Rate of the stochastic (drift) component: one number for all traits, which
#'   then evolve independently; a symmetric positive semidefinite rate matrix with one row and
#'   column per trait, for correlated traits; or a named list or vector with a matrix or a
#'   number per regime (names as in `regimes`). See [rate_matrix()]. Regimes in the list but
#'   not in the tree are ignored.
#' @param theta Optimum of each regime: one unnamed value for a single optimum, or a
#'   named vector with one value per regime (names as in `regimes`). For optima that differ
#'   between traits, a matrix with one column per trait and one row per regime, with the
#'   regimes as row names (one unnamed row for a single optimum).
#' @param regimes Selective regime of each branch, in the order of the edges
#'   of `tree` reordered "cladewise" (`ape::reorder.phylo(tree, "cladewise")`),
#'   as returned by [paint_regimes()]; or a tree with the same edges as `tree` and a mapped
#'   history (class `simmap`, from [trait_mk()] or `phytools::make.simmap()`), whose states are
#'   the regimes, in any order of edges. The mapped tree can also be passed as `tree` itself.
#'   `NULL` puts every branch in one regime.
#' @param root Trait value at the root: one number for all traits, or one per trait.
#' @return Matrix with one row per tip (tip labels as row names) and one
#'   column per trait, named as the columns of the rate matrix, if any. The attribute `nodes` holds the values at the internal nodes, one row
#'   per node named by its number in `tree$edge`, the root first.
#' @references
#'   Butler, M. A., & King, A. A. (2004). Phylogenetic comparative analysis: a
#'   modeling approach for adaptive evolution. *The American Naturalist*, 164,
#'   683--695. \doi{10.1086/426002}
#'
#'   Clavel, J., Escarguel, G., & Merceron, G. (2015). mvMORPH: an R package for fitting
#'   multivariate evolutionary models to morphometric data. *Methods in Ecology and
#'   Evolution*, 6, 1311--1319. \doi{10.1111/2041-210X.12420}
#'
#'   Hansen, T. F. (1997). Stabilizing selection and the comparative analysis
#'   of adaptation. *Evolution*, 51, 1341--1351.
#'   \doi{10.1111/j.1558-5646.1997.tb01457.x}
#' @seealso [paint_regimes()] to assign regimes to clades, [rate_matrix()] for rate matrices,
#'   [trait_mk()] for a history of regimes that change within branches.
#' @family species traits
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(20)
#' # Brownian motion
#' trait_ou(tr, alpha = 0)[1:3, ]
#' # two optima: the clade of node 25 evolves towards 10, the rest towards 0
#' reg <- paint_regimes(tr, nodes = 25)
#' x <- trait_ou(tr, alpha = 2, theta = c(`0` = 0, `1` = 10), regimes = reg)
#' round(x[, 1], 1)
#' # two correlated traits
#' r <- rate_matrix(c(1, 0.5), correlation = 0.8)
#' y <- trait_ou(tr, sigma2 = r)
#' cor(y)
#' # the shift to regime 1 moves the optimum of trait 1 but not of trait 2
#' th <- rbind(`0` = c(0, 0), `1` = c(10, 0))
#' z <- trait_ou(tr, alpha = 2, sigma2 = r, theta = th, regimes = reg)
#' plot_traits(z, tr, regimes = reg)
#' @export
trait_ou <- function(tree, n_traits = 1, alpha = 0, sigma2 = 1, theta = 0,
                      regimes = NULL, root = 0) {
  check_tree(tree)
  check_count(n_traits, "n_traits")
  if (!is.numeric(alpha) || any(alpha < 0)) stop("`alpha` must be 0 or more", call. = FALSE)
  order <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  maps <- NULL
  if (inherits(regimes, "phylo")) {
    # match the edges by parent and child, whatever the order of the edges of each tree
    at <- match(paste(tree$edge[order, 1], tree$edge[order, 2]),
                paste(regimes$edge[, 1], regimes$edge[, 2]))
    if (is.null(regimes$maps) || nrow(regimes$edge) != nrow(tree$edge) || anyNA(at)) {
      stop("a tree given as `regimes` needs a mapped history (class `simmap`) and the same ",
           "edges as `tree`", call. = FALSE)
    }
    maps <- regimes$maps[at]
    regimes <- NULL
  }
  tree <- ape::reorder.phylo(tree, "cladewise")
  ne <- nrow(tree$edge)
  t <- tree$edge.length
  if (is.null(maps)) {
    if (is.null(regimes) && (if (is.matrix(theta)) nrow(theta) else length(theta)) != 1) {
      stop("several optima need `regimes`", call. = FALSE)
    }
    if (is.null(regimes)) {
      reg <- rep("0", ne)
    } else {
      stopifnot(length(regimes) == ne)
      reg <- as.character(regimes)
    }
    maps <- lapply(seq_len(ne), function(e) stats::setNames(t[e], reg[e]))
  }
  states <- unique(unlist(lapply(maps, names), use.names = FALSE))
  by_regime <- function(value, what) {
    if (length(value) == 1 && is.null(names(value)) && !is.list(value)) {
      return(stats::setNames(rep(list(value[[1]]), length(states)), states))
    }
    if (is.null(names(value)) || !all(states %in% names(value))) {
      stop("`", what, "` must be named with every regime: ", paste(states, collapse = ", "),
           call. = FALSE)
    }
    as.list(value)[states]
  }
  sigma2 <- strip_rate_matrix(sigma2)
  th <- if (is.matrix(theta)) {  # one row per regime, one column per trait
    if (nrow(theta) == 1 && is.null(rownames(theta))) {
      stats::setNames(rep(list(theta[1, ]), length(states)), states)
    } else {
      by_regime(stats::setNames(lapply(seq_len(nrow(theta)), function(i) theta[i, ]),
                                rownames(theta)), "theta")
    }
  } else {
    by_regime(theta, "theta")
  }
  rates <- if (is.matrix(sigma2)) {
    stats::setNames(rep(list(sigma2), length(states)), states)
  } else {
    by_regime(sigma2, "sigma2")
  }
  is_matrix <- vapply(rates, is.matrix, logical(1))
  multivariate <- any(is_matrix)
  if (multivariate) {
    if (!all(is_matrix)) {
      stop("`sigma2` must be numbers or rate matrices for all regimes, not both", call. = FALSE)
    }
    p <- unique(vapply(rates, ncol, integer(1)))
    if (length(p) != 1) stop("the rate matrices in `sigma2` must have the same size", call. = FALSE)
    if (n_traits != 1 && n_traits != p) {
      stop("`n_traits` is ", n_traits, " but the rate matrices have ", p, " traits",
           call. = FALSE)
    }
    n_traits <- p
    trait_names <- colnames(rates[[1]])
    rates <- lapply(rates, matrix_root)
  } else {
    trait_names <- NULL
    if (any(!vapply(rates, function(s) is.numeric(s) && length(s) == 1 && is.finite(s) &&
                      s >= 0, logical(1)))) {
      stop("`sigma2` must be zero or more", call. = FALSE)
    }
  }
  if (is.matrix(theta)) {
    if (n_traits == 1 && !multivariate) n_traits <- ncol(theta)
    if (ncol(theta) != n_traits) {
      stop("`theta` has ", ncol(theta), " columns but there are ", n_traits, " traits",
           call. = FALSE)
    }
  }
  if (!is.numeric(root) || anyNA(root) || !length(root) %in% c(1, n_traits)) {
    stop("`root` needs 1 or ", n_traits, " numbers", call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  x <- matrix(0, n_tip + tree$Nnode, n_traits)
  x[n_tip + 1, ] <- root
  for (e in seq_len(ne)) {  # cladewise: parents before children
    value <- x[tree$edge[e, 1], ]
    m <- maps[[e]]
    for (k in seq_along(m)) {
      s <- names(m)[k]
      len <- m[[k]]
      # same arithmetic as before rate matrices, so scalar rates give the same values
      rate <- if (multivariate) 1 else rates[[s]]
      if (alpha > 0) {
        decay <- exp(-alpha * len)
        sd_e <- sqrt(rate / (2 * alpha) * (1 - exp(-2 * alpha * len)))
      } else {
        decay <- 1
        sd_e <- sqrt(rate * len)
      }
      noise <- if (multivariate) {
        sd_e * drop(stats::rnorm(n_traits) %*% rates[[s]])
      } else {
        stats::rnorm(n_traits, 0, sd_e)
      }
      value <- th[[s]] + (value - th[[s]]) * decay + noise
    }
    x[tree$edge[e, 2], ] <- value
  }
  out <- x[seq_len(n_tip), , drop = FALSE]
  if (is.null(trait_names)) trait_names <- paste0("trait", seq_len(n_traits))
  dimnames(out) <- list(tree$tip.label, trait_names)
  attr(out, "nodes") <- x[-seq_len(n_tip), , drop = FALSE]
  dimnames(attr(out, "nodes")) <- list(n_tip + seq_len(tree$Nnode), colnames(out))
  out
}

# upper triangular U with crossprod(U) = m, for a symmetric positive semidefinite matrix
matrix_root <- function(m) {
  if (!is.numeric(m) || nrow(m) != ncol(m) || anyNA(m) || !isSymmetric(unname(m))) {
    stop("a rate matrix must be square, symmetric and without NA", call. = FALSE)
  }
  u <- tryCatch(chol(m), error = function(e) NULL)
  if (!is.null(u)) return(u)
  ev <- eigen(m, symmetric = TRUE)
  if (min(ev$values) < -sqrt(.Machine$double.eps) * max(1, abs(ev$values))) {
    stop("a rate matrix must be positive semidefinite", call. = FALSE)
  }
  t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(m)))
}

#' Assign selective regimes to clades
#'
#' Paints the branches of a tree with selective regimes for [trait_ou()]:
#' every branch starts in regime `"0"`, and the branches descending from each
#' node in `nodes` take a new regime (`"1"`, `"2"`, ...). Nodes are painted in
#' order, so a later node nested in an earlier one overrides it.
#'
#' @param tree A rooted `phylo` tree.
#' @param nodes Internal node numbers whose descendant branches get a new
#'   regime.
#' @param stem Should the branch leading to each node also take the new
#'   regime (the shift happens at the start of the stem rather than at the
#'   node)?
#' @return Character vector with the regime of each branch, in the order of
#'   `tree$edge` after reordering it "cladewise" (the order [trait_ou()]
#'   uses).
#' @family species traits
#' @examples
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(8)
#' table(paint_regimes(tr, nodes = 10))
#' @export
paint_regimes <- function(tree, nodes, stem = FALSE) {
  check_tree(tree, branch_lengths = FALSE)
  tree <- ape::reorder.phylo(tree, "cladewise")
  n_tip <- ape::Ntip(tree)
  stopifnot(all(nodes > n_tip), all(nodes <= n_tip + tree$Nnode))
  reg <- rep("0", nrow(tree$edge))
  for (k in seq_along(nodes)) {
    desc <- c(nodes[k], phangorn_descendants(tree, nodes[k]))
    in_clade <- tree$edge[, 1] %in% desc
    if (stem) in_clade <- in_clade | tree$edge[, 2] == nodes[k]
    reg[in_clade] <- as.character(k)
  }
  structure(reg, order = "cladewise")
}

# all nodes (tips and internal) below `node`
phangorn_descendants <- function(tree, node) {
  out <- integer(0)
  front <- node
  while (length(front)) {
    kids <- tree$edge[tree$edge[, 1] %in% front, 2]
    out <- c(out, kids)
    front <- kids
  }
  out
}
