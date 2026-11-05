#' Traits evolving by an Ornstein-Uhlenbeck process with one or several optima
#'
#' Simulates traits along a tree under the Ornstein-Uhlenbeck (OU) model of
#' Hansen (1997), with the multi-optimum ("multi-peak") extension of Butler &
#' King (2004). Along each branch the trait is pulled towards the optimum
#' `theta` of the selective regime of that branch with strength `alpha`, while
#' drifting with rate `sigma2`. Values are drawn exactly from the transition
#' distribution of the process: over a branch of length t, a trait with value
#' x at the start has, at the end, a normal distribution with mean
#' `theta + (x - theta) * exp(-alpha * t)` and variance
#' `sigma2 / (2 * alpha) * (1 - exp(-2 * alpha * t))`.
#'
#' With `alpha = 0` the process is Brownian motion with rate `sigma2`, and the
#' optima play no role. Each branch belongs to a single regime, so a shift in
#' optimum happens at the start of the branch where the regime changes.
#'
#' @param tree A rooted `phylo` tree with branch lengths.
#' @param n_traits Number of independent traits.
#' @param alpha Strength of selection towards the optimum (0 = Brownian
#'   motion).
#' @param sigma2 Rate of the stochastic (drift) component.
#' @param theta Optimum of each regime: one value for a single optimum, or a
#'   named vector with one value per regime (names as in `regimes`).
#' @param regimes Selective regime of each branch, in the order of the edges
#'   of `tree` reordered "cladewise" (`ape::reorder.phylo(tree, "cladewise")`),
#'   as returned by [paint_regimes()]. `NULL` puts every branch in one
#'   regime.
#' @param root Trait value at the root.
#' @return Matrix with one row per tip (tip labels as row names) and one
#'   column per trait.
#' @references
#'   Butler, M. A., & King, A. A. (2004). Phylogenetic comparative analysis: a
#'   modeling approach for adaptive evolution. *The American Naturalist*, 164,
#'   683--695. \doi{10.1086/426002}
#'
#'   Hansen, T. F. (1997). Stabilizing selection and the comparative analysis
#'   of adaptation. *Evolution*, 51, 1341--1351.
#'   \doi{10.1111/j.1558-5646.1997.tb01457.x}
#' @seealso [paint_regimes()] to assign regimes to clades.
#' @examples
#' set.seed(1)
#' tr <- ape::rcoal(20)
#' # Brownian motion
#' traits_ou(tr, alpha = 0)[1:3, ]
#' # two optima: the clade of node 25 evolves towards 10, the rest towards 0
#' reg <- paint_regimes(tr, nodes = 25)
#' x <- traits_ou(tr, alpha = 2, theta = c(`0` = 0, `1` = 10), regimes = reg)
#' round(x[, 1], 1)
#' @export
traits_ou <- function(tree, n_traits = 1, alpha = 0, sigma2 = 1, theta = 0,
                      regimes = NULL, root = 0) {
  stopifnot(inherits(tree, "phylo"), !is.null(tree$edge.length), n_traits >= 1,
            alpha >= 0, sigma2 >= 0)
  tree <- ape::reorder.phylo(tree, "cladewise")
  ne <- nrow(tree$edge)
  if (is.null(regimes)) {
    if (length(theta) != 1) stop("several optima need `regimes`")
    th <- rep(unname(theta), ne)
  } else {
    stopifnot(length(regimes) == ne)
    reg <- as.character(regimes)
    if (length(theta) == 1 && is.null(names(theta))) {
      th <- rep(unname(theta), ne)
    } else {
      if (is.null(names(theta)) || !all(reg %in% names(theta)))
        stop("`theta` must be named with every regime in `regimes`")
      th <- unname(theta[reg])
    }
  }
  t <- tree$edge.length
  if (alpha > 0) {
    decay <- exp(-alpha * t)
    sd_e <- sqrt(sigma2 / (2 * alpha) * (1 - exp(-2 * alpha * t)))
  } else {
    decay <- rep(1, ne)
    sd_e <- sqrt(sigma2 * t)
  }
  n_tip <- ape::Ntip(tree)
  x <- matrix(0, n_tip + tree$Nnode, n_traits)
  x[n_tip + 1, ] <- root
  for (e in seq_len(ne)) {  # cladewise: parents before children
    pa <- tree$edge[e, 1]
    ch <- tree$edge[e, 2]
    x[ch, ] <- th[e] + (x[pa, ] - th[e]) * decay[e] + stats::rnorm(n_traits, 0, sd_e[e])
  }
  out <- x[seq_len(n_tip), , drop = FALSE]
  dimnames(out) <- list(tree$tip.label, paste0("trait", seq_len(n_traits)))
  out
}

#' Assign selective regimes to clades
#'
#' Paints the branches of a tree with selective regimes for [traits_ou()]:
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
#'   `tree$edge` after reordering it "cladewise" (the order [traits_ou()]
#'   uses).
#' @examples
#' tr <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
#' paint_regimes(tr, nodes = 6)
#' @export
paint_regimes <- function(tree, nodes, stem = FALSE) {
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
