#' Plot traits on the tree: phylomorphospace or traitgram
#'
#' Draws simulated traits with the tree that generated them. With two traits, a
#' phylomorphospace: the tips and internal nodes in the plane of the two traits, joined by the
#' branches of the tree. With one trait, a traitgram: the trait against time from the root.
#' Branches are coloured by selective regime, with the change of colour at the point of each
#' shift in a mapped history.
#'
#' The values of the internal nodes are those kept by [traits_ou()] in the attribute `nodes`,
#' so the branches show the simulated path between nodes, not a reconstruction. Without that
#' attribute, only the tips are drawn.
#'
#' @param x Traits of the tips from [traits_ou()]: a matrix with one row per tip, named by tip
#'   label.
#' @param tree The tree of `x`, a `phylo` object; or `NULL` to draw only the tips.
#' @param traits The two columns of `x` to draw (names or numbers), or one for a traitgram.
#' @param regimes Regimes of the branches, as in [traits_ou()]: a vector in the cladewise order
#'   of the edges (from [paint_regimes()]) or a tree with a mapped history (class `simmap`, from
#'   [trait_mk()]). `NULL` draws every branch in one colour.
#' @param colours Colours of the regimes, named by regime or in the order of the sorted
#'   regimes; one colour when there are no regimes.
#' @param labels Should the tips be labelled?
#' @param legend Where to put the legend of the regimes (a keyword of
#'   [graphics::legend()]), or `NULL` for none.
#' @param ... Further arguments to [graphics::plot()], such as `main`, `xlim` or `asp`.
#' @return `NULL`, invisibly; called for the plot.
#' @seealso [traits_ou()], [rate_matrix()] and [plot.rate_matrix()].
#' @family species traits
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(2)
#' tr <- treeshapesim::sim_tree_shape(30)
#' habitat <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
#' r <- rate_matrix(c(1, 1), correlation = list(forest = 0.8, open = -0.8),
#'                  constraint = "variance")
#' x <- traits_ou(habitat, sigma2 = r, regimes = habitat)
#' plot_traits(x, habitat, regimes = habitat)
#' # traitgram of the first trait
#' plot_traits(x, habitat, traits = 1, regimes = habitat)
#' @export
plot_traits <- function(x, tree = NULL, traits = c(1, 2), regimes = NULL, colours = NULL,
                        labels = FALSE, legend = "topleft", ...) {
  x <- as.matrix(x)
  if (!is.numeric(x) || is.null(rownames(x))) {
    stop("`x` must be a numeric matrix with tip labels as row names", call. = FALSE)
  }
  if (!length(traits) %in% 1:2) stop("`traits` must name one or two columns", call. = FALSE)
  if (length(traits) == 2 && ncol(x) < 2) {
    stop("`x` has one trait: use `traits = 1` for a traitgram", call. = FALSE)
  }
  if (is.character(traits) && !all(traits %in% colnames(x)) ||
      is.numeric(traits) && any(traits > ncol(x) | traits < 1)) {
    stop("`traits` must name columns of `x`", call. = FALSE)
  }
  lab <- colnames(x)[traits] %||% paste("trait", traits)
  nodes <- attr(x, "nodes")
  if (is.null(tree) || is.null(nodes)) {
    if (length(traits) == 1) stop("a traitgram needs the tree and the values of the nodes",
                                  call. = FALSE)
    cols <- colours %||% "#0072B2"
    graphics::plot(x[, traits[1]], x[, traits[2]], xlab = lab[1], ylab = lab[2], pch = 19,
                   col = cols[1], ...)
    if (labels) graphics::text(x[, traits[1]], x[, traits[2]], rownames(x), pos = 3, cex = 0.7)
    return(invisible(NULL))
  }
  check_tree(tree, branch_lengths = FALSE)
  n_tip <- ape::Ntip(tree)
  if (!all(tree$tip.label %in% rownames(x)) || nrow(nodes) != tree$Nnode) {
    stop("`x` does not match `tree`: give the tree that `traits_ou()` simulated on",
         call. = FALSE)
  }
  values <- rbind(x[tree$tip.label, , drop = FALSE], nodes[as.character(n_tip +
                                                                          seq_len(tree$Nnode)), ,
                                                          drop = FALSE])
  maps <- branch_maps(tree, regimes)
  states <- sort(unique(unlist(lapply(maps, names), use.names = FALSE)))
  cols <- regime_colours(states, colours)
  edge <- tree$edge
  # coordinates of the start and end of each piece of each branch
  if (length(traits) == 2) {
    xy <- values[, traits, drop = FALSE]
  } else {
    depth <- ape::node.depth.edgelength(tree)
    xy <- cbind(depth, values[, traits])
  }
  args <- list(...)
  plot_args <- utils::modifyList(list(x = range(xy[, 1]), y = range(xy[, 2]), type = "n",
                                      xlab = if (length(traits) == 2) lab[1] else "time",
                                      ylab = if (length(traits) == 2) lab[2] else lab[1]),
                                 args)
  do.call(graphics::plot, plot_args)
  end_state <- character(nrow(edge))
  for (e in seq_len(nrow(edge))) {
    m <- maps[[e]]
    from <- xy[edge[e, 1], ]
    to <- xy[edge[e, 2], ]
    frac <- if (sum(m) > 0) cumsum(c(0, m)) / sum(m) else c(0, rep(1, length(m)))
    for (k in seq_along(m)) {
      a <- from + (to - from) * frac[k]
      b <- from + (to - from) * frac[k + 1]
      graphics::segments(a[1], a[2], b[1], b[2], col = cols[[names(m)[k]]], lwd = 1.5)
    }
    end_state[e] <- names(m)[length(m)]
  }
  tip_edge <- match(seq_len(n_tip), edge[, 2])
  graphics::points(xy[seq_len(n_tip), 1], xy[seq_len(n_tip), 2], pch = 19,
                   col = unlist(cols[end_state[tip_edge]]))
  if (length(traits) == 2) {
    graphics::points(xy[n_tip + 1, 1], xy[n_tip + 1, 2], pch = 15, cex = 1.2)  # the root
  }
  if (labels) {
    graphics::text(xy[seq_len(n_tip), 1], xy[seq_len(n_tip), 2], tree$tip.label,
                   pos = if (length(traits) == 2) 3 else 4, cex = 0.7, xpd = NA)
  }
  if (!is.null(legend) && length(states) > 1) {
    graphics::legend(legend, legend = states, col = unlist(cols), lwd = 2, bty = "n",
                     title = "regime")
  }
  invisible(NULL)
}

#' Plot evolutionary rate matrices as ellipses
#'
#' Draws each rate matrix as the ellipse of the changes of two traits per unit time (one
#' standard deviation by default), one ellipse per regime. The shape shows the evolutionary
#' correlation, the size the rates: proportional matrices give nested ellipses of the same
#' shape, shared eigenvectors give ellipses with the same axes, equal variances give ellipses
#' inscribed in the same box.
#'
#' @param x A `rate_matrix` object from [rate_matrix()].
#' @param traits The two traits to draw (names or numbers); with more traits, the ellipse is
#'   that of the marginal covariance of the two.
#' @param level Number of standard deviations of the ellipse.
#' @param colours Colours of the regimes, in the order of the regimes.
#' @param legend Where to put the legend of the regimes, or `NULL` for none.
#' @param ... Further arguments to [graphics::plot()].
#' @return `NULL`, invisibly; called for the plot.
#' @seealso [rate_matrix()], [plot_traits()].
#' @family species traits
#' @examples
#' r <- rate_matrix(c(1, 2), correlation = 0.6, constraint = "shared",
#'                  eigenvalues = list(slow = c(1, 0.2), fast = c(1, 2)))
#' plot(r)
#' @export
plot.rate_matrix <- function(x, traits = c(1, 2), level = 1, colours = NULL,
                             legend = "topleft", ...) {
  ms <- if (is.list(x)) lapply(unclass(x), unclass) else list(rates = unclass(x))
  if (ncol(ms[[1]]) < 2) stop("a rate matrix of one trait has no ellipse", call. = FALSE)
  stopifnot(length(traits) == 2, is.numeric(level), length(level) == 1, level > 0)
  ells <- lapply(ms, function(m) ellipse_points(m[traits, traits, drop = FALSE], level))
  cols <- regime_colours(names(ms), colours)
  lim <- max(abs(unlist(ells)))
  lab <- colnames(ms[[1]])[traits] %||% paste("trait", traits)
  plot_args <- utils::modifyList(list(x = c(-lim, lim), y = c(-lim, lim), type = "n",
                                      asp = 1, xlab = paste("change in", lab[1]),
                                      ylab = paste("change in", lab[2])), list(...))
  do.call(graphics::plot, plot_args)
  graphics::abline(h = 0, v = 0, col = "grey85")
  for (k in seq_along(ells)) graphics::lines(ells[[k]], col = cols[[k]], lwd = 2)
  if (!is.null(legend) && length(ms) > 1) {
    graphics::legend(legend, legend = names(ms), col = unlist(cols), lwd = 2, bty = "n",
                     title = "regime")
  }
  invisible(NULL)
}

#' Plot a community assembled by rules in trait space
#'
#' Draws the regional pool in the plane of two traits, with the species kept by
#' [assemble_community()] filled and those removed open, in the metric of the assembly:
#'
#' * filtering: the optimum (a cross) and the boundary of the filter, the contour at the
#'   distance from the optimum of the farthest species kept (a circle for Euclidean distances,
#'   an ellipse for standardized or Mahalanobis distances);
#' * limiting similarity: around each species kept, the contour at half the smallest distance
#'   between species kept, so that contours touch but do not overlap;
#' * neutral: the species only.
#'
#' With more than two traits, the contours are the projections of the ellipsoids on the plane
#' of the two traits, and species that look inside a contour may lie outside it in the other
#' traits.
#'
#' @param x An `mc_assembly` object from [assemble_community()].
#' @param traits The two traits to draw (names or numbers).
#' @param colours Colours of the species kept and removed, and of the contours.
#' @param labels Should the species be labelled?
#' @param ... Further arguments to [graphics::plot()]. The aspect ratio is 1 by default
#'   (`asp = 1`), so that the contours of Euclidean distances are circles.
#' @return `NULL`, invisibly; called for the plot.
#' @seealso [assemble_community()].
#' @family rule-based assembly
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(40)
#' y <- traits_ou(tr, sigma2 = rate_matrix(c(1, 4), correlation = 0.9))
#' op <- par(mfrow = c(1, 2))
#' plot_assembly(assemble_community(y, 15, rule = "filtering", optimum = c(0, 0)),
#'               main = "Euclidean")
#' plot_assembly(assemble_community(y, 15, rule = "filtering", optimum = c(0, 0),
#'                                  distance = "mahalanobis"), main = "Mahalanobis")
#' par(op)
#' @export
plot_assembly <- function(x, traits = c(1, 2), colours = c("#0072B2", "grey60", "#D55E00"),
                          labels = FALSE, ...) {
  if (!inherits(x, "mc_assembly")) {
    stop("`x` must be the output of assemble_community()", call. = FALSE)
  }
  check_colours(colours, 3)
  tr <- x$traits
  if (ncol(tr) < 2) stop("the assembly used one trait: there is no plane to draw",
                         call. = FALSE)
  stopifnot(length(traits) == 2)
  kept <- rownames(tr) %in% x$community
  s <- x$metric[traits, traits, drop = FALSE]
  contours <- list()
  if (x$rule == "filtering") {
    w <- solve(chol(x$metric))
    d <- sqrt(rowSums((sweep(tr, 2, x$optimum) %*% w)^2))
    radius <- max(d[kept])
    contours <- list(sweep(ellipse_points(s, radius), 2, x$optimum[traits], "+"))
  } else if (x$rule == "limiting_similarity" && sum(kept) > 1) {
    w <- solve(chol(x$metric))
    radius <- min(stats::dist(tr[kept, , drop = FALSE] %*% w)) / 2
    unit <- ellipse_points(s, radius)
    contours <- lapply(which(kept), function(i) sweep(unit, 2, tr[i, traits], "+"))
  }
  pts <- tr[, traits, drop = FALSE]
  all_xy <- do.call(rbind, c(list(pts), contours))
  lab <- colnames(tr)[traits] %||% paste("trait", traits)
  plot_args <- utils::modifyList(list(x = range(all_xy[, 1]), y = range(all_xy[, 2]),
                                      type = "n", asp = 1, xlab = lab[1], ylab = lab[2]),
                                 list(...))
  do.call(graphics::plot, plot_args)
  for (cn in contours) graphics::lines(cn, col = colours[3])
  graphics::points(pts[!kept, 1], pts[!kept, 2], pch = 1, col = colours[2])
  graphics::points(pts[kept, 1], pts[kept, 2], pch = 19, col = colours[1])
  if (x$rule == "filtering") {
    graphics::points(x$optimum[traits[1]], x$optimum[traits[2]], pch = 3, cex = 1.5, lwd = 2,
                     col = colours[3])
  }
  if (labels) graphics::text(pts[, 1], pts[, 2], rownames(pts), pos = 3, cex = 0.7)
  invisible(NULL)
}

# points on the ellipse {d : d' S^-1 d = r^2}
ellipse_points <- function(s, r, n = 100) {
  ev <- eigen(s, symmetric = TRUE)
  a <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(a), sin(a))
  r * circle %*% t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), 2))
}

# the regimes along each edge of `tree` (in the order of tree$edge) as named lengths
branch_maps <- function(tree, regimes) {
  if (is.null(regimes)) {
    return(lapply(tree$edge.length, function(len) stats::setNames(len, "all")))
  }
  if (inherits(regimes, "phylo")) {
    at <- match(paste(tree$edge[, 1], tree$edge[, 2]),
                paste(regimes$edge[, 1], regimes$edge[, 2]))
    if (is.null(regimes$maps) || anyNA(at)) {
      stop("a tree given as `regimes` needs a mapped history (class `simmap`) and the same ",
           "edges as `tree`", call. = FALSE)
    }
    return(regimes$maps[at])
  }
  if (length(regimes) != nrow(tree$edge)) {
    stop("`regimes` needs one value per branch, in cladewise order", call. = FALSE)
  }
  order <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  reg <- character(nrow(tree$edge))
  reg[order] <- as.character(regimes)
  lapply(seq_len(nrow(tree$edge)), function(e) stats::setNames(tree$edge.length[e], reg[e]))
}

# a named list of colours, one per state
regime_colours <- function(states, colours) {
  okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442")
  if (is.null(colours)) {
    if (length(states) > length(okabe_ito)) {
      colours <- grDevices::hcl.colors(length(states), "Dark 3")
    } else {
      colours <- okabe_ito[seq_along(states)]
    }
  }
  if (!is.null(names(colours))) {
    if (!all(states %in% names(colours))) {
      stop("`colours` must be named with every regime: ", paste(states, collapse = ", "),
           call. = FALSE)
    }
    colours <- colours[states]
  } else {
    if (length(colours) < length(states)) {
      stop("`colours` needs ", length(states), " colours", call. = FALSE)
    }
    colours <- colours[seq_along(states)]
  }
  check_colours(unname(colours), length(states))
  stats::setNames(as.list(unname(colours)), states)
}
