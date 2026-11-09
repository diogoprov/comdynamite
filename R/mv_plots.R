#' Plot traits on the tree: phylomorphospace or traitgram
#'
#' Draws simulated traits with the tree that generated them. With two traits, a
#' phylomorphospace: the tips and internal nodes in the plane of the two traits, joined by the
#' branches of the tree. With one trait, a traitgram: the trait against time from the root.
#' Branches are coloured by selective regime, with the change of colour at the point of each
#' shift in a mapped history.
#'
#' The values of the internal nodes are those kept by [trait_ou()] in the attribute `nodes`,
#' so the branches show the simulated path between nodes, not a reconstruction. Without that
#' attribute, only the tips are drawn.
#'
#' @param x Traits of the tips from [trait_ou()]: a matrix with one row per tip, named by tip
#'   label.
#' @param tree The tree of `x`, a `phylo` object; or `NULL` to draw only the tips.
#' @param traits The two columns of `x` to draw (names or numbers), or one for a traitgram.
#' @param regimes Regimes of the branches, as in [trait_ou()]: a vector in the cladewise order
#'   of the edges (from [paint_regimes()]) or a tree with a mapped history (class `simmap`, from
#'   [trait_mk()]). `NULL` draws every branch in one colour.
#' @param colours Colours of the regimes, named by regime or in the order of the sorted
#'   regimes; one colour when there are no regimes.
#' @param tip_labels Should the tips be labelled?
#' @param legend Draw the legend of the regimes?
#' @param point_size Size of the tips.
#' @param linewidth Width of the branches.
#' @param labels,... Deprecated: `labels` is now `tip_labels`; graphical arguments of base
#'   graphics are no longer used (modify the ggplot object instead).
#' @return A ggplot object.
#' @seealso [trait_ou()], [rate_matrix()] and [plot.rate_matrix()].
#' @family species traits
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(2)
#' tr <- treeshapesim::sim_tree_shape(30)
#' habitat <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
#' r <- rate_matrix(c(1, 1), correlation = list(forest = 0.8, open = -0.8),
#'                  constraint = "variance")
#' x <- trait_ou(habitat, sigma2 = r, regimes = habitat)
#' plot_traits(x, habitat, regimes = habitat)
#' # traitgram of the first trait
#' plot_traits(x, habitat, traits = 1, regimes = habitat)
#' @export
plot_traits <- function(x, tree = NULL, traits = c(1, 2), regimes = NULL, colours = NULL,
                        tip_labels = FALSE, legend = TRUE, point_size = 1.8, linewidth = 0.7,
                        ..., labels = NULL) {
  need_ggplot("plot_traits()")
  tip_labels <- renamed_arg(tip_labels, labels, "tip_labels", "labels")
  legend <- legend_flag(legend)
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
  .data <- ggplot2::.data
  if (is.null(tree) || is.null(nodes)) {
    if (length(traits) == 1) stop("a traitgram needs the tree and the values of the nodes",
                                  call. = FALSE)
    tips <- data.frame(a = x[, traits[1]], b = x[, traits[2]], label = rownames(x))
    p <- ggplot2::ggplot(tips, ggplot2::aes(.data$a, .data$b)) +
      ggplot2::geom_point(colour = (colours %||% "#0072B2")[1], size = point_size)
    if (tip_labels) p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.5,
                                                vjust = -0.6)
    return(p + ggplot2::labs(x = lab[1], y = lab[2]) + ggplot2::theme_minimal())
  }
  check_tree(tree, branch_lengths = FALSE)
  n_tip <- ape::Ntip(tree)
  if (!all(tree$tip.label %in% rownames(x)) || nrow(nodes) != tree$Nnode) {
    stop("`x` does not match `tree`: give the tree that `trait_ou()` simulated on",
         call. = FALSE)
  }
  values <- rbind(x[tree$tip.label, , drop = FALSE],
                  nodes[as.character(n_tip + seq_len(tree$Nnode)), , drop = FALSE])
  maps <- branch_maps(tree, regimes)
  states <- sort(unique(unlist(lapply(maps, names), use.names = FALSE)))
  cols <- unlist(regime_colours(states, colours))
  edge <- tree$edge
  xy <- if (length(traits) == 2) {
    values[, traits, drop = FALSE]
  } else {
    cbind(ape::node.depth.edgelength(tree), values[, traits])
  }
  # one segment per piece of each branch, split where a mapped history changes regime
  pieces <- do.call(rbind, lapply(seq_len(nrow(edge)), function(e) {
    m <- maps[[e]]
    from <- xy[edge[e, 1], ]
    to <- xy[edge[e, 2], ]
    frac <- if (sum(m) > 0) cumsum(c(0, m)) / sum(m) else c(0, rep(1, length(m)))
    data.frame(x = from[1] + (to[1] - from[1]) * frac[-length(frac)],
               y = from[2] + (to[2] - from[2]) * frac[-length(frac)],
               xend = from[1] + (to[1] - from[1]) * frac[-1],
               yend = from[2] + (to[2] - from[2]) * frac[-1],
               regime = names(m), end = c(rep(FALSE, length(m) - 1), TRUE))
  }))
  pieces$edge <- rep(seq_len(nrow(edge)), lengths(maps))
  tip_edge <- match(seq_len(n_tip), edge[, 2])
  last <- pieces[pieces$end, ]
  tips <- data.frame(x = xy[seq_len(n_tip), 1], y = xy[seq_len(n_tip), 2],
                     regime = last$regime[match(tip_edge, last$edge)], label = tree$tip.label)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = pieces, ggplot2::aes(x = .data$x, y = .data$y,
                                                      xend = .data$xend, yend = .data$yend,
                                                      colour = .data$regime),
                          linewidth = linewidth) +
    ggplot2::geom_point(data = tips, ggplot2::aes(.data$x, .data$y, colour = .data$regime),
                        size = point_size)
  if (length(traits) == 2) {
    root <- data.frame(x = xy[n_tip + 1, 1], y = xy[n_tip + 1, 2])
    p <- p + ggplot2::geom_point(data = root, ggplot2::aes(.data$x, .data$y), shape = 15,
                                 size = point_size * 1.4)
  }
  if (tip_labels) {
    p <- p + ggplot2::geom_text(data = tips, ggplot2::aes(.data$x, .data$y, label = .data$label),
                                size = 2.5, hjust = if (length(traits) == 2) 0.5 else -0.15,
                                vjust = if (length(traits) == 2) -0.7 else 0.5)
  }
  p + ggplot2::scale_colour_manual(values = cols, name = "Regime",
                                   guide = if (legend && length(states) > 1) "legend" else "none") +
    ggplot2::labs(x = if (length(traits) == 2) lab[1] else "Time from the root",
                  y = if (length(traits) == 2) lab[2] else lab[1]) +
    ggplot2::theme_minimal()
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
#' @param legend Draw the legend of the regimes?
#' @param linewidth Width of the ellipses.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [rate_matrix()], [plot_traits()].
#' @family species traits
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' r <- rate_matrix(c(1, 2), correlation = 0.6, constraint = "shared",
#'                  eigenvalues = list(slow = c(1, 0.2), fast = c(1, 2)))
#' plot(r)
#' @export
plot.rate_matrix <- function(x, traits = c(1, 2), level = 1, colours = NULL, legend = TRUE,
                             linewidth = 0.9, ...) {
  need_ggplot("plot() of rate_matrix()")
  legend <- legend_flag(legend)
  ms <- if (is.list(x)) lapply(unclass(x), unclass) else list(rates = unclass(x))
  if (ncol(ms[[1]]) < 2) stop("a rate matrix of one trait has no ellipse", call. = FALSE)
  if (length(traits) != 2) stop("`traits` must name two traits", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || !(level > 0)) {
    stop("`level` must be one positive number", call. = FALSE)
  }
  cols <- unlist(regime_colours(names(ms), colours))
  ells <- do.call(rbind, lapply(names(ms), function(k) {
    e <- ellipse_points(ms[[k]][traits, traits, drop = FALSE], level)
    data.frame(x = e[, 1], y = e[, 2], regime = k)
  }))
  ells$regime <- factor(ells$regime, names(ms))
  lim <- max(abs(c(ells$x, ells$y)))
  lab <- colnames(ms[[1]])[traits] %||% paste("trait", traits)
  .data <- ggplot2::.data
  ggplot2::ggplot(ells, ggplot2::aes(.data$x, .data$y, colour = .data$regime)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_path(linewidth = linewidth) +
    ggplot2::scale_colour_manual(values = cols, name = "Regime",
                                 guide = if (legend && length(ms) > 1) "legend" else "none") +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = paste("Change in", lab[1]), y = paste("Change in", lab[2])) +
    ggplot2::theme_minimal()
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
#' @param species_labels Should the species be labelled?
#' @param point_size Size of the species.
#' @param labels,... Deprecated: `labels` is now `species_labels`; graphical arguments of base
#'   graphics are no longer used. The aspect ratio is always 1, so that the contours of
#'   Euclidean distances are circles.
#' @return A ggplot object.
#' @seealso [assemble_community()].
#' @family rule-based assembly
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(40)
#' y <- trait_ou(tr, sigma2 = rate_matrix(c(1, 4), correlation = 0.9))
#' plot_assembly(assemble_community(y, 15, rule = "filtering", optimum = c(0, 0))) +
#'   ggplot2::ggtitle("Euclidean")
#' plot_assembly(assemble_community(y, 15, rule = "filtering", optimum = c(0, 0),
#'                                  distance = "mahalanobis")) +
#'   ggplot2::ggtitle("Mahalanobis")
#' @export
plot_assembly <- function(x, traits = c(1, 2), colours = c("#0072B2", "grey60", "#D55E00"),
                          species_labels = FALSE, point_size = 2, ..., labels = NULL) {
  need_ggplot("plot_assembly()")
  species_labels <- renamed_arg(species_labels, labels, "species_labels", "labels")
  if (!inherits(x, "mc_assembly")) {
    stop("`x` must be the output of assemble_community()", call. = FALSE)
  }
  check_colours(colours, 3)
  tr <- x$traits
  if (ncol(tr) < 2) stop("the assembly used one trait: there is no plane to draw",
                         call. = FALSE)
  if (length(traits) != 2) stop("`traits` must name two traits", call. = FALSE)
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
  lab <- colnames(tr)[traits] %||% paste("trait", traits)
  species <- data.frame(x = tr[, traits[1]], y = tr[, traits[2]], label = rownames(tr),
                        status = factor(ifelse(kept, "Kept", "Removed"), c("Kept", "Removed")))
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(species, ggplot2::aes(.data$x, .data$y))
  if (length(contours) > 0) {
    lines <- do.call(rbind, lapply(seq_along(contours), function(k) {
      data.frame(x = contours[[k]][, 1], y = contours[[k]][, 2], group = k)
    }))
    p <- p + ggplot2::geom_path(data = lines, ggplot2::aes(group = .data$group),
                                colour = colours[3])
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$status, shape = .data$status),
                               size = point_size) +
    ggplot2::scale_colour_manual(values = c(Kept = colours[1], Removed = colours[2]),
                                 name = NULL, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(Kept = 19, Removed = 1), name = NULL, drop = FALSE)
  if (x$rule == "filtering") {
    optimum <- data.frame(x = x$optimum[traits[1]], y = x$optimum[traits[2]])
    p <- p + ggplot2::geom_point(data = optimum, shape = 3, size = point_size * 2, stroke = 1.2,
                                 colour = colours[3])
  }
  if (species_labels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.5, vjust = -0.7)
  }
  p + ggplot2::coord_equal() + ggplot2::labs(x = lab[1], y = lab[2]) +
    ggplot2::theme_minimal()
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

# stops with a message naming the function when ggplot2 (Suggests) is not installed
need_ggplot <- function(what) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(what, " needs the package ggplot2: install.packages(\"ggplot2\")", call. = FALSE) # nocov
  }
}

# `legend` as TRUE or FALSE; a position keyword or NULL, as in base graphics before 0.7.0,
# still works
legend_flag <- function(legend) {
  if (is.null(legend)) return(FALSE)
  if (is.character(legend)) return(TRUE)
  if (!is.logical(legend) || length(legend) != 1 || is.na(legend)) {
    stop("`legend` must be TRUE or FALSE", call. = FALSE)
  }
  legend
}
