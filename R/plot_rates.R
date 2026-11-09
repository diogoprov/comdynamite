# Rates of trait evolution painted on the branches of a tree, with ggplot2 (Suggests).

#' Plot rates that vary among branches
#'
#' Draws the tree as a rectangular phylogram with each branch coloured by its rate. Works with
#' any rate per branch: the truth from [trait_rates()], [history_rates()], [turnover_rates()]
#' or [food_web_rates()], or estimates such as those of RRphylo, as long as they follow the
#' order of the rows of `tree$edge`.
#'
#' @param tree A `phylo` object with branch lengths, or a result of [trait_rates()], which
#'   carries its tree and rates.
#' @param rates Rate of each branch, in the order of the rows of `tree$edge`, or a result of
#'   [trait_rates()], whose `edge_rate` is used. Left out when `tree` is a result of
#'   [trait_rates()].
#' @param log Colour by the logarithm of the rate (the default), since rates are multiplicative.
#'   Needs positive rates; use `FALSE` for estimates that can be negative.
#' @param palette Name of a palette of [grDevices::hcl.colors()].
#' @param tip_labels Draw the names of the tips.
#' @param linewidth Width of the branches.
#' @param title Title of the plot, or `NULL`.
#' @return A ggplot object.
#' @seealso [trait_rates()], [history_rates()], [turnover_rates()].
#' @family species traits
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(30)
#' x <- trait_rates(tree, "shift", shifts = 2, shift_size = c(8, 1 / 8))
#' plot_rates(x)
#' habitat <- trait_mk(tree, q = 0.5, states = c("forest", "open"))
#' plot_rates(tree, history_rates(habitat, "state", state_rates = c(forest = 1, open = 5)))
#' @export
plot_rates <- function(tree, rates = NULL, log = TRUE, palette = "viridis", tip_labels = FALSE,
                       linewidth = 1, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_rates() needs the package ggplot2", call. = FALSE) # nocov
  }
  if (!inherits(tree, "phylo") && is.list(tree) && inherits(tree$tree, "phylo") &&
      !is.null(tree$edge_rate)) {
    if (is.null(rates)) rates <- tree$edge_rate
    tree <- tree$tree
  }
  check_tree(tree)
  if (is.null(rates)) stop("`rates` is missing", call. = FALSE)
  if (is.list(rates) && !is.null(rates$edge_rate)) rates <- rates$edge_rate
  if (!is.numeric(rates) || length(rates) != nrow(tree$edge)) {
    stop("`rates` must have one value per branch (", nrow(tree$edge), "), in the order of ",
         "the rows of `tree$edge`", call. = FALSE)
  }
  stopifnot(is.logical(log), length(log) == 1, !is.na(log),
            is.logical(tip_labels), length(tip_labels) == 1, !is.na(tip_labels),
            is.character(palette), length(palette) == 1)
  rates <- as.numeric(rates)
  if (log && any(rates <= 0, na.rm = TRUE)) {
    stop("`log = TRUE` needs positive rates; use `log = FALSE`", call. = FALSE)
  }
  layout <- tree_layout(tree)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  branches <- data.frame(x = layout$x[parent], xend = layout$x[child], y = layout$y[child],
                         yend = layout$y[child], rate = rates)
  connectors <- data.frame(x = layout$x[parent], xend = layout$x[parent], y = layout$y[parent],
                           yend = layout$y[child], rate = rates)
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(rbind(connectors, branches),
                       ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                                    yend = .data$yend, colour = .data$rate)) +
    ggplot2::geom_segment(linewidth = linewidth, lineend = "square") +
    ggplot2::scale_colour_gradientn(colours = grDevices::hcl.colors(30, palette),
                                    transform = if (log) "log10" else "identity",
                                    labels = function(b) format(signif(b, 2), trim = TRUE),
                                    name = "Rate") +
    ggplot2::labs(x = "Time from the root", y = NULL, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line())
  if (tip_labels) {
    n <- length(tree$tip.label)
    tips <- data.frame(x = layout$x[seq_len(n)], y = layout$y[seq_len(n)],
                       label = tree$tip.label)
    p <- p + ggplot2::geom_text(data = tips, ggplot2::aes(x = .data$x, y = .data$y,
                                                          label = .data$label),
                                inherit.aes = FALSE, hjust = -0.1, size = 3) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.15)))
  }
  p
}

# coordinates of the nodes of a rectangular phylogram: x is the distance from the root, y puts
# the tips at 1..n in the order of a cladewise traversal and each internal node at the middle
# of its children
tree_layout <- function(tree) {
  n <- length(tree$tip.label)
  x <- ape::node.depth.edgelength(tree)
  y <- numeric(n + tree$Nnode)
  cladewise <- ape::reorder.phylo(tree, "cladewise")$edge[, 2]
  tips <- cladewise[cladewise <= n]
  y[tips] <- seq_along(tips)
  # internal nodes from the smallest clades up, so children come before parents
  internal <- n + seq_len(tree$Nnode)
  size <- ape::node.depth(tree)[internal]
  for (node in internal[order(size)]) {
    y[node] <- mean(range(y[tree$edge[tree$edge[, 1] == node, 2]]))
  }
  list(x = x, y = y)
}
