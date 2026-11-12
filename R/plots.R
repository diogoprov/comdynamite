#' Plot the results of a node test on the phylogeny
#'
#' Draws the tree with the associations found by [node_test()] at each node, as in
#' Figure 3 of Leibold et al. (2010). Blue circles are predictors and red triangles are
#' biogeographic events, each with the names of the variables. Solid symbols remain
#' significant after the false discovery rate, and open symbols are significant only as
#' isolated tests.
#'
#' The article does not say which tests form a family for the correction. Here, the
#' correction of Benjamini & Hochberg (1995) uses all nodes and variables of each set, so
#' predictors and events are corrected separately.
#'
#' As the smallest p-value is 1 / (`n_perm` + 1), solid symbols need many permutations.
#' With few species, `test = "max"` may leave no symbols at all (see [node_test()]).
#'
#' @param x Result of [node_test()].
#' @param tree The `phylo` object that gave the contrasts of [node_contrasts()].
#' @param test Permutation test to show: `"max"` for the larger p-value of sites and
#'   species, `"site"` for sites alone, or `"joint"` for the permutation of the article.
#'   See [node_test()] for what each one tests.
#' @param type `"phylogram"`, as in the article, or `"fan"`, for very large trees, with the
#'   names of the variables along the radii so that they stay clear of the species.
#' @param level Significance level.
#' @param alpha Deprecated, use `level`.
#' @param colours Colours of predictors and biogeographic events.
#' @param text_size Size of the text, in mm (as in [ggplot2::geom_text()]).
#' @param legend Draw the legend?
#' @param col,cex,... Deprecated: `col` is now `colours` and `cex` is `text_size` (in mm,
#'   about 3.5 times the old value); arguments of [ape::plot.phylo()] are no longer used.
#' @return A ggplot object, with the table of symbols in the attribute `marks`: one row per
#'   symbol, with `node`, `set`, `variable`, `p`, `p_adjusted` and `corrected`, which is `TRUE`
#'   for solid symbols.
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing. *Journal of the Royal Statistical
#' Society: Series B*, 57, 289--300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @family clades and biogeography
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' nodes <- node_contrasts(daphniids$tree_leibold)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
#'                  n_perm = 49)
#' plot_node_test(out, daphniids$tree_leibold, test = "site")
#' \donttest{
#' # solid symbols need more permutations
#' out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
#'                  n_perm = 999)
#' plot_node_test(out, daphniids$tree_leibold, test = "site")
#' }
#' @export
plot_node_test <- function(x, tree, test = c("max", "site", "joint"),
                           type = c("phylogram", "fan"), level = 0.05,
                           colours = c("#0072B2", "#D55E00"), text_size = 2.8, legend = TRUE,
                           ..., alpha = NULL, col = NULL, cex = NULL) {
  need_ggplot("plot_node_test()")
  level <- deprecated_level(level, alpha)
  colours <- renamed_arg(colours, col, "colours", "col")
  if (!is.null(cex)) text_size <- 3.5 * renamed_arg(text_size, cex, "text_size", "cex")
  legend <- legend_flag(legend)
  if (!inherits(x, "mc_node_test")) stop("`x` must come from node_test()", call. = FALSE)
  check_tree(tree, branch_lengths = FALSE)
  check_colours(colours, 2)
  test <- match.arg(test)
  type <- match.arg(type)
  check_node_tree(x, tree)
  marks <- node_marks(x, test, level)
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  layout <- tree_layout(tree)
  n_tip <- length(tree$tip.label)
  depth <- max(layout$x)
  edges <- tree_segments(tree, layout)
  tips <- data.frame(x = layout$x[seq_len(n_tip)], y = layout$y[seq_len(n_tip)],
                     label = tree$tip.label)
  fan <- type == "fan"
  # text along the radius of a fan, turned on its left half so that it is not upside down
  turn <- function(y) {
    angle <- 90 - 360 * (y - 0.5) / n_tip
    flip <- angle < -90
    list(angle = ifelse(flip, angle + 180, angle), flip = flip)
  }
  .data <- ggplot2::.data
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = .data$x, xend = .data$xend,
                                                     y = .data$y, yend = .data$yend),
                          linewidth = 0.4)
  if (fan) {
    t <- turn(tips$y)
    tips$angle <- t$angle
    tips$hjust <- ifelse(t$flip, 1, 0)
    p <- p + ggplot2::geom_text(data = tips, ggplot2::aes(.data$x, .data$y, label = .data$label,
                                                          angle = .data$angle,
                                                          hjust = .data$hjust),
                                size = text_size, nudge_x = 0.02 * depth)
  } else {
    p <- p + ggplot2::geom_text(data = tips, ggplot2::aes(.data$x, .data$y, label = .data$label),
                                size = text_size, hjust = 0, nudge_x = 0.02 * depth)
  }
  if (!is.null(marks) && nrow(marks) > 0) {
    groups <- split(marks, list(marks$node, marks$set == "biogeography", !marks$corrected),
                    drop = TRUE)
    symbols <- do.call(rbind, lapply(groups, function(g) {
      data.frame(node = as.integer(g$node[1]), set = g$set[1], corrected = g$corrected[1],
                 label = paste(g$variable, collapse = ", "))
    }))
    # one symbol per set and level at each node: predictors above, events below, solid ones
    # next to the node
    symbols <- symbols[order(symbols$node, symbols$set, !symbols$corrected), ]
    side <- ifelse(symbols$set == "biogeography", -1, 1)
    rank <- stats::ave(side, symbols$node, side, FUN = seq_along)
    symbols$x <- layout$x[symbols$node]
    symbols$y <- layout$y[symbols$node] + side * (rank - 0.4) * 0.45
    symbols$key <- factor(paste(ifelse(symbols$set == "biogeography", "Biogeography",
                                       "Predictor"),
                                ifelse(symbols$corrected, "after correction", "uncorrected"),
                                sep = ", "),
                          c("Predictor, after correction", "Predictor, uncorrected",
                            "Biogeography, after correction", "Biogeography, uncorrected"))
    symbols$colour <- colours[1 + (symbols$set == "biogeography")]
    key_colour <- stats::setNames(rep(colours, each = 2), levels(symbols$key))
    key_fill <- stats::setNames(c(colours[1], "white", colours[2], "white"), levels(symbols$key))
    key_shape <- stats::setNames(c(21, 21, 24, 24), levels(symbols$key))
    if (fan) {
      t <- turn(symbols$y)
      symbols$angle <- t$angle
      symbols$hjust <- ifelse(t$flip, 0, 1)
    } else {
      symbols$angle <- 0
      symbols$hjust <- 1
    }
    p <- p +
      ggplot2::geom_label(data = symbols,
                          ggplot2::aes(.data$x, .data$y, label = .data$label,
                                       angle = .data$angle, hjust = .data$hjust),
                          colour = symbols$colour, fill = "white", linewidth = 0,
                          size = text_size * 0.9, nudge_x = -0.015 * depth,
                          label.padding = ggplot2::unit(0.05, "lines")) +
      ggplot2::geom_point(data = symbols, ggplot2::aes(.data$x, .data$y, colour = .data$key,
                                                       fill = .data$key, shape = .data$key),
                          size = 2.2, stroke = 0.9) +
      ggplot2::scale_colour_manual(values = key_colour, name = NULL, drop = TRUE) +
      ggplot2::scale_fill_manual(values = key_fill, name = NULL, drop = TRUE) +
      ggplot2::scale_shape_manual(values = key_shape, name = NULL, drop = TRUE)
  }
  p <- p + ggplot2::theme_void() +
    ggplot2::theme(legend.position = if (legend) "bottom" else "none")
  p <- if (fan) {
    p + ggplot2::coord_polar(theta = "y", clip = "off") +
      ggplot2::scale_y_continuous(limits = c(0.5, n_tip + 0.5), expand = c(0, 0)) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.2, 0.35)))
  } else {
    p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.25, 0.3))) +
      ggplot2::coord_cartesian(clip = "off")
  }
  attr(p, "marks") <- marks
  p
}

#' Plot the semipartial correlations of a node test
#'
#' Draws a table of nodes by variables with the semipartial correlations of [node_test()],
#' as in Table S1 of Leibold et al. (2010). Bold values with an asterisk remain significant
#' after the false discovery rate, as in [plot_node_test()]; other bold values are
#' significant only as isolated tests. Blank cells could not be calculated.
#'
#' Rows are the node numbers of the tree, as in `tree$edge`; `ape::nodelabels()` on a plot of
#' [ape::plot.phylo()] shows where they are.
#'
#' @inheritParams plot_node_test
#' @param colours Colours of negative and positive values.
#' @param col,cex,... Deprecated: `col` is now `colours` and `cex` is `text_size`; arguments of
#'   [graphics::image()] are no longer used.
#' @return A ggplot object, with the matrix of semipartial correlations (predictors followed
#'   by biogeographic events) in the attribute `values`.
#' @references
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @family clades and biogeography
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' nodes <- node_contrasts(daphniids$tree_leibold)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
#'                  n_perm = 49)
#' plot_node_heatmap(out, test = "site")
#' @export
plot_node_heatmap <- function(x, test = c("max", "site", "joint"), level = 0.05,
                              colours = c("#0072B2", "#D55E00"), text_size = 2.5, ...,
                              alpha = NULL, col = NULL, cex = NULL) {
  need_ggplot("plot_node_heatmap()")
  level <- deprecated_level(level, alpha)
  colours <- renamed_arg(colours, col, "colours", "col")
  if (!is.null(cex)) text_size <- 3.5 * renamed_arg(text_size, cex, "text_size", "cex")
  if (!inherits(x, "mc_node_test")) stop("`x` must come from node_test()", call. = FALSE)
  check_colours(colours, 2)
  test <- match.arg(test)
  values <- cbind(x$semipartial, x$biogeography$semipartial)
  if (all(is.na(values))) stop("no node was tested, so there is nothing to show", call. = FALSE)
  marks <- node_marks(x, test, level)
  significant <- matrix(0, nrow(values), ncol(values))
  significant[cbind(match(marks$node, rownames(values)), marks$column)] <- 1 + marks$corrected
  cells <- data.frame(node = factor(rep(rownames(values), ncol(values)), rev(rownames(values))),
                      # by position: a predictor and an event may share a name
                      variable = rep(seq_len(ncol(values)), each = nrow(values)),
                      value = as.vector(values), strength = as.vector(significant))
  cells$label <- ifelse(is.na(cells$value), "",
                        paste0(format(round(cells$value, 2), nsmall = 2),
                               c("", "", "*")[cells$strength + 1]))
  limit <- max(abs(values), na.rm = TRUE)
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(cells, ggplot2::aes(.data$variable, .data$node)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label,
                                    fontface = ifelse(.data$strength > 0, "bold", "plain")),
                       size = text_size) +
    ggplot2::scale_fill_gradient2(low = colours[1], mid = "white", high = colours[2],
                                  limits = c(-limit, limit), na.value = "white",
                                  name = "Semipartial\ncorrelation") +
    ggplot2::scale_x_continuous(breaks = seq_len(ncol(values)), labels = colnames(values),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Node") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid = ggplot2::element_blank())
  attr(p, "values") <- values
  p
}

# the tree must be able to hold the species of each node
check_node_tree <- function(x, tree) {
  sizes <- colSums(node_contrasts(tree) != 0)
  size <- sizes[as.character(x$nodes$node)]
  if (anyNA(size) || any(size < x$nodes$n_species)) {
    stop("`tree` is not the tree of the contrasts in `x`", call. = FALSE)
  }
}

# one row for each significant test, in the sets of predictors and events
node_marks <- function(x, test, alpha) {
  field <- c(max = "p_partial_max", site = "p_partial_site", joint = "p_partial_joint")[[test]]
  one_set <- function(p, set, offset = 0) {
    if (is.null(p)) return(NULL)
    table <- data.frame(node = rep(rownames(p), ncol(p)), set = set,
                        variable = rep(colnames(p), each = nrow(p)),
                        column = rep(seq_len(ncol(p)), each = nrow(p)) + offset, p = as.vector(p),
                        p_adjusted = stats::p.adjust(as.vector(p), "BH"))
    table <- table[!is.na(table$p) & table$p <= alpha, ]
    table$corrected <- table$p_adjusted <= alpha
    table
  }
  marks <- rbind(one_set(x[[field]], "environment"),
                 one_set(x$biogeography[[field]], "biogeography", ncol(x[[field]])))
  rownames(marks) <- NULL
  marks
}

#' Plot the fractions of a trait partition
#'
#' Draws the total variation as a line divided into the fractions of [trait_partition()],
#' with one bar for the environment and one for space, as in Figure 1 of Desdevises et al.
#' (2003). The two bars overlap in `b`, the variation they share. `a` is explained by the
#' environment alone, `c` by space alone, and `d` is unexplained. Dark segments are
#' significant, and `b` is grey when it was not tested (see `shared_null` in
#' [trait_partition()]).
#'
#' Negative fractions are shown as zero, following Legendre (2008). The result of
#' [trait_partition()] keeps them, because all fractions sum to one only with them.
#'
#' @param x Result of [trait_partition()].
#' @param value Column of the fractions to show: `"r2_adjusted"`, `"r2"`, or
#'   `"r2_component"` for the share of the mean or variance component.
#' @param test Permutation test: `"max"` for the larger p-value of sites and species,
#'   `"site"` or `"species"`.
#' @param level Significance level.
#' @param alpha Deprecated, use `level`.
#' @param colours Colour of significant segments.
#' @param labels Names of the two sets of predictors.
#' @param text_size Size of the text, in mm.
#' @param col,cex,... Deprecated: `col` is now `colours` and `cex` is `text_size`; arguments of
#'   [graphics::plot.window()] are no longer used.
#' @return A ggplot object, with the rows of `x$fractions` for `a`, `b`, `c` and `d` in the
#'   attribute `fractions`, with `shown`, the values drawn.
#' @references
#' Desdevises, Y., Legendre, P., Azouzi, L., & Morand, S. (2003). Quantifying
#' phylogenetically structured environmental variation. *Evolution*, 57, 2647--2652.
#' \doi{10.1111/j.0014-3820.2003.tb01508.x}
#'
#' Legendre, P. (2008). Studying beta diversity: ecological variation partitioning by
#' multiple regression and canonical analysis. *Journal of Plant Ecology*, 1, 3--8.
#' \doi{10.1093/jpe/rtm001}
#' @family environment and space
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' trait <- stats::setNames(seq(0, 100, 10), letters[1:11])
#' env <- runif(40, 0, 100)
#' space <- cbind(runif(40), runif(40))
#' comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
#' colnames(comm) <- names(trait)
#' out <- trait_partition(comm, trait, env, space, n_perm = 49)
#' plot_trait_partition(out)
#' @export
plot_trait_partition <- function(x, value = c("r2_adjusted", "r2", "r2_component"),
                                 test = c("max", "site", "species"), level = 0.05,
                                 colours = "#0072B2", labels = c("Environment", "Space"),
                                 text_size = 3.2, ..., alpha = NULL, col = NULL, cex = NULL) {
  need_ggplot("plot_trait_partition()")
  level <- deprecated_level(level, alpha)
  colours <- renamed_arg(colours, col, "colours", "col")
  if (!is.null(cex)) text_size <- 3.5 * renamed_arg(text_size, cex, "text_size", "cex")
  if (!inherits(x, "mc_trait_partition")) {
    stop("`x` must come from trait_partition()", call. = FALSE)
  }
  if (length(labels) != 2) stop("`labels` must give two names", call. = FALSE)
  check_colours(colours, 1)
  value <- match.arg(value)
  if (identical(test, "taxa")) {
    warning("`test = \"taxa\"` is deprecated, use `\"species\"`", call. = FALSE)
    test <- "species"
  }
  test <- match.arg(test)
  table <- x$fractions
  shown <- table[match(c("a", "b", "c", "d"), table$fraction), ]
  shown$shown <- pmax(shown[[value]], 0)
  p_of <- stats::setNames(shown[[paste0("p_", test)]], shown$fraction)
  status <- function(fraction) {
    if (is.na(p_of[[fraction]])) "Not tested"
    else if (p_of[[fraction]] <= level) "Significant" else "Not significant"
  }
  edges <- c(0, cumsum(shown$shown))
  total <- edges[5]
  bars <- data.frame(xmin = edges[c(1, 2, 2, 3)], xmax = edges[c(2, 3, 3, 4)],
                     ymin = c(2, 2, 1.1, 1.1),
                     status = vapply(c("a", "b", "b", "c"), status, character(1)))
  bars$ymax <- bars$ymin + 0.7
  bars$status <- factor(bars$status, c("Significant", "Not significant", "Not tested"))
  fills <- c(Significant = colours, `Not significant` = grDevices::adjustcolor(colours, 0.3),
             `Not tested` = "grey70")
  # beside the bar of space, on the side with more room
  on_right <- edges[5] - edges[4] >= edges[2]
  gap <- 0.01 * total
  texts <- data.frame(x = c(edges[1], if (on_right) edges[4] + gap else edges[2] - gap,
                            mean(edges[4:5])),
                      y = c(2.85, 1.45, 2.35), label = c(labels, "Unexplained"),
                      hjust = c(0, 1 - on_right, 0.5))
  fractions <- data.frame(x = fraction_positions(edges, 0.12 * total), y = 0.45,
                          label = paste0("[", shown$fraction, "] ",
                                         format(round(shown$shown, 2), nsmall = 2)))
  .data <- ggplot2::.data
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                                 ymin = .data$ymin, ymax = .data$ymax,
                                                 fill = .data$status), colour = "white") +
    ggplot2::geom_text(data = texts, ggplot2::aes(.data$x, .data$y, label = .data$label,
                                                  hjust = .data$hjust), size = text_size,
                       vjust = 0) +
    ggplot2::annotate("segment", x = 0, xend = total, y = 0.9, yend = 0.9, linewidth = 1.5) +
    ggplot2::annotate("segment", x = edges, xend = edges, y = 0.75, yend = 1.05) +
    ggplot2::geom_text(data = fractions, ggplot2::aes(.data$x, .data$y, label = .data$label),
                       size = text_size) +
    ggplot2::scale_fill_manual(values = fills, name = NULL, drop = FALSE) +
    ggplot2::labs(title = paste(x$component, "component")) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.margin = ggplot2::margin(10, 20, 10, 20))
  attr(p, "fractions") <- shown
  p
}

# centres of the fractions, moved apart so that their labels do not overlap
fraction_positions <- function(edges, width) {
  centres <- (edges[-1] + edges[-length(edges)]) / 2
  for (i in seq_along(centres)[-1]) centres[i] <- max(centres[i], centres[i - 1] + width)
  excess <- centres[length(centres)] + width / 2 - edges[length(edges)]
  if (excess > 0) {
    centres[length(centres)] <- centres[length(centres)] - excess
    for (i in rev(seq_along(centres)[-length(centres)])) {
      centres[i] <- min(centres[i], centres[i + 1] - width)
    }
  }
  centres
}

#' Community matrix ordered by phylogenetic beta diversity
#'
#' The heatmap that Peres-Neto & Kembel (2015) suggest for [phylo_beta_test()]: sites follow
#' an axis of the matrix predicted by the environment and species follow the tree, so that a
#' link between environment and phylogeny appears as clades that replace each other along the
#' sites.
#'
#' @param x Result of [phylo_beta_test()].
#' @param comm The community matrix given to [phylo_beta_test()].
#' @param tree The `phylo` object given to [phylo_beta_test()].
#' @param axis Axis of `x$scores` that orders the sites.
#' @param colours Colour of the largest value, with white for absences.
#' @param text_size Size of the names of species, in mm.
#' @param col,cex,... Deprecated: `col` is now `colours` and `cex` is `text_size`; arguments of
#'   [graphics::image()] are no longer used.
#' @return A ggplot object, with the ordered matrix in the attribute `ordered`.
#' @references
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724.
#' \doi{10.1007/s11258-014-0405-0}
#' @family the whole tree
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
#'                          rmax_range = c(0, 30))
#' env <- runif(50, 0, 100)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' keep <- rowSums(comm) > 0
#' out <- phylo_beta_test(comm[keep, ], tree, env[keep], n_perm = 9)
#' plot_phylo_beta(out, comm[keep, ], tree)
#' @export
plot_phylo_beta <- function(x, comm, tree, axis = 1, colours = "#0072B2", text_size = 2.5,
                            ..., col = NULL, cex = NULL) {
  need_ggplot("plot_phylo_beta()")
  colours <- renamed_arg(colours, col, "colours", "col")
  if (!is.null(cex)) text_size <- 3.5 * renamed_arg(text_size, cex, "text_size", "cex")
  if (!inherits(x, "mc_phylo_beta")) stop("`x` must come from phylo_beta_test()", call. = FALSE)
  check_tree(tree, branch_lengths = FALSE)
  check_colours(colours, 1)
  matched <- present_species(check_phylo_comm(comm, tree, allow_empty = TRUE))
  comm <- matched$comm[rowSums(matched$comm) > 0, , drop = FALSE]
  scores <- as.matrix(x$scores)
  if (nrow(scores) != nrow(comm)) {
    stop("`comm` has ", nrow(comm), " occupied sites and `x` has ", nrow(scores),
         ": use the community given to `phylo_beta_test()`", call. = FALSE)
  }
  if (axis < 1 || axis > ncol(scores)) {
    stop("`axis` must be between 1 and ", ncol(scores), call. = FALSE)
  }
  tree <- ape::ladderize(matched$tree)
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  layout <- tree_layout(tree)
  n_species <- length(tree$tip.label)
  # columns follow the tips in the order of the tree
  ordered <- comm[order(scores[, axis]), tree$tip.label[order(layout$y[seq_len(n_species)])],
                  drop = FALSE]
  n_sites <- nrow(ordered)
  cells <- data.frame(species = rep(seq_len(n_species), each = n_sites),
                      site = rep(seq_len(n_sites), n_species), value = as.vector(ordered))
  # the tree hangs over the columns, root at the top, over a third of the height
  height <- max(3, n_sites / 3)
  depth <- max(layout$x)
  to_y <- function(d) n_sites + 0.5 + height * (1 - d / depth)
  edges <- tree_segments(tree, layout)
  branches <- data.frame(x = edges$y, xend = edges$yend, y = to_y(edges$x), yend = to_y(edges$xend))
  .data <- ggplot2::.data
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = cells, ggplot2::aes(.data$species, .data$site, fill = .data$value)) +
    ggplot2::geom_segment(data = branches, ggplot2::aes(x = .data$x, xend = .data$xend,
                                                        y = .data$y, yend = .data$yend),
                          linewidth = 0.4) +
    ggplot2::scale_fill_gradient(low = "white", high = colours, name = "Abundance") +
    ggplot2::scale_x_continuous(breaks = seq_len(n_species), labels = colnames(ordered),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = paste("Sites, ordered by axis", axis)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                       size = text_size * ggplot2::.pt),
                   panel.grid = ggplot2::element_blank())
  attr(p, "ordered") <- ordered
  p
}
