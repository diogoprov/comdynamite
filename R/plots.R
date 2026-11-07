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
#' @param col Colours of predictors and biogeographic events.
#' @param cex Size of the text.
#' @param legend Add a legend?
#' @param ... Other arguments of [ape::plot.phylo()].
#' @return Invisibly, a data frame with one row per symbol, with `node`, `set`, `variable`,
#'   `p`, `p_adjusted` and `corrected`, which is `TRUE` for solid symbols.
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing. *Journal of the Royal Statistical
#' Society: Series B*, 57, 289--300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @examples
#' nodes <- node_contrasts(daphniids$tree)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
#'                  n_perm = 49)
#' plot_node_test(out, daphniids$tree, test = "site")
#' \donttest{
#' # solid symbols need more permutations
#' out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
#'                  n_perm = 999)
#' plot_node_test(out, daphniids$tree, test = "site")
#' }
#' @export
plot_node_test <- function(x, tree, test = c("max", "site", "joint"),
                           type = c("phylogram", "fan"), level = 0.05,
                           col = c("#0072B2", "#D55E00"), cex = 0.8, legend = TRUE, ..., alpha = NULL) {
  level <- deprecated_level(level, alpha)
  stopifnot(inherits(x, "mc_node_test"), inherits(tree, "phylo"), length(col) == 2)
  test <- match.arg(test)
  type <- match.arg(type)
  check_node_tree(x, tree)
  marks <- node_marks(x, test, level)
  # names of species and of variables can go beyond the plot region
  old <- graphics::par(xpd = NA)
  on.exit(graphics::par(old))
  ape::plot.phylo(tree, type = type, cex = cex, ...)
  last <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  for (node in unique(marks$node)) {
    here <- marks[marks$node == node, ]
    # one symbol per set and level, with solid ones next to the node
    groups <- split(here, list(here$set == "biogeography", !here$corrected), drop = TRUE)
    side <- vapply(groups, function(g) if (g$set[1] == "biogeography") -1 else 1, numeric(1))
    at <- as.integer(node)
    angle <- if (type == "fan") atan2(last$yy[at], last$xx[at]) else 0
    steps <- symbol_steps(side, graphics::strheight("M", cex = cex))
    x <- last$xx[at] - steps * sin(angle)
    y <- last$yy[at] + steps * cos(angle)
    room <- if (type == "fan") sqrt(x^2 + y^2) else x - graphics::par("usr")[1]
    for (i in seq_along(groups)) draw_symbol(x[i], y[i], groups[[i]], col, cex, angle, room[i])
  }
  if (legend) {
    graphics::legend("bottomleft", c("Predictor", "Biogeography", "After correction", "Uncorrected"),
                     pch = c(19, 17, 19, 21), col = c(col, "black", "black"), bty = "n", cex = cex)
  }
  invisible(marks)
}

#' Plot the semipartial correlations of a node test
#'
#' Draws a table of nodes by variables with the semipartial correlations of [node_test()],
#' as in Table S1 of Leibold et al. (2010). Bold values with an asterisk remain significant
#' after the false discovery rate, as in [plot_node_test()]; other bold values are
#' significant only as isolated tests. Blank cells could not be calculated.
#'
#' Rows are the node numbers of the tree. To find them, call `ape::nodelabels()` after
#' [plot_node_test()].
#'
#' @inheritParams plot_node_test
#' @param col Colours of negative and positive values.
#' @param ... Other arguments of [graphics::image()].
#' @return Invisibly, the matrix of semipartial correlations, with predictors followed by
#'   biogeographic events.
#' @references
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @examples
#' nodes <- node_contrasts(daphniids$tree)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
#'                  n_perm = 49)
#' plot_node_heatmap(out, test = "site")
#' @export
plot_node_heatmap <- function(x, test = c("max", "site", "joint"), level = 0.05,
                              col = c("#0072B2", "#D55E00"), cex = 0.8, ..., alpha = NULL) {
  level <- deprecated_level(level, alpha)
  stopifnot(inherits(x, "mc_node_test"), length(col) == 2)
  test <- match.arg(test)
  values <- cbind(x$semipartial, x$biogeography$semipartial)
  if (all(is.na(values))) stop("no node was tested, so there is nothing to show", call. = FALSE)
  marks <- node_marks(x, test, level)
  n_nodes <- nrow(values)
  n_variables <- ncol(values)
  rows <- rev(seq_len(n_nodes))
  limit <- max(abs(values), na.rm = TRUE)
  ramp <- grDevices::colorRampPalette(c(col[1], "white", col[2]))(101)
  old <- graphics::par(mar = c(6, 4, 1, 1))
  on.exit(graphics::par(old))
  graphics::image(seq_len(n_variables), seq_len(n_nodes), t(values[rows, , drop = FALSE]),
                  zlim = c(-limit, limit), col = ramp, axes = FALSE, xlab = "", ylab = "Node", ...)
  graphics::axis(1, seq_len(n_variables), colnames(values), las = 2, tick = FALSE, cex.axis = cex)
  graphics::axis(2, seq_len(n_nodes), rownames(values)[rows], las = 1, tick = FALSE, cex.axis = cex)
  graphics::box()
  significant <- matrix(0, n_nodes, n_variables)
  significant[cbind(match(marks$node, rownames(values)), marks$column)] <- 1 + marks$corrected
  cells <- which(!is.na(values), arr.ind = TRUE)
  strength <- significant[cells]
  graphics::text(cells[, 2], n_nodes + 1 - cells[, 1],
                 paste0(format(round(values[cells], 2), nsmall = 2), c("", "", "*")[strength + 1]),
                 font = ifelse(strength > 0, 2, 1), cex = cex * 0.85)
  invisible(values)
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

# symbols in a column at the node, predictors above and events below
symbol_steps <- function(side, height) {
  side * (stats::ave(side, side, FUN = seq_along) - 0.4) * height * 1.3
}

# labels go towards the root, over a white box, along the radius of a fan
draw_symbol <- function(x, y, group, col, cex, angle, room) {
  colour <- col[1 + (group$set[1] == "biogeography")]
  label <- paste(group$variable, collapse = ", ")
  size <- cex * 0.85
  width <- graphics::strwidth(label, cex = size)
  height <- graphics::strheight(label, cex = size)
  gap <- graphics::strwidth("n", cex = size)
  place <- label_place(x, y, width, gap, angle, room)
  along <- c(cos(angle), sin(angle)) * width / 2
  across <- c(-sin(angle), cos(angle)) * height / 2
  graphics::polygon(place$x + c(-1, 1, 1, -1) * along[1] + c(-1, -1, 1, 1) * across[1],
                    place$y + c(-1, 1, 1, -1) * along[2] + c(-1, -1, 1, 1) * across[2],
                    col = "white", border = NA)
  graphics::text(place$x, place$y, label, srt = place$srt, col = colour, cex = size)
  shape <- if (group$set[1] == "biogeography") c(17, 24) else c(19, 21)
  graphics::points(x, y, pch = shape[2 - group$corrected[1]], col = colour, bg = "white")
}

# centre and rotation of a label; it goes away from the root only if there is no room,
# and turns on the left of a fan so that it is not upside down
label_place <- function(x, y, width, gap, angle, room) {
  direction <- if (room > gap + width) -1 else 1
  shift <- direction * (gap + width / 2)
  list(x = x + shift * cos(angle), y = y + shift * sin(angle),
       srt = angle * 180 / pi + if (cos(angle) < 0) 180 else 0)
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
#' @param col Colour of significant segments.
#' @param labels Names of the two sets of predictors.
#' @param cex Size of the text.
#' @param ... Other arguments of [graphics::plot.window()].
#' @return Invisibly, the rows of `x$fractions` for `a`, `b`, `c` and `d`, with `shown`, the
#'   values drawn.
#' @references
#' Desdevises, Y., Legendre, P., Azouzi, L., & Morand, S. (2003). Quantifying
#' phylogenetically structured environmental variation. *Evolution*, 57, 2647--2652.
#' \doi{10.1111/j.0014-3820.2003.tb01508.x}
#'
#' Legendre, P. (2008). Studying beta diversity: ecological variation partitioning by
#' multiple regression and canonical analysis. *Journal of Plant Ecology*, 1, 3--8.
#' \doi{10.1093/jpe/rtm001}
#' @examples
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
                                 col = "#0072B2", labels = c("Environment", "Space"),
                                 cex = 0.9, ..., alpha = NULL) {
  level <- deprecated_level(level, alpha)
  stopifnot(inherits(x, "mc_trait_partition"), length(labels) == 2)
  value <- match.arg(value)
  if (identical(test, "taxa")) {
    warning("`test = \"taxa\"` is deprecated, use `\"species\"`", call. = FALSE)
    test <- "species"
  }
  test <- match.arg(test)
  table <- x$fractions
  shown <- table[match(c("a", "b", "c", "d"), table$fraction), ]
  shown$shown <- pmax(shown[[value]], 0)
  p <- stats::setNames(shown[[paste0("p_", test)]], shown$fraction)
  fill <- function(fraction) {
    if (is.na(p[[fraction]])) return("grey70")
    if (p[[fraction]] <= level) col else grDevices::adjustcolor(col, 0.3)
  }
  edges <- c(0, cumsum(shown$shown))
  old <- graphics::par(mar = c(3, 1, 3, 1), xpd = NA)
  on.exit(graphics::par(old))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, edges[5]), ylim = c(0, 3), ...)
  segment <- function(i, y, colour) {
    graphics::rect(edges[i], y, edges[i + 1], y + 0.7, col = colour, border = "white")
  }
  segment(1, 2, fill("a"))
  segment(2, 2, fill("b"))
  segment(2, 1.1, fill("b"))
  segment(3, 1.1, fill("c"))
  graphics::text(edges[1], 2.85, labels[1], adj = c(0, 0), cex = cex)
  # beside its bar, on the side with more room
  on_right <- edges[5] - edges[4] >= edges[2]
  gap <- graphics::strwidth("n", cex = cex)
  graphics::text(if (on_right) edges[4] + gap else edges[2] - gap, 1.45, labels[2],
                 adj = c(1 - on_right, 0.5), cex = cex)
  graphics::text(mean(edges[4:5]), 2.35, "Unexplained", cex = cex)
  graphics::segments(edges[1], 0.9, edges[5], 0.9, lwd = 4)
  graphics::segments(edges, 0.75, edges, 1.05)
  graphics::text(fraction_positions(edges, graphics::strwidth("[a] 0.00", cex = cex)), 0.45,
                 paste0("[", shown$fraction, "] ", format(round(shown$shown, 2), nsmall = 2)),
                 cex = cex)
  graphics::title(main = paste(x$component, "component"))
  graphics::legend("bottom", c("Significant", "Not significant", "Not tested"), bty = "n",
                   fill = c(col, grDevices::adjustcolor(col, 0.3), "grey70"), border = NA,
                   horiz = TRUE, cex = cex * 0.9, inset = c(0, -0.08))
  invisible(shown)
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
#' @param col Colour of the largest value, with white for absences.
#' @param cex Size of the text.
#' @param ... Passed to [graphics::image()].
#' @return The ordered matrix, invisibly.
#' @references
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724.
#' \doi{10.1007/s11258-014-0405-0}
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
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
plot_phylo_beta <- function(x, comm, tree, axis = 1, col = "#0072B2", cex = 0.7, ...) {
  stopifnot(inherits(x, "mc_phylo_beta"), inherits(tree, "phylo"))
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
  tips <- tree$edge[tree$edge[, 2] <= length(tree$tip.label), 2]
  ordered <- comm[order(scores[, axis]), tree$tip.label[tips], drop = FALSE]
  n_species <- ncol(ordered)
  n_sites <- nrow(ordered)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  graphics::layout(matrix(1:2), heights = c(1, 3))
  # without the extra 4% of the axis, tips line up with the columns of the image
  graphics::par(mar = c(0, 4, 1, 1), xaxs = "i")
  ape::plot.phylo(tree, direction = "downwards", show.tip.label = FALSE,
                  x.lim = c(0.5, n_species + 0.5), no.margin = FALSE)
  graphics::par(mar = c(7, 4, 0, 1))
  ramp <- grDevices::colorRampPalette(c("white", col))(max(2, min(max(ordered) + 1, 50)))
  graphics::image(seq_len(n_species), seq_len(n_sites), t(ordered), col = ramp, axes = FALSE,
                  xlab = "", ylab = paste("Sites, ordered by axis", axis), ...)
  graphics::axis(1, seq_len(n_species), colnames(ordered), las = 2, tick = FALSE, cex.axis = cex)
  graphics::box()
  invisible(ordered)
}
