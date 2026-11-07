#' Results as data frames, summaries and plots
#'
#' `as.data.frame()` gives the main table of a result. `summary()` gives the same table, with
#' the tests of [node_test()] in order of `p_adjusted`. `plot()` calls [plot_node_test()] or
#' [plot_trait_partition()].
#'
#' @param x,object Result of [node_test()], [pcc_test()], [phylo_beta_test()],
#'   [phylo_vector()] or [trait_partition()].
#' @param row.names,optional Not used.
#' @param tree The `phylo` object that gave the contrasts of [node_contrasts()].
#' @param ... Passed to the plot function, and not used otherwise.
#' @return A data frame, or for `plot()` the value of the plot function.
#' @examples
#' nodes <- node_contrasts(daphniids$tree)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, n_perm = 9)
#' head(as.data.frame(out), 3)
#' head(summary(out), 3)
#' plot(out, daphniids$tree, test = "site")
#' @name result-methods
NULL

#' @rdname result-methods
#' @export
as.data.frame.mc_node_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$nodes
}

#' @rdname result-methods
#' @export
summary.mc_node_test <- function(object, ...) {
  table <- object$nodes[!is.na(object$nodes$p_max), ]
  table <- table[order(table$p_adjusted, table$p_max), ]
  rownames(table) <- NULL
  table
}

#' @rdname result-methods
#' @export
plot.mc_node_test <- function(x, tree, ...) {
  plot_node_test(x, tree, ...)
}

#' @rdname result-methods
#' @export
as.data.frame.mc_pcc <- function(x, row.names = NULL, optional = FALSE, ...) {
  whole_tree_row(x)
}

#' @rdname result-methods
#' @export
summary.mc_pcc <- function(object, ...) {
  whole_tree_row(object)
}

#' @rdname result-methods
#' @export
as.data.frame.mc_phylo_beta <- function(x, row.names = NULL, optional = FALSE, ...) {
  whole_tree_row(x)
}

#' @rdname result-methods
#' @export
summary.mc_phylo_beta <- function(object, ...) {
  whole_tree_row(object)
}

#' @rdname result-methods
#' @export
as.data.frame.mc_phylo_vector <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(vector = names(x$correlation), correlation = unname(x$correlation),
             p_site = x$p_site, p_taxa = x$p_taxa, p_max = x$p_max)
}

#' @rdname result-methods
#' @export
summary.mc_phylo_vector <- function(object, ...) {
  as.data.frame(object)
}

#' @rdname result-methods
#' @export
as.data.frame.mc_trait_partition <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$fractions
}

#' @rdname result-methods
#' @export
summary.mc_trait_partition <- function(object, ...) {
  object$fractions
}

#' @rdname result-methods
#' @export
plot.mc_trait_partition <- function(x, ...) {
  plot_trait_partition(x, ...)
}

whole_tree_row <- function(x) {
  data.frame(r2 = x$r2, statistic = x$statistic, p_site = x$p_site, p_taxa = x$p_taxa,
             p_max = x$p_max)
}
