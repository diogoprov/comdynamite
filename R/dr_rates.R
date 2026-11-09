#' Tip diversification rates (DR) of Jetz et al. (2012)
#'
#' `dr_rates()` computes the DR statistic of Jetz et al. (2012), a measure of the recent rate
#' of speciation of each tip from the tree alone: the inverse of its equal-splits measure of
#' evolutionary distinctiveness, `DR_i = 1 / sum_j l_j / 2^(j - 1)`, where `l_j` is the length of
#' the `j`-th edge on the path from tip `i` to the root, the terminal edge first. Tips in
#' recently and rapidly diversifying clades have short, heavily shared paths and a high DR.
#'
#' With `per = "edge"`, each edge gets the geometric mean of the DR of the tips below it (the
#' DR of its tip for a terminal edge), a value per edge in the order of `tree$edge` that can be
#' given as the `edge_covariate` of [fit_damocles()] or [damocles_samples()]. DR is defined for
#' tips; the value of an internal edge summarises its descendants and is not an estimate of the
#' rate along that edge, unlike the branch-specific rates of model-based methods such as ClaDS.
#' `data-raw/damocles_clads.R` checks the tip values against the equal-splits measure of
#' `picante::evol.distinct()` and compares DR with ClaDS as covariates.
#'
#' @param tree A rooted `phylo` tree with branch lengths.
#' @param per `"tip"` for one value per tip, named by the tips, or `"edge"` for one value per
#'   edge, in the order of `tree$edge`.
#' @return A numeric vector of DR values.
#' @references Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The
#'   global diversity of birds in space and time. *Nature*, 491, 444--448.
#'   \doi{10.1038/nature11631}
#' @seealso [fit_damocles()], [damocles_samples()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(10)
#' tr <- treeshapesim::sim_tree_shape(20)
#' head(dr_rates(tr))
#' head(dr_rates(tr, per = "edge"))
#' @export
dr_rates <- function(tree, per = c("tip", "edge")) {
  per <- match.arg(per)
  if (!inherits(tree, "phylo")) stop("`tree` must be a `phylo` object", call. = FALSE)
  if (is.null(tree$edge.length) || anyNA(tree$edge.length) || any(tree$edge.length < 0)) {
    stop("`tree` needs branch lengths, zero or more", call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  parent_edge <- integer(n_tip + tree$Nnode)
  parent_edge[tree$edge[, 2]] <- seq_len(nrow(tree$edge))
  dr <- vapply(seq_len(n_tip), function(tip) {
    total <- 0
    weight <- 1
    node <- tip
    while (parent_edge[node] > 0) {
      e <- parent_edge[node]
      total <- total + tree$edge.length[e] * weight
      weight <- weight / 2
      node <- tree$edge[e, 1]
    }
    1 / total
  }, numeric(1))
  names(dr) <- tree$tip.label
  if (per == "tip") return(dr)
  vapply(seq_len(nrow(tree$edge)), function(e) {
    below <- damocles_below(tree, tree$edge[e, 2])[seq_len(n_tip)]
    exp(mean(log(dr[below])))
  }, numeric(1))
}
