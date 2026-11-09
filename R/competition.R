#' Competition matrix
#'
#' Intraspecific coefficients are 1 and interspecific coefficients lie in
#' \[0, `inter_max`), all multiplied by `scale`. With the default `scale` of
#' 0.05, as in Thompson et al. (2020), a species alone at its optimum with
#' `r_max = 5` reaches an equilibrium of \eqn{(5 - 1) / 0.05 = 80} individuals,
#' and `inter_max` is the strength of interspecific relative to intraspecific
#' competition. The matrix is symmetric.
#'
#' * `"none"`: no interspecific competition.
#' * `"random"`: interspecific coefficients drawn from U(0, `inter_max`).
#' * `"phylo"`: the same uniform values, assigned to pairs of species by rank
#'   of phylogenetic similarity (1 - d / max(d), with d the patristic
#'   distance), so that close relatives compete more strongly while the
#'   distribution of coefficients stays the same as with `"random"`.
#'
#' @param species A `phylo` tree, or the number of species (not for
#'   `type = "phylo"`).
#' @param type Type of interspecific competition.
#' @param inter_max Maximum interspecific coefficient, relative to
#'   intraspecific competition.
#' @param scale Multiplier of all coefficients.
#' @return Square matrix: element \[j, k\] is the effect of species k on
#'   species j.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(5)
#' round(competition_matrix(tr, "phylo", inter_max = 0.5), 3)
#' @export
competition_matrix <- function(species, type = c("none", "random", "phylo"),
                               inter_max = 0.1, scale = 0.05) {
  type <- match.arg(type)
  stopifnot(inter_max >= 0, scale > 0)
  is_tree <- inherits(species, "phylo")
  if (type == "phylo" && !is_tree) stop("type = \"phylo\" needs a 'phylo' tree as `species`", call. = FALSE)
  sp <- if (is_tree) species$tip.label else paste0("sp", seq_len(species))
  k <- length(sp)
  a <- matrix(0, k, k, dimnames = list(sp, sp))
  up <- upper.tri(a)
  vals <- if (type == "none") numeric(sum(up)) else stats::runif(sum(up), 0, inter_max)
  if (type == "phylo") {
    d <- ape::cophenetic.phylo(species)[sp, sp]
    s <- (1 - d / max(d))[up]
    vals <- sort(vals)[rank(s, ties.method = "random")]
  }
  a[up] <- vals
  a <- a + t(a)
  diag(a) <- 1
  scale * a
}
