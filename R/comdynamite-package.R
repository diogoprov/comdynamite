#' @section Analyse data:
#' \itemize{
#'   \item Which clades respond to the environment or to biogeographic events?
#'     [node_contrasts()] codes the nodes, [biogeography_matrix()] the events and
#'     [node_test()] tests them. [plot_node_test()] and [plot_node_heatmap()] show the results.
#'   \item Does phylogenetic composition follow the environment across the whole tree?
#'     [pcc_test()] and [pcps_test()] work on phylogenetic eigenvectors and [phylo_beta_test()] on the
#'     phylogenetic distances between communities, from [phylo_beta()], shown by
#'     [plot_phylo_beta()].
#'   \item Which phylogenetic axis relates to a predictor? [phylo_vector()].
#'   \item How much of trait or phylogenetic variation is due to environment and space?
#'     [trait_partition()] and [plot_trait_partition()].
#'   \item Are sites neighbours? [torus_permutations()] gives restricted null models for a
#'     torus, and [rigid_permutations()] and [patch_permutations()] for the events of a
#'     landscape with edges.
#' }
#'
#' @section Simulate data:
#' \itemize{
#'   \item Dynamics with growth, competition and dispersal: [simulate_metacommunity()], from
#'     [landscape_implicit()] or [landscape_grid()], [species_traits()] and
#'     [competition_matrix()].
#'   \item Static communities along gradients: [simulate_static()] and [simulate_logistic()].
#'   \item Assembly by rules, without dynamics: [assemble_community()].
#'   \item Traits on a tree: [species_traits()], [trait_bm()] and [traits_ou()].
#' }
#'
#' @section Where to start:
#' The articles at <https://provetelab.org/comdynamite/articles/> go through each question.
#' "Preparing your data" shows the format that the analyses expect.
#' @keywords internal
"_PACKAGE"
