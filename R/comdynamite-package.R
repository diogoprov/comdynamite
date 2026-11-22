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
#'   \item Traits on a tree: [species_traits()], [trait_bm()] and [trait_ou()], with
#'     correlated traits from [rate_matrix()].
#' }
#'
#' @section Arguments:
#' Arguments that choose among options, such as `method`, `null` or `model`, are matched by
#' [match.arg()] and are case-sensitive: `null = "variogram"` works and `null = "Variogram"`
#' stops with the list of valid options. Each scalar argument is checked for its length, type
#' and range before any computation, and missing or infinite values in `comm`, `env` and the
#' other data inputs stop the function with a message that names the argument; remove or
#' impute missing values before the analysis. Factors used as labels (the site and species
#' columns of the data on individuals, `group` of [drift_test()], `regimes` of the DAMOCLES
#' functions, `strata` of [env_surrogates()]) are taken as unordered: the order of the levels
#' of an ordered factor is ignored. Factors and character columns are refused as predictors in
#' `env`; code them as dummy variables.
#'
#' Inputs are converted once, at the start: community tables become double matrices (integer
#' counts included, so the `comm` of [mc_data()] is double), and site and species labels in
#' results are character, even when the input columns were factors. Attributes and classes of
#' numeric columns, such as units of the units package, are dropped without a message: values
#' are taken in the units given, and results are in those units. The analyses take either an
#' [mc_data()] object, which keeps the parts matched by site and species names, or the parts
#' as separate arguments (`comm`, `tree`, `traits`, `env`), which are matched by the same names
#' inside each function.
#'
#' @section Citation:
#' Cite the package with `citation("comdynamite")` and, with it, the original papers of the
#' methods you use. Each function that implements a published method shows these papers the
#' first time it runs in a session, and [method_references()] lists them by function; see
#' [method_references()] for how to turn the message off.
#'
#' @section Audience:
#' The package is written for community ecologists, biogeographers and evolutionary biologists
#' who analyse species composition, traits and phylogenies across sites, and for
#' methodologists who want to check, by simulation, whether a method recovers a known process.
#' It assumes familiarity with R and with the methods it implements, whose original papers each
#' help page cites.
#'
#' @section Where to start:
#' The articles at <https://provetelab.org/comdynamite/articles/> go through each question.
#' "Preparing your data" shows the format that the analyses expect.
#' @keywords internal
"_PACKAGE"

# computed variable of ggplot2::geom_histogram(), used through ggplot2::after_stat()
utils::globalVariables("density")
