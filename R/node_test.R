#' Node contrasts of a phylogeny
#'
#' Codes each node of a tree as a contrast between its two descendant clades, as
#' in Leibold et al. (2010). One clade starts at -1 and the other at 1, and the
#' value is halved at each node on the way to the tips. Species outside the node
#' get zero. Nodes that join two tips are left out, as in the original code,
#' because they contrast two single species.
#'
#' @param tree A `phylo` object with dichotomous nodes.
#' @return Matrix with one row per tip and one column per node, named by the node
#'   numbers of `tree`.
#' @references
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @examples
#' node_contrasts(ape::read.tree(text = "((a,b),(c,(d,e)));"))
#' @export
node_contrasts <- function(tree) {
  stopifnot(inherits(tree, "phylo"))
  n_tips <- length(tree$tip.label)
  children <- split(tree$edge[, 2], tree$edge[, 1])
  if (any(lengths(children) != 2)) stop("the tree must be dichotomous and rooted", call. = FALSE)
  # value of each tip below a node, starting from the value of the node
  spread <- function(node, value) {
    if (node <= n_tips) return(stats::setNames(value, node))
    unlist(lapply(children[[as.character(node)]], spread, value = value / 2),
           use.names = TRUE)
  }
  informative <- Filter(function(node) any(children[[node]] > n_tips), names(children))
  contrasts <- vapply(informative, function(node) {
    pair <- children[[node]]
    out <- numeric(n_tips)
    for (side in 1:2) {
      values <- spread(pair[side], c(-1, 1)[side])
      out[as.integer(sub(".*\\.", "", names(values)))] <- values
    }
    out
  }, numeric(n_tips))
  matrix(contrasts, n_tips, dimnames = list(tree$tip.label, informative))
}

#' Test clades against the environment, node by node
#'
#' Implements the regression of Leibold et al. (2010). For each node, sites without
#' species of that node are removed, the contrast is standardised with species
#' weights and summed within sites, and the result is regressed on the predictors
#' standardised with site weights, using site richness as weights. Each node shows
#' whether the two clades that split there differ in their response to the predictors.
#'
#' The test statistic is the explained over the residual sum of squares. It is
#' tested by permuting sites, and by replacing the species of the node with a random
#' set of the same size from all species. Following Peres-Neto et al. (2012) and
#' ter Braak et al. (2012), the p-value to report is the larger of the two. With many
#' nodes, these p-values need a correction for multiple tests.
#'
#' R2, adjusted R2 and the statistic equal those of the MATLAB function
#' `MetacommunityPhylogenetics` of the authors, checked in Octave. Three parts differ:
#' \itemize{
#'   \item In the original, the species permutations shuffle only the first columns
#'     of the community matrix, as many as there are species in the node, so they
#'     do not draw from all species. Here they do.
#'   \item The original keeps the total of the observed data after permuting
#'     species. Here it is recalculated.
#'   \item The semipartial values of the original use the raw predictors, without an
#'     intercept, whereas the full model uses standardised predictors. Here
#'     `partial_r2` is the loss of R2 when a predictor is removed, with standardised
#'     predictors in both models.
#' }
#'
#' As in the original, R2 sums squares without site weights, although the fit
#' uses them, so it is not the R2 of a weighted regression, and the adjusted R2
#' uses the degrees of freedom of an unweighted model.
#'
#' Nodes for which a predictor has no variance among the remaining sites, or with too
#' few sites, get `NA`. The original stops with an error.
#'
#' @param comm Matrix with one row per site and one column per species, named
#'   as the rows of `nodes`.
#' @param nodes Matrix of contrasts with one row per species, from
#'   [node_contrasts()].
#' @param env Predictor, or matrix or data frame of predictors, with one row per
#'   site.
#' @param n_perm Number of permutations of sites and of species.
#' @return List of class `mc_node_test` with `nodes`, a data frame with one row per
#'   node with `node`, `n_species`, `n_sites`, `r2`, `r2_adjusted`, `statistic`,
#'   `p_site`, `p_taxa` and `p_max`; `partial_r2`, a matrix of nodes by predictors; and
#'   `p_partial`, with the larger of the two p-values for each of them.
#' @references
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#'
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of
#' spatial contingency and environmental filtering on metacommunity phylogenetics.
#' *Ecology*, 93, S14--S30. \doi{10.1890/11-0494.1}
#'
#' ter Braak, C. J. F., Cormont, A., & Dray, S. (2012). Improved testing of species
#' traits-environment relationships in the fourth-corner problem. *Ecology*, 93,
#' 1525--1526.
#' @seealso [pcc_test()] and [phylo_beta_test()] test the whole tree at once.
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(12)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
#'                          rmax_range = c(0, 30))
#' env <- runif(40, 0, 100)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' keep <- rowSums(comm) > 0
#' node_test(comm[keep, ], node_contrasts(tree), env[keep], n_perm = 19)
#' @export
node_test <- function(comm, nodes, env, n_perm = 999) {
  comm <- as.matrix(comm) * 1
  nodes <- as.matrix(nodes)
  env <- as.matrix(env)
  stopifnot(setequal(colnames(comm), rownames(nodes)), nrow(env) == nrow(comm),
            n_perm >= 1, all(comm >= 0))
  comm <- comm[, rownames(nodes), drop = FALSE]
  names_env <- colnames(env) %||% paste0("predictor", seq_len(ncol(env)))
  results <- lapply(seq_len(ncol(nodes)), function(i) one_node(comm, nodes[, i], env, n_perm))
  matrix_of <- function(field) {
    out <- do.call(rbind, lapply(results, `[[`, field))
    dimnames(out) <- list(colnames(nodes), names_env)
    out
  }
  table <- do.call(rbind, lapply(results, `[[`, "summary"))
  structure(list(nodes = data.frame(node = colnames(nodes) %||% seq_len(ncol(nodes)), table),
                 partial_r2 = matrix_of("partial_r2"), p_partial = matrix_of("p_partial")),
            class = "mc_node_test")
}

#' @export
print.mc_node_test <- function(x, ...) {
  tested <- x$nodes[!is.na(x$nodes$p_max), ]
  cat("<mc_node_test>", nrow(x$nodes), "nodes,", nrow(tested), "tested\n")
  cat("nodes with p (max) <= 0.05, before correction:", sum(tested$p_max <= 0.05), "\n")
  invisible(x)
}

one_node <- function(comm, contrast, env, n_perm) {
  members <- which(contrast != 0)
  values <- contrast[members]
  observed <- node_fit(comm[, members, drop = FALSE], values, env)
  count <- function(null, value) {
    (sum(null >= value, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
  }
  blank <- rep(NA_real_, ncol(env))
  if (is.na(observed$statistic)) {
    return(list(summary = data.frame(n_species = length(members), n_sites = observed$n_sites,
                                     r2 = NA_real_, r2_adjusted = NA_real_, statistic = NA_real_,
                                     p_site = NA_real_, p_taxa = NA_real_, p_max = NA_real_),
                partial_r2 = blank, p_partial = blank))
  }
  null_of <- function(draw) {
    t(replicate(n_perm, {
      fit <- draw()
      c(fit$statistic, fit$partial_f)
    }))
  }
  by_site <- null_of(function() {
    node_fit(comm[sample.int(nrow(comm)), members, drop = FALSE], values, env)
  })
  by_taxa <- null_of(function() {
    node_fit(comm[, sample.int(ncol(comm), length(members)), drop = FALSE], values, env)
  })
  p_site <- count(by_site[, 1], observed$statistic)
  p_taxa <- count(by_taxa[, 1], observed$statistic)
  p_partial <- vapply(seq_len(ncol(env)), function(j) {
    max(count(by_site[, j + 1], observed$partial_f[j]),
        count(by_taxa[, j + 1], observed$partial_f[j]))
  }, numeric(1))
  list(summary = data.frame(n_species = length(members), n_sites = observed$n_sites,
                            r2 = observed$r2, r2_adjusted = observed$r2_adjusted,
                            statistic = observed$statistic, p_site = p_site, p_taxa = p_taxa,
                            p_max = max(p_site, p_taxa)),
       partial_r2 = observed$partial_r2, p_partial = p_partial)
}

# statistics of one node; all NA when the model cannot be fitted
node_fit <- function(incidence, values, env) {
  occupied <- rowSums(incidence) > 0
  incidence <- incidence[occupied, , drop = FALSE]
  env <- env[occupied, , drop = FALSE]
  n <- nrow(incidence)
  p <- ncol(env)
  failed <- list(statistic = NA_real_, partial_f = rep(NA_real_, p), n_sites = n)
  site_weights <- rowSums(incidence)
  species_weights <- colSums(incidence)
  if (n < p + 2 || any(apply(env, 2, stats::var) == 0)) return(failed)
  standard <- weighted_standardise(matrix(values), species_weights)
  if (any(!is.finite(standard))) return(failed)
  response <- drop(incidence %*% standard)
  predictors <- weighted_standardise(env, site_weights)
  fit <- weighted_fit(predictors, response, site_weights)
  if (is.null(fit)) return(failed)
  total <- sum(response^2)
  explained <- sum(fit^2)
  residual <- sum((response - fit)^2)
  r2 <- explained / total
  partial_r2 <- vapply(seq_len(p), function(j) {
    if (p == 1) return(r2)
    r2 - sum(weighted_fit(predictors[, -j, drop = FALSE], response, site_weights)^2) / total
  }, numeric(1))
  list(statistic = explained / residual, partial_f = partial_r2 / residual, r2 = r2,
       r2_adjusted = 1 - (n - 1) / (n - p - 1) * (1 - r2), partial_r2 = partial_r2, n_sites = n)
}

weighted_fit <- function(predictors, response, weights) {
  cross <- crossprod(predictors, weights * predictors)
  if (rcond(cross) < .Machine$double.eps) return(NULL)
  drop(predictors %*% solve(cross, crossprod(predictors, weights * response)))
}
