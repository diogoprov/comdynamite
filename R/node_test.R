#' Node contrasts of a phylogeny
#'
#' Codes each node of a tree as a contrast between its two descendant clades, as
#' in Leibold et al. (2010). One clade starts at -1 and the other at 1, and the
#' value is halved at each node on the way to the tips. Species outside the node
#' get zero.
#'
#' @param tree A `phylo` object with dichotomous nodes.
#' @param cherries Keep nodes that join two tips? The article tests them (its Figure 3),
#'   but the function `CreateNodeContrasts` of the authors leaves them out.
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
node_contrasts <- function(tree, cherries = TRUE) {
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
  informative <- Filter(function(node) cherries || any(children[[node]] > n_tips), names(children))
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
#' set of the same size from all species. The article permutes both at once; following
#' Peres-Neto et al. (2012) and ter Braak et al. (2012), the p-value here is the larger
#' of the two tests. With many nodes these p-values need a correction, and
#' `p_adjusted` gives the false discovery rate that the article uses. This correction
#' needs many permutations, as the smallest p-value is 1 / (`n_perm` + 1), so the default
#' is the 9999 of the article.
#'
#' With few species at a node, permutations of species have little power, so `p_max`
#' is conservative: with the [daphniids] of the article, of 13 species, no predictor
#' passes it, whereas sites alone recover most of the environmental results of its
#' Figure 3 (see the article on validation).
#'
#' Each predictor gets a semipartial correlation, the square root of the loss of R2
#' when it is removed, with the sign of its slope when it is alone in the model.
#'
#' Biogeographic events follow the article: a column of `biogeography` has 1 for sites
#' in one area, -1 for sites in the other and 0 for sites in neither. Each event is
#' tested only in its sites, controlling for `env`. Other contrasts are allowed, for
#' example a column without zeros, which compares one area with all the others. Events
#' with few sites in one area give unstable values.
#'
#' R2, adjusted R2 and the statistic equal those of the MATLAB function
#' `MetacommunityPhylogenetics` of the authors, checked in Octave. Three parts of that
#' function differ from this one:
#' \itemize{
#'   \item In the original, the species permutations shuffle only the first columns
#'     of the community matrix, as many as there are species in the node, so they
#'     do not draw from all species. Here they do.
#'   \item The original keeps the total of the observed data after permuting
#'     species. Here it is recalculated.
#'   \item The semipartial values of the original use the raw predictors, without an
#'     intercept, and multiply by the absolute slope. The article uses standardised
#'     predictors and divides by it, to keep only the sign, which is what we do.
#' }
#'
#' As in the original, R2 sums squares without site weights, although the fit
#' uses them, so it is not the R2 of a weighted regression, and the adjusted R2
#' uses the degrees of freedom of an unweighted model.
#'
#' Nodes for which a predictor has no variance among the remaining sites, or with too
#' few sites, get `NA`. The original stops with an error.
#'
#' @param comm Matrix with one row per site and one column per species. Species are
#'   matched to the rows of `nodes` by name, in any order, and those in only one of the
#'   two are dropped with a warning.
#' @param nodes Matrix of contrasts with one row per species, from
#'   [node_contrasts()].
#' @param env Predictor, or matrix or data frame of predictors, with one row per
#'   site.
#' @param biogeography Optional matrix of biogeographic events with one row per site
#'   and values of 1, -1 and 0.
#' @param n_perm Number of permutations of sites and of species.
#' @return List of class `mc_node_test` with `nodes`, a data frame with one row per
#'   node with `node`, `n_species`, `n_sites`, `r2`, `r2_adjusted`, `statistic`,
#'   `p_site`, `p_taxa`, `p_max` and `p_adjusted`; `semipartial`, a matrix of nodes by
#'   predictors; `p_partial`, with the larger of the two p-values for each of them;
#'   `p_partial_site`, with the p-values of sites alone; and `biogeography`, a list with
#'   the same three matrices for the events, or `NULL`.
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
#' # daphniids of Leibold et al. (2010)
#' nodes <- node_contrasts(daphniids$tree)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, n_perm = 19)
#' out
#' round(out$semipartial, 2)
#' @export
node_test <- function(comm, nodes, env, biogeography = NULL, n_perm = 9999) {
  comm <- as.matrix(comm) * 1
  nodes <- as.matrix(nodes)
  env <- as.matrix(env)
  stopifnot(nrow(env) == nrow(comm), n_perm >= 1, all(comm >= 0))
  species <- match_species(comm, rownames(nodes), "`nodes`")
  comm <- comm[, species, drop = FALSE]
  nodes <- nodes[species, , drop = FALSE]
  node_names <- colnames(nodes) %||% seq_len(ncol(nodes))
  env_names <- colnames(env) %||% paste0("predictor", seq_len(ncol(env)))
  per_node <- function(run) lapply(seq_len(ncol(nodes)), function(i) run(nodes[, i]))
  matrix_of <- function(results, field, names) {
    matrix(unlist(lapply(results, `[[`, field)), length(results), byrow = TRUE,
           dimnames = list(node_names, names))
  }
  results <- per_node(function(contrast) one_node(comm, contrast, env, n_perm))
  table <- data.frame(node = node_names, do.call(rbind, lapply(results, `[[`, "summary")))
  table$p_adjusted <- stats::p.adjust(table$p_max, "BH")
  structure(list(nodes = table, semipartial = matrix_of(results, "semipartial", env_names),
                 p_partial = matrix_of(results, "p_partial", env_names),
                 p_partial_site = matrix_of(results, "p_partial_site", env_names),
                 biogeography = if (!is.null(biogeography)) {
                   biogeography <- as.matrix(biogeography)
                   stopifnot(nrow(biogeography) == nrow(comm), all(biogeography %in% c(-1, 0, 1)))
                   events <- colnames(biogeography) %||% paste0("event", seq_len(ncol(biogeography)))
                   by_event <- lapply(seq_len(ncol(biogeography)), function(e) {
                     sites <- biogeography[, e] != 0
                     predictors <- cbind(biogeography[sites, e], env[sites, , drop = FALSE])
                     per_node(function(contrast) {
                       one_node(comm[sites, , drop = FALSE], contrast, predictors, n_perm)
                     })
                   })
                   first <- function(field) {
                     vapply(by_event, function(event) {
                       vapply(event, function(node) node[[field]][1], numeric(1))
                     }, numeric(ncol(nodes)))
                   }
                   shaped <- function(x) {
                     matrix(x, ncol(nodes), length(events), dimnames = list(node_names, events))
                   }
                   list(semipartial = shaped(first("semipartial")),
                        p_partial = shaped(first("p_partial")),
                        p_partial_site = shaped(first("p_partial_site")))
                 }),
            class = "mc_node_test")
}

#' @export
print.mc_node_test <- function(x, ...) {
  tested <- x$nodes[!is.na(x$nodes$p_max), ]
  cat("<mc_node_test>", nrow(x$nodes), "nodes,", nrow(tested), "tested\n")
  cat("nodes with p (max) <= 0.05:", sum(tested$p_max <= 0.05), "before correction,",
      sum(tested$p_adjusted <= 0.05), "after\n")
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
                semipartial = blank, p_partial = blank, p_partial_site = blank))
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
  partial_site <- vapply(seq_len(ncol(env)), function(j) {
    count(by_site[, j + 1], observed$partial_f[j])
  }, numeric(1))
  partial_taxa <- vapply(seq_len(ncol(env)), function(j) {
    count(by_taxa[, j + 1], observed$partial_f[j])
  }, numeric(1))
  list(summary = data.frame(n_species = length(members), n_sites = observed$n_sites,
                            r2 = observed$r2, r2_adjusted = observed$r2_adjusted,
                            statistic = observed$statistic, p_site = p_site, p_taxa = p_taxa,
                            p_max = max(p_site, p_taxa)),
       semipartial = observed$semipartial, p_partial = pmax(partial_site, partial_taxa),
       p_partial_site = partial_site)
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
  slopes <- colSums(site_weights * predictors * response)
  list(statistic = explained / residual, partial_f = partial_r2 / residual, r2 = r2,
       r2_adjusted = 1 - (n - 1) / (n - p - 1) * (1 - r2),
       semipartial = sign(slopes) * sqrt(pmax(partial_r2, 0)), n_sites = n)
}

weighted_fit <- function(predictors, response, weights) {
  cross <- crossprod(predictors, weights * predictors)
  if (rcond(cross) < .Machine$double.eps) return(NULL)
  drop(predictors %*% solve(cross, crossprod(predictors, weights * response)))
}
