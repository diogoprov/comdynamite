#' Regimes of the edges of a tree from a list of clades
#'
#' `clade_regimes()` assigns each edge of `tree` to a regime named after the clade it lies in,
#' for the `regimes` of [fit_damocles()] and [damocles_test()]: a partition of the tree defined
#' a priori, such as families or biogeographic groups, against which partitions from methods
#' that estimate shifts in rates of diversification can be compared.
#'
#' A clade is the crown group of the most recent common ancestor of its tips: the edges below
#' that node belong to the clade, and its stem edge too when `stem = TRUE`; a clade of one tip
#' is its terminal edge. If the tips given
#' are not monophyletic the clade covers every descendant of their common ancestor, with a
#' warning. Clades later in the list override earlier ones, so nested clades go after the
#' clades that hold them. Edges in no clade take the regime `background`.
#'
#' @param tree A rooted `phylo` tree.
#' @param clades A named list of character vectors, the tips of each clade.
#' @param background Name of the regime of the edges outside every clade.
#' @param stem If `TRUE`, the stem edge of each clade belongs to it.
#' @return A factor with one regime per edge of `tree`, in the order of `tree$edge`, with the
#'   regimes that occur as levels: `background` first, then the clades.
#' @seealso [fit_damocles()], [damocles_lrt()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(4)
#' tr <- treeshapesim::sim_tree_shape(20)
#' # the clade of the first node after the root
#' reg <- clade_regimes(tr, list(A = ape::extract.clade(tr, ape::Ntip(tr) + 2)$tip.label))
#' table(reg)
#' @export
clade_regimes <- function(tree, clades, background = "background", stem = FALSE) {
  if (!inherits(tree, "phylo")) stop("`tree` must be a `phylo` object", call. = FALSE)
  if (!is.list(clades) || length(clades) == 0 || is.null(names(clades)) ||
      any(names(clades) == "") || anyDuplicated(names(clades))) {
    stop("`clades` must be a list of tip labels with different names", call. = FALSE)
  }
  if (length(background) != 1 || background %in% names(clades)) {
    stop("`background` must be one name, different from the names of `clades`", call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  preorder <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  out <- rep(background, nrow(tree$edge))
  for (name in names(clades)) {
    tips <- clades[[name]]
    unknown <- setdiff(tips, tree$tip.label)
    if (length(unknown) > 0) {
      stop("tips of clade ", name, " not in `tree`: ", species_list(unknown), call. = FALSE)
    }
    node <- if (length(tips) == 1) match(tips, tree$tip.label) else ape::getMRCA(tree, tips)
    inside <- logical(n_tip + tree$Nnode)
    inside[node] <- TRUE
    for (e in preorder) {
      if (inside[tree$edge[e, 1]]) inside[tree$edge[e, 2]] <- TRUE
    }
    if (sum(inside[seq_len(n_tip)]) > length(unique(tips))) {
      warning("the tips of clade ", name, " are not monophyletic; it covers every ",
              "descendant of their common ancestor", call. = FALSE)
    }
    # a clade of one tip is its terminal edge
    edges <- inside[tree$edge[, 1]] | ((stem || node <= n_tip) & tree$edge[, 2] == node)
    out[edges] <- name
  }
  factor(out, levels = intersect(c(background, names(clades)), out))
}

#' Likelihood-ratio test between nested DAMOCLES fits
#'
#' `damocles_lrt()` compares two fits of [fit_damocles()] to the same communities, where `null`
#' is a special case of `alternative`: constant rates against rates that change among
#' `regimes` or with an `edge_covariate` (where `beta = 0` gives the constant model), rates
#' shared by all communities against separate rates, or `sympatry = 0` against an estimated
#' `sympatry`. Twice the difference of the total log-likelihoods is compared with a
#' chi-squared distribution with as many degrees of freedom as extra parameters.
#'
#' The function counts parameters but cannot check that the models are nested; that is up to
#' the user. When the null value lies on the boundary of the parameter space, as
#' `sympatry = 0` or a rate of zero in one regime, the chi-squared p-value is conservative.
#' The test answers whether the partition or covariate improves the fit to the local presences;
#' with partitions or covariates estimated from the same tree by methods of diversification,
#' a significant result is an association between those rates and the local dynamics, not
#' evidence of a mechanism.
#'
#' @param null,alternative Fits of class `mc_damocles_fit` to the same communities.
#' @return A data frame with one row per model: number of free parameters, total
#'   log-likelihood, AIC, and, for the alternative, the statistic, degrees of freedom and
#'   p-value.
#' @seealso [fit_damocles()], [clade_regimes()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(5)
#' tr <- treeshapesim::sim_tree_shape(40)
#' reg <- clade_regimes(tr, list(A = ape::extract.clade(tr, ape::Ntip(tr) + 2)$tip.label))
#' comm <- simulate_damocles(tr, gamma = 1, mu = 1, n = 5)
#' constant <- fit_damocles(comm, tr, shared = TRUE)
#' by_clade <- fit_damocles(comm, tr, shared = TRUE, regimes = reg, vary = "gamma")
#' damocles_lrt(constant, by_clade)
#' @export
damocles_lrt <- function(null, alternative) {
  if (!inherits(null, "mc_damocles_fit") || !inherits(alternative, "mc_damocles_fit")) {
    stop("`null` and `alternative` must be fits of `fit_damocles()`", call. = FALSE)
  }
  same <- identical(attr(null, "damocles")$presence, attr(alternative, "damocles")$presence)
  if (!same) stop("the two fits must be to the same communities", call. = FALSE)
  loglik <- c(damocles_total_loglik(null), damocles_total_loglik(alternative))
  n_par <- as.numeric(c(damocles_total_par(null), damocles_total_par(alternative)))
  df <- n_par[2] - n_par[1]
  if (df <= 0) {
    stop("`alternative` must have more free parameters than `null`", call. = FALSE)
  }
  statistic <- max(0, 2 * (loglik[2] - loglik[1]))
  data.frame(n_par = n_par, loglik = loglik, AIC = 2 * n_par - 2 * loglik,
             statistic = c(NA, statistic), df = c(NA, df),
             p_value = c(NA, stats::pchisq(statistic, df, lower.tail = FALSE)),
             row.names = c("null", "alternative"))
}

#' Profile-likelihood intervals of DAMOCLES parameters
#'
#' Confidence intervals for the parameters of [fit_damocles()] from the profile likelihood: the
#' interval holds the values at which the log-likelihood, maximised over the other parameters,
#' is within \eqn{\chi^2_1(\mathrm{level}) / 2}{qchisq(level, 1) / 2} of its maximum, half the
#' `level` quantile of a chi-squared distribution with one degree of freedom. With shared
#' parameters there is one interval
#' per parameter; otherwise one per community. A bound the profile does not reach within a
#' factor of \eqn{e^{10}}{exp(10)} of the estimate (rates), 10 units (slopes `beta`) or 10 units
#' on the logit
#' scale (`sympatry`) is reported as the edge of the parameter space (0 or `Inf` for rates,
#' `-Inf` or `Inf` for slopes, 0 or 1 for `sympatry`): the data do not bound it.
#'
#' @param object A fit of class `mc_damocles_fit`.
#' @param parm Names of the parameters, among the columns of the fit; by default all free
#'   parameters.
#' @param level Confidence level.
#' @param ... Not used.
#' @return A data frame with the community (`"shared"` for shared parameters), the parameter,
#'   its estimate and the lower and upper bounds.
#' @seealso [fit_damocles()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(6)
#' tr <- treeshapesim::sim_tree_shape(40)
#' comm <- simulate_damocles(tr, gamma = 1, mu = 1, n = 10)
#' fit <- fit_damocles(comm, tr, shared = TRUE)
#' confint(fit)
#' @exportS3Method stats::confint
confint.mc_damocles_fit <- function(object, parm = NULL, level = 0.95, ...) {
  info <- attr(object, "damocles")
  design <- info$design
  parm <- parm %||% design$names
  unknown <- setdiff(parm, design$names)
  if (length(unknown) > 0) {
    stop("not free parameters of the fit: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (length(level) != 1 || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one number between 0 and 1", call. = FALSE)
  }
  post <- damocles_postorder(info$tree)
  groups <- if (info$shared) list(seq_len(nrow(info$presence))) else
    as.list(seq_len(nrow(info$presence)))
  labels <- if (info$shared) "shared" else rownames(info$presence)
  drop <- stats::qchisq(level, 1) / 2
  rows <- list()
  for (g in seq_along(groups)) {
    objective <- damocles_objective(info$presence[groups[[g]], , drop = FALSE], info$tree,
                                    post, info$root, design)
    best <- info$par[[g]]
    top <- info$loglik[g]
    for (name in parm) {
      j <- match(name, design$names)
      bounds <- damocles_profile_bounds(objective, best, j, top - drop)
      rows[[length(rows) + 1]] <- data.frame(
        community = labels[g], parameter = name,
        estimate = damocles_natural(best[j], j, design),
        lower = damocles_natural(bounds[1], j, design),
        upper = damocles_natural(bounds[2], j, design))
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# --- internals ---------------------------------------------------------------------------

damocles_total_loglik <- function(fit) attr(fit, "loglik_total") %||% sum(fit$loglik)

damocles_total_par <- function(fit) {
  if (attr(fit, "damocles")$shared) attr(fit, "n_par") else attr(fit, "n_par") * nrow(fit)
}

# where the profile log-likelihood of working parameter `j` crosses `threshold`, on each side
# of the estimate; -Inf or Inf when it does not within 10 units
damocles_profile_bounds <- function(objective, best, j, threshold) {
  profile <- function(value) {
    fixed <- function(rest) objective(append(rest, value, after = j - 1))
    -damocles_minimise(fixed, best[-j])$value
  }
  gap <- function(value) profile(value) - threshold
  vapply(c(-1, 1), function(side) {
    inner <- best[j]
    for (step in seq_len(10)) {
      outer <- best[j] + side * step
      if (gap(outer) < 0) {
        return(stats::uniroot(gap, sort(c(inner, outer)), tol = 1e-6)$root)
      }
      inner <- outer
    }
    side * Inf
  }, numeric(1))
}
