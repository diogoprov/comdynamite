#' DAMOCLES fits over many configurations of regimes or covariates
#'
#' `damocles_samples()` fits [damocles_fit()], with parameters shared by all communities, once
#' for each configuration in a list of `regimes` or of `edge_covariate`, and compares each fit
#' with constant rates by [damocles_lrt()]. The configurations are typically samples from the
#' posterior distribution of a method that estimates rates of diversification along the tree,
#' so that the uncertainty in where rates shift, or in branch-specific rates, carries over to
#' the conclusions about local colonisation and extinction instead of resting on a single
#' configuration.
#'
#' With BAMM output read by `BAMMtools::getEventData()`, for example, the regimes of sample `i`
#' are `factor(edata$eventVectors[[i]])` and the mean branch rates of speciation are the columns
#' of `BAMMtools::getMarginalBranchRateMatrix(edata)$lambda_branch_matrix`, both in the order of
#' `tree$edge` when `edata` was built on `tree` (`data-raw/damocles_bamm.R` shows the whole
#' path). With ClaDS (`RPANDA::fit_ClaDS()`), the branch rates of one posterior sample are
#' `exp(l_0 + alpha * sampler$alpha_effect + x)`, where `l_0`, `alpha` and the values `x` of the
#' edges are columns of the chains, transformed as in `RPANDA::getMAPS_ClaDS()`
#' (`data-raw/damocles_clads.R`). Any other method works the same way once its output is a
#' vector per edge; [dr_rates()] gives one from the tree alone.
#'
#' Configurations of regimes and of continuous rates behave differently. A posterior sample of
#' regimes is one plausible partition, and the spread of the fits across samples shows how
#' much the conclusions depend on where the shifts are. A posterior sample of branch rates also
#' carries the sampling noise of each rate, and a noisy covariate attenuates its slope
#' (regression dilution): in `data-raw/damocles_clads.R`, with colonisation simulated from the
#' MAP rates of ClaDS with `beta = 2`, the MAP rates gave a slope of 1.6 and rejected constant
#' rates, while single posterior samples gave a median slope of 0.22 and rejected constant
#' rates in 65 of 100.
#'
#' `integrate = TRUE` adds one fit that uses all configurations at once: the parameters
#' maximise the likelihood averaged over the configurations, \eqn{L(\theta) = K^{-1} \sum_k
#' L(\theta \mid c_k)}{L(theta) = mean over k of L(theta | configuration k)}, so that the
#' uncertainty in the configuration enters the likelihood instead of a summary of separate
#' fits, and this fit is compared with constant rates (which do not depend on the
#' configuration) by a likelihood-ratio test. In simulations with posteriors of branch rates
#' from ClaDS (`data-raw/s3_clads_validation.R`), this test had the highest power of the ways
#' of summarising the samples, with the size of the test near or below its nominal level,
#' while its slope, like that of the posterior mean, was attenuated by the noise in the
#' rates. For continuous rates, then, use the integrated fit for the test, and report the
#' slope of the posterior mean (or median, or MAP) of each rate with [damocles_fit()] as a
#' lower bound of the effect, knowing that single samples understate it further. The
#' parameters must mean the same in every configuration: always true for `edge_covariate`,
#' and true for `regimes` only when every configuration has the same regimes with the same
#' meaning (for example, `"0"` for the background and `"1"` for the clade below a shift), not
#' for the arbitrary labels of the regimes of each BAMM sample.
#'
#' `summary()` reports, across configurations, how often the varying rates fit better than
#' constant rates (p-value below 0.05, and AIC lower by more than 2), quantiles of the
#' parameters that mean the same in every fit (such as `mu` when only `gamma` varies, or
#' `beta_gamma`),
#' and quantiles of the rate of each edge, which are comparable whatever the regimes of each
#' sample. Regimes differ among posterior samples, so the parameters of each regime are only
#' comparable within a sample.
#'
#' The configurations are only as good as the method that produced them. The reliability of
#' the rate shifts inferred by BAMM has been debated (Moore et al. 2016; Rabosky et al. 2017);
#' comparing several sources of regimes, including partitions defined a priori with
#' [clade_regimes()], and testing a configuration against random partitions with
#' [damocles_partition_test()] guard against relying on one of them.
#'
#' @inheritParams damocles_fit
#' @param regimes,edge_covariate A list of configurations, each a vector with one value per edge
#'   of `tree` in the order of `tree$edge` (see [damocles_fit()]). Give one of the two.
#' @param integrate `TRUE` or `FALSE` (the default). `TRUE` also fits the parameters that
#'   maximise the likelihood averaged over the configurations and tests them against constant
#'   rates (see Details).
#' @param cores Number of cores, with configurations fitted in parallel. One whole number, 1 or
#'   more. Not available on Windows, where one core is used.
#' @return A list of class `mc_damocles_samples` with `samples`, a data frame with one row per
#'   configuration (`n_par`, `loglik`, `delta_aic` against constant rates, and the
#'   likelihood-ratio `statistic`, `df` and `p_value`), `constant`, the fit with constant rates,
#'   `fits`, the fit of each configuration, `gamma` and `mu`, matrices of the rate of each
#'   edge (rows) in each configuration (columns), and `integrated`, `NULL` unless `integrate =
#'   TRUE`, when it is a list with the `estimates` of the integrated fit, its `loglik` and
#'   `n_par`, `delta_aic` against constant rates, the likelihood-ratio `statistic`, `df` and
#'   `p_value`, and the `convergence` code of the optimiser.
#' @references Moore, B. R., Höhna, S., May, M. R., Rannala, B., & Huelsenbeck, J. P. (2016).
#'   Critically evaluating the theory and performance of Bayesian analysis of macroevolutionary
#'   mixtures. *Proceedings of the National Academy of Sciences*, 113, 9569--9574.
#'   \doi{10.1073/pnas.1518659113}
#'
#'   Rabosky, D. L., Mitchell, J. S., & Chang, J. (2017). Is BAMM flawed? Theoretical and
#'   practical concerns in the analysis of multi-rate diversification models. *Systematic
#'   Biology*, 66, 477--498. \doi{10.1093/sysbio/syx037}
#' @seealso [damocles_fit()], [damocles_lrt()], [damocles_partition_test()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(7)
#' tr <- treeshapesim::sim_tree_shape(30)
#' crown <- ape::Ntip(tr) + 2
#' comm <- simulate_damocles(tr, ifelse(clade_regimes(tr, list(A = crown)) == "A", 3, 0.5),
#'                           mu = 1, n = 5)
#' # three configurations, as if sampled from a posterior
#' configurations <- lapply(c(crown, crown + 1, crown + 2), function(node) {
#'   clade_regimes(tr, list(A = node))
#' })
#' res <- damocles_samples(comm, tr, regimes = configurations, vary = "gamma",
#'                         integrate = TRUE)
#' res
#' summary(res)
#' @export
damocles_samples <- function(comm, tree, regimes = NULL, edge_covariate = NULL,
                             vary = c("both", "gamma", "mu"),
                             root = c("both", "absent", "present"), sympatry = 0,
                             integrate = FALSE, cores = 1) {
  vary <- match.arg(vary)
  root <- match.arg(root)
  check_flag(integrate, "integrate")
  check_count(cores, "cores")
  if (is.null(regimes) == is.null(edge_covariate)) {
    stop("give a list of `regimes` or a list of `edge_covariate`", call. = FALSE)
  }
  configurations <- regimes %||% edge_covariate
  if (!is.list(configurations) || length(configurations) == 0) {
    stop("`regimes` or `edge_covariate` must be a list of configurations", call. = FALSE)
  }
  check_damocles_tree(tree)
  presence <- damocles_presence(comm, tree)
  constant <- damocles_fit(presence, tree, root = root, shared = TRUE, sympatry = sympatry)
  run <- function(i) {
    damocles_fit(presence, tree, root = root, shared = TRUE, sympatry = sympatry,
                 regimes = if (!is.null(regimes)) configurations[[i]],
                 edge_covariate = if (!is.null(edge_covariate)) configurations[[i]],
                 vary = vary)
  }
  fits <- damocles_apply(seq_along(configurations), run, cores)
  rows <- lapply(fits, function(fit) {
    if (attr(fit, "n_par") <= attr(constant, "n_par")) {
      return(data.frame(n_par = attr(fit, "n_par"), loglik = attr(fit, "loglik_total"),
                        delta_aic = 0, statistic = NA_real_, df = NA_real_,
                        p_value = NA_real_))
    }
    lrt <- damocles_lrt(constant, fit)
    data.frame(n_par = lrt$n_par[2], loglik = lrt$loglik[2],
               delta_aic = lrt$AIC[2] - lrt$AIC[1], statistic = lrt$statistic[2],
               df = lrt$df[2], p_value = lrt$p_value[2])
  })
  samples <- cbind(sample = seq_along(fits), do.call(rbind, rows))
  rates <- function(rate) {
    vapply(fits, function(fit) attr(fit, "edge_rates")[[1]][[rate]], numeric(nrow(tree$edge)))
  }
  integrated <- if (integrate) {
    damocles_integrated(presence, tree, configurations, !is.null(regimes), vary, root,
                        sympatry, constant)
  }
  structure(list(samples = samples, constant = constant, fits = fits,
                 gamma = rates("gamma"), mu = rates("mu"), integrated = integrated),
            class = "mc_damocles_samples")
}

#' log of the mean of exp(ll), with values at or below -1e10 (failed likelihoods) taken as
#' -Inf, and -1e10 when none is finite
#' @noRd
log_mean_exp <- function(ll) {
  ll[ll <= -1e10] <- -Inf
  if (!any(is.finite(ll))) return(-1e10)
  top <- max(ll)
  top + log(mean(exp(ll - top)))
}

#' fit maximising the likelihood averaged over the configurations, L(theta) = mean over k of
#' L(theta | configuration k), and its likelihood-ratio test against constant rates
#' @noRd
damocles_integrated <- function(presence, tree, configurations, is_regimes, vary, root,
                                sympatry, constant) {
  designs <- lapply(configurations, function(cfg) {
    damocles_design(tree, if (is_regimes) cfg, if (!is_regimes) cfg, vary, sympatry)
  })
  if (length(unique(lapply(designs, `[[`, "names"))) > 1) {
    stop("`integrate = TRUE` in `damocles_samples()` needs configurations with the same ",
         "parameters: the same regimes in every configuration", call. = FALSE)
  }
  design <- designs[[1]]
  post <- damocles_postorder(tree)
  objectives <- lapply(designs, function(d) damocles_objective(presence, tree, post, root, d))
  objective <- function(par) {
    -log_mean_exp(-vapply(objectives, function(o) o(par), numeric(1)))
  }
  base <- rbind(c(constant$gamma[1], constant$mu[1]), damocles_starts(post$tree, NULL))
  best <- damocles_optim(objective, damocles_par_starts(base, design))
  loglik <- -best$value
  loglik[loglik <= -1e10] <- -Inf
  loglik0 <- attr(constant, "loglik_total")
  df <- as.numeric(design$n_par - attr(constant, "n_par"))
  statistic <- if (df > 0) max(2 * (loglik - loglik0), 0) else NA_real_
  list(estimates = damocles_par_values(best$par, design), loglik = loglik,
       n_par = design$n_par, delta_aic = 2 * df - 2 * (loglik - loglik0),
       statistic = statistic, df = if (df > 0) df else NA_real_,
       p_value = if (df > 0) stats::pchisq(statistic, df, lower.tail = FALSE) else NA_real_,
       convergence = best$convergence)
}

#' @export
print.mc_damocles_samples <- function(x, ...) {
  s <- x$samples
  cat("DAMOCLES fits over", nrow(s), "configurations, parameters shared by",
      nrow(attr(x$constant, "damocles")$presence), "communities\n")
  cat("Better than constant rates: p < 0.05 in", sum(s$p_value < 0.05, na.rm = TRUE),
      "and AIC lower by more than 2 in", sum(s$delta_aic < -2), "of", nrow(s), "\n")
  if (!is.null(x$integrated)) {
    cat("Integrated over configurations: likelihood-ratio statistic",
        format(x$integrated$statistic, digits = 3), "on", x$integrated$df, "df, p =",
        format(x$integrated$p_value, digits = 3), "\n")
  }
  invisible(x)
}

#' @rdname damocles_samples
#' @param object A result of `damocles_samples()`.
#' @param probs Probabilities of the quantiles.
#' @param ... Not used.
#' @export
summary.mc_damocles_samples <- function(object, probs = c(0.05, 0.5, 0.95), ...) {
  s <- object$samples
  quantiles <- function(x) stats::quantile(x, probs, na.rm = TRUE, names = TRUE)
  comparison <- data.frame(
    measure = c("p_value < 0.05", "delta_aic < -2"),
    proportion = c(mean(!is.na(s$p_value) & s$p_value < 0.05), mean(s$delta_aic < -2)))
  # parameters of a regime are named after it but differ among configurations
  shared <- intersect(Reduce(intersect, lapply(object$fits, attr, "parameters")),
                      c("gamma", "mu", "beta_gamma", "beta_mu"))
  parameters <- if (length(shared) > 0) {
    t(vapply(shared, function(name) {
      quantiles(vapply(object$fits, function(fit) fit[[name]][1], numeric(1)))
    }, numeric(length(probs))))
  }
  edges <- data.frame(edge = seq_len(nrow(object$gamma)),
                      gamma = t(apply(object$gamma, 1, quantiles)),
                      mu = t(apply(object$mu, 1, quantiles)), check.names = FALSE)
  structure(list(comparison = comparison, statistic = quantiles(s$statistic),
                 delta_aic = quantiles(s$delta_aic), parameters = parameters, edges = edges,
                 n = nrow(s), integrated = object$integrated),
            class = "summary.mc_damocles_samples")
}

#' @export
print.summary.mc_damocles_samples <- function(x, ...) {
  cat("DAMOCLES fits over", x$n, "configurations\n\n")
  cat("Proportion of configurations better than constant rates:\n")
  print(x$comparison, row.names = FALSE, digits = 3)
  cat("\nLikelihood-ratio statistic:\n")
  print(x$statistic, digits = 3)
  cat("Difference in AIC (varying - constant):\n")
  print(x$delta_aic, digits = 3)
  if (!is.null(x$parameters)) {
    cat("\nParameters shared by every fit:\n")
    print(x$parameters, digits = 3)
  }
  if (!is.null(x$integrated)) {
    i <- x$integrated
    cat("\nIntegrated over configurations (log-likelihood ", format(i$loglik, digits = 6),
        "):\n", sep = "")
    print(i$estimates, digits = 3)
    cat("Against constant rates: statistic", format(i$statistic, digits = 3), "on", i$df,
        "df, p =", format(i$p_value, digits = 3), "\n")
  }
  cat("\nRates of each edge: see `$edges`\n")
  invisible(x)
}

#' Random partitions of a tree like a given one
#'
#' `random_regimes()` draws a partition of the edges of `tree` with the same regimes as
#' `regimes` and clades of similar size, placed at random: the null partitions of
#' [damocles_partition_test()]. Each regime other than that of the root is taken as one or more
#' clades, each starting at an edge where the regime changes (as a shift on that branch in
#' BAMM, or a clade with its stem in [clade_regimes()]); each clade is moved to an edge, drawn at
#' random, whose subtree holds between `1 / tolerance` and `tolerance` times as many edges, from
#' the largest clade to the smallest, so that nested clades stay nested when sizes allow it.
#' Draws in which a later clade covers an earlier one completely are discarded.
#'
#' Edges whose subtree largely coincides with the original clade, such as the edge just above
#' it or those of its large daughter clades, give partitions that are nearly the original one.
#' `max_overlap` excludes them: the Jaccard index between the edges of the original clade and
#' those of its replacement (shared edges over edges in either) must not exceed it. With
#' `max_overlap = 1` any edge of the right size can be drawn.
#'
#' A clade can also be replaced by its complement: when the clade of a regime is one of the two
#' clades at the root, the other one has no edge in common with it but splits the tree in the
#' same way, and the random partition is the original one with the regimes swapped (in trees with
#' a shift in the rate of speciation, the faster clade often holds most of the tree). So
#' `max_overlap` also bounds the agreement between the whole partitions, measured by the
#' adjusted Rand index between the regimes of the edges (1 for the same partition, whatever the
#' labels; about 0 for unrelated ones); draws above it are discarded.
#'
#' @param tree A rooted `phylo` tree.
#' @param regimes Regime of each edge, in the order of `tree$edge`.
#' @param tolerance Largest ratio between the sizes of a clade and its random replacement. One
#'   number, 1 or more.
#' @param max_overlap Largest Jaccard index between the edges of a clade and those of its
#'   random replacement, and largest adjusted Rand index between `regimes` and the random
#'   partition. One number between 0 and 1.
#' @return A factor with one regime per edge, with the levels of `regimes`.
#' @seealso [damocles_partition_test()], [clade_regimes()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(8)
#' tr <- treeshapesim::sim_tree_shape(30)
#' reg <- clade_regimes(tr, list(A = ape::Ntip(tr) + 2))
#' table(reg)
#' table(random_regimes(tr, reg))
#' @export
random_regimes <- function(tree, regimes, tolerance = 2, max_overlap = 0.25) {
  check_tree(tree, branch_lengths = FALSE)
  n_edge <- nrow(tree$edge)
  if (length(regimes) != n_edge || anyNA(regimes)) {
    stop("random_regimes(): `regimes` needs one regime for each edge of `tree`, in the order ",
         "of `tree$edge`, ",
         "with no missing values", call. = FALSE)
  }
  check_number(tolerance, "tolerance", lower = 1)
  check_number(max_overlap, "max_overlap", 0, 1)
  regimes <- droplevels(factor(regimes))
  labels <- as.character(regimes)
  parent_edge <- match(tree$edge[, 1], tree$edge[, 2])
  at_root <- is.na(parent_edge)
  background <- names(which.max(table(labels[at_root])))
  above <- ifelse(at_root, background, labels[parent_edge])
  starts <- which(labels != above)
  size <- damocles_subtree_size(tree)
  starts <- starts[order(size[starts], decreasing = TRUE)]
  subtree <- damocles_subtree_edges(tree)
  # edges each clade may move to: not overlapping the original clade more than allowed (the
  # subtree of its sister edge never does, so some edge is always allowed)
  allowed <- lapply(starts, function(start) {
    shared <- colSums(subtree & subtree[, start])
    shared / (size + size[start] - shared) <= max_overlap + 1e-12
  })
  for (attempt in seq_len(100)) {
    out <- rep(background, n_edge)
    used <- logical(n_edge)
    for (k in seq_along(starts)) {
      start <- starts[k]
      ratio <- abs(log(size / size[start]))
      ratio[used | !allowed[[k]]] <- Inf
      candidates <- which(ratio <= log(tolerance))
      if (length(candidates) == 0) candidates <- which(ratio == min(ratio))
      pick <- candidates[sample.int(length(candidates), 1)]
      used[pick] <- TRUE
      out[subtree[, pick]] <- labels[start]
    }
    if (setequal(unique(out), levels(regimes)) &&
        adjusted_rand(out, labels) <= max_overlap + 1e-12) {
      return(factor(out, levels = levels(regimes)))
    }
  }
  stop("could not place every regime without one covering another or the partition matching ",
       "`regimes`; increase `tolerance` or `max_overlap`", call. = FALSE)
}

#' adjusted Rand index (Hubert & Arabie 1985, J. Classif. 2: 193-218, doi:10.1007/BF01908075)
#' between two partitions of the same items
#' @noRd
adjusted_rand <- function(a, b) {
  pairs <- function(x) sum(x * (x - 1) / 2)
  counts <- table(a, b)
  index <- pairs(counts)
  rows <- pairs(rowSums(counts))
  cols <- pairs(colSums(counts))
  expected <- rows * cols / pairs(length(a))
  top <- (rows + cols) / 2
  if (top == expected) return(1)
  (index - expected) / (top - expected)
}

#' Test a partition of the tree against random partitions
#'
#' `damocles_partition_test()` asks whether rates that change among `regimes` fit the local
#' presences better than rates that change among partitions of the same complexity placed at
#' random on the tree ([random_regimes()]). The statistic is the likelihood ratio of the fit
#' with `regimes` against constant rates, with parameters shared by all communities
#' ([damocles_lrt()]); its p-value is the proportion of random partitions with a statistic as
#' large or larger.
#'
#' The chi-squared test of [damocles_lrt()] asks whether rates vary at all along the tree in the
#' way the partition allows. When they vary somewhere, any partition captures part of that and
#' the chi-squared test rejects constant rates even for a partition unrelated to where the rates
#' change. The random partitions keep the number of parameters and the sizes of the clades, so
#' this test asks instead whether the partition given is better placed than arbitrary ones:
#' whether the regimes from a method of diversification, or a partition defined a priori,
#' carry information on local colonisation and extinction beyond rate heterogeneity in general.
#' `data-raw/validation_damocles.R` shows the difference between the two tests: with 40
#' replicates of 10 communities on trees of 60 tips and a clade holding about a quarter of the
#' edges, the chi-squared test rejected constant rates in 3 of 40 without an effect, in 40 of
#' 40 with colonisation faster in the clade, and in 29 of 40 with it faster in another clade of
#' similar size; the test against random partitions rejected in 4, 37 and 0 of 40. The 4 of 40
#' without an effect is within sampling error of 5% but may be slightly liberal.
#'
#' Random partitions that nearly coincide with the one tested, such as the clade one edge
#' above it, capture the same rates and leave the test without power: with `max_overlap = 1`
#' the test against random partitions rejected in none of 40 replicates with the effect in the
#' clade. The default `max_overlap = 0.25` keeps only placements elsewhere on the tree, so the
#' test asks whether the partition beats partitions of the same complexity placed somewhere
#' else. With the whales and the best BAMM configuration of BAMMtools (a regime on 69 of 172
#' edges) and colonisation three times faster in it, both tests rejected constant rates (p =
#' 0.007 and 0.01; `data-raw/damocles_bamm.R`).
#'
#' The test needs room on the tree for clades like those of `regimes` elsewhere. When the
#' partition tested is one of the clades at the root, or a clade holding most of the tree, the
#' only placements allowed may split the tree as `regimes` does (its complement) or be few, and
#' the random statistics then repeat a handful of values: the p-value falls at its smallest or
#' at one depending on numerical noise. In simulations with one shift in the rate of speciation
#' and no effect on local rates (`data-raw/s1_type1_shifts.R`), this made the test reject in 12%
#' to 40% of replicates before random partitions matching `regimes` were excluded
#' ([random_regimes()]) and before the p-value counted each distinct partition once. The
#' `n_rand` draws are a way to find distinct placements: the p-value is the proportion of the
#' distinct random partitions, plus the one tested, with a statistic as large or larger, so with
#' fewer than 19 distinct ones it cannot fall below 0.05 (a warning says so). Trees with a strong
#' shift in diversification often leave that little room. In the same simulations after both
#' changes (50 replicates per setting, 99 draws), replicates with one shift had a median of 3
#' distinct random partitions and 79 of 500 had none; replicates with 19 or more rejected in 8%
#' (37 of 452).
#'
#' `n_boot` adds a parametric bootstrap of the same likelihood ratio: communities are simulated
#' with the constant rates estimated from `comm` ([simulate_damocles()]), as many as in `comm`,
#' and each set is fitted with and without `regimes`; `p_boot` is the proportion of simulated
#' statistics, plus the observed one, as large or larger. It asks the question of the
#' chi-squared test (do rates vary along the tree as the partition allows?) without relying on
#' the asymptotic distribution of the likelihood ratio, which is doubtful with few communities
#' and rates near zero, and it needs no room on the tree, so it is available for any partition.
#' It does not ask whether the partition is better placed than others: with rate heterogeneity
#' elsewhere on the tree it rejects, like the chi-squared test. Each simulation costs two fits,
#' twice a random partition. Set `n_rand = 0` to run the bootstrap alone.
#'
#' @inheritParams damocles_fit
#' @param regimes Regime of each edge of `tree`, in the order of `tree$edge`, with two or more
#'   regimes.
#' @param n_rand Number of random partitions drawn; each distinct one is fitted once. One whole
#'   number, 0 or more; 0 skips the test against random partitions.
#' @param n_boot Number of parametric bootstrap simulations under the constant rates estimated.
#'   One whole number, 0 or more; 0, the default, skips the bootstrap.
#' @param tolerance,max_overlap Passed to [random_regimes()]: random clades hold between
#'   `1 / tolerance` and `tolerance` times as many edges as the clades of `regimes`, share
#'   at most `max_overlap` of their edges with them (Jaccard index), and the whole random
#'   partition agrees with `regimes` at most as much (adjusted Rand index). `tolerance` is one
#'   number, 1 or more; `max_overlap` one number between 0 and 1.
#' @param cores Number of cores, with random partitions and bootstrap simulations fitted in
#'   parallel. One whole number, 1 or more. Not available on Windows, where one core is used.
#'   The random partitions and the simulated communities are drawn before the fits, so the
#'   p-values do not depend on the number of cores.
#' @return A list of class `mc_damocles_partition` with the observed `statistic`, its degrees
#'   of freedom `df`, the chi-squared p-value `p_chisq`, the p-value against random partitions
#'   `p_random` (`NA` when no random partition could be drawn), the statistics of the distinct
#'   random partitions `null`, their number `n_distinct`, `n_rand`, the p-value of the
#'   parametric bootstrap `p_boot` (`NA` without it), the simulated statistics `null_boot`, and
#'   `fit`, the fit with `regimes`.
#' @seealso [random_regimes()], [damocles_samples()], [damocles_lrt()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(9)
#' tr <- treeshapesim::sim_tree_shape(30)
#' reg <- clade_regimes(tr, list(A = ape::Ntip(tr) + 2))
#' comm <- simulate_damocles(tr, ifelse(reg == "A", 3, 0.5), mu = 1, n = 3)
#' # few sites and random partitions, to keep the example quick
#' damocles_partition_test(comm, tr, reg, n_rand = 2, vary = "gamma")
#' \donttest{
#' # random partitions and the parametric bootstrap refit the model many times
#' comm <- simulate_damocles(tr, ifelse(reg == "A", 3, 0.5), mu = 1, n = 5)
#' damocles_partition_test(comm, tr, reg, n_rand = 19, vary = "gamma")
#' damocles_partition_test(comm, tr, reg, n_rand = 0, n_boot = 19, vary = "gamma")
#' }
#' @export
damocles_partition_test <- function(comm, tree, regimes, n_rand = 99, n_boot = 0,
                                    vary = c("both", "gamma", "mu"),
                                    root = c("both", "absent", "present"), sympatry = 0,
                                    tolerance = 2, max_overlap = 0.25, cores = 1) {
  vary <- match.arg(vary)
  root <- match.arg(root)
  check_count(n_rand, "n_rand", min = 0)
  check_count(n_boot, "n_boot", min = 0)
  if (n_rand + n_boot == 0) {
    stop("`n_rand` or `n_boot` must be 1 or more", call. = FALSE)
  }
  # random_regimes() runs inside tryCatch() below, so its arguments are checked here
  check_number(tolerance, "tolerance", lower = 1)
  check_number(max_overlap, "max_overlap", 0, 1)
  check_count(cores, "cores")
  check_damocles_tree(tree)
  if (nlevels(droplevels(factor(regimes))) < 2) {
    stop("`regimes` needs two or more regimes", call. = FALSE)
  }
  presence <- damocles_presence(comm, tree)
  constant <- damocles_fit(presence, tree, root = root, shared = TRUE, sympatry = sympatry)
  fit_with <- function(reg) {
    damocles_fit(presence, tree, root = root, shared = TRUE, sympatry = sympatry,
                 regimes = reg, vary = vary)
  }
  fit <- fit_with(regimes)
  observed <- damocles_lrt(constant, fit)
  random <- partition_random(tree, regimes, n_rand, tolerance, max_overlap, constant, fit_with,
                             observed$statistic[2], cores)
  boot <- partition_bootstrap(presence, tree, regimes, n_boot, vary, root, sympatry, constant,
                              observed$statistic[2], cores)
  structure(list(statistic = observed$statistic[2], df = observed$df[2],
                 p_chisq = observed$p_value[2], p_random = random$p, null = random$null,
                 n_distinct = length(random$null), n_rand = n_rand, p_boot = boot$p,
                 null_boot = boot$null,
                 fit = fit),
            class = "mc_damocles_partition")
}

#' the test of damocles_partition_test() against random partitions: a list with the
#' likelihood-ratio statistics of the distinct ones (`null`) and the p-value of `statistic`
#' (`p`, NA when none could be drawn). `fit_with` fits the model with a partition
#' @noRd
partition_random <- function(tree, regimes, n_rand, tolerance, max_overlap, constant, fit_with,
                             statistic, cores) {
  # Draw the random partitions first, then fit each distinct one once: the p-value does not
  # depend on the number of cores.
  draws <- lapply(seq_len(n_rand), function(i) {
    tryCatch(random_regimes(tree, regimes, tolerance, max_overlap), error = function(e) e)
  })
  failed <- vapply(draws, inherits, logical(1), what = "error")
  if (n_rand > 0 && all(failed)) {
    warning("no random partition could be drawn (", conditionMessage(draws[[1]]), "); ",
            "`p_random` is NA", call. = FALSE)
  }
  draws <- draws[!failed]
  keys <- vapply(draws, function(r) paste(as.integer(r), collapse = ","), character(1))
  distinct <- unique(keys)
  stats <- unlist(damocles_apply(match(distinct, keys), function(i) {
    damocles_lrt(constant, fit_with(draws[[i]]))$statistic[2]
  }, cores))
  null <- if (length(draws)) unname(stats) else numeric(0)
  # Each distinct partition counts once: repeated draws of a few partitions would otherwise let
  # the p-value fall below what the number of distinct alternatives allows.
  p_random <- if (length(null)) {
    (sum(null >= statistic - 1e-8) + 1) / (length(null) + 1)
  } else {
    NA_real_
  }
  if (length(null) && length(null) < 19) {
    warning("only ", length(null), " distinct random partitions in ", n_rand, " draws: the ",
            "tree has few places for clades like those of `regimes`, and `p_random` cannot be ",
            "below ", format(1 / (length(null) + 1), digits = 2), call. = FALSE)
  }
  list(null = null, p = p_random)
}

#' the parametric bootstrap of damocles_partition_test(): communities simulated with the
#' `constant` rates estimated, refitted with and without `regimes`. A list with the simulated
#' statistics (`null`) and the p-value of `statistic` (`p`, NA without simulations)
#' @noRd
partition_bootstrap <- function(presence, tree, regimes, n_boot, vary, root, sympatry,
                                constant, statistic, cores) {
  if (n_boot == 0) return(list(null = numeric(0), p = NA_real_))
  # All communities are simulated before the fits, so the p-value does not depend on the
  # number of cores.
  sims <- lapply(seq_len(n_boot), function(i) {
    simulate_damocles(tree, constant$gamma[1], constant$mu[1], n = nrow(presence),
                      root = root, sympatry = constant$sympatry[1])
  })
  null_boot <- unlist(damocles_apply(sims, function(sim) {
    tryCatch({
      base <- damocles_fit(sim, tree, root = root, shared = TRUE, sympatry = sympatry)
      alt <- damocles_fit(sim, tree, root = root, shared = TRUE, sympatry = sympatry,
                          regimes = regimes, vary = vary)
      damocles_lrt(base, alt)$statistic[2]
    }, error = function(e) NA_real_)
  }, cores))
  null_boot <- unname(null_boot[is.finite(null_boot)])
  if (length(null_boot) < n_boot) {
    warning(n_boot - length(null_boot), " of ", n_boot, " bootstrap fits failed and were ",
            "left out", call. = FALSE)
  }
  p_boot <- NA_real_
  if (length(null_boot)) {
    p_boot <- (sum(null_boot >= statistic - 1e-8) + 1) / (length(null_boot) + 1)
  }
  list(null = null_boot, p = p_boot)
}

#' @export
print.mc_damocles_partition <- function(x, ...) {
  cat("Partition of the tree against constant rates (DAMOCLES)\n")
  cat("Likelihood ratio against constant rates:", format(x$statistic, digits = 4), "on",
      x$df, "df\n")
  cat("p (chi-squared):", format(x$p_chisq, digits = 3), "\n")
  if ((x$n_rand %||% 1) > 0) {
    cat("p (random partitions): ", format(x$p_random, digits = 3), " (", x$n_distinct,
        " distinct random partitions)\n", sep = "")
  }
  if (!is.null(x$p_boot) && (length(x$null_boot) || !is.na(x$p_boot))) {
    cat("p (parametric bootstrap): ", format(x$p_boot, digits = 3), " (",
        length(x$null_boot), " simulations under constant rates)\n", sep = "")
  }
  invisible(x)
}

# --- internals ---------------------------------------------------------------------------

#' number of edges in the subtree of each edge, the edge included
#' @noRd
damocles_subtree_size <- function(tree) {
  below <- numeric(ape::Ntip(tree) + tree$Nnode)
  for (e in ape::reorder.phylo(tree, "postorder", index.only = TRUE)) {
    child <- tree$edge[e, 2]
    below[tree$edge[e, 1]] <- below[tree$edge[e, 1]] + 1 + below[child]
  }
  1 + below[tree$edge[, 2]]
}

#' edges (rows) in the subtree of each edge (columns), the edge included
#' @noRd
damocles_subtree_edges <- function(tree) {
  n_edge <- nrow(tree$edge)
  out <- diag(n_edge) == 1
  # postorder: an edge's subtree is itself plus the subtrees of the edges below it
  child_edges <- split(seq_len(n_edge), tree$edge[, 1])
  for (e in ape::reorder.phylo(tree, "postorder", index.only = TRUE)) {
    below <- child_edges[[as.character(tree$edge[e, 2])]]
    if (!is.null(below)) out[, e] <- out[, e] | rowSums(out[, below, drop = FALSE]) > 0
  }
  out
}

damocles_apply <- function(x, f, cores) {
  check_count(cores, "cores")
  workers <- if (.Platform$OS.type == "windows") 1L else as.integer(cores)
  if (workers > 1) parallel::mclapply(x, f, mc.cores = workers) else lapply(x, f)
}
