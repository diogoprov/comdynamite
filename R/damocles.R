#' Local presence along a tree: the dynamic null model of DAMOCLES
#'
#' `simulate_damocles()` simulates the presence and absence of each lineage of a clade in one
#' local community, as it diversifies along `tree`, under the dynamic null model of community
#' assembly of Pigot & Etienne (2015), DAMOCLES. Along each branch the lineage is either locally
#' absent (state 0) or present (state 1); an absent lineage colonises the community at rate
#' `gamma` and a present one goes locally extinct at rate `mu`, the same for every lineage
#' unless rates are given for each edge. At each speciation, which is allopatric, a locally
#' present ancestor leaves one daughter present and the other absent, chosen at random, and an absent ancestor leaves two absent daughters.
#' The states at the tips form the community.
#'
#' `sympatry` relaxes the rule of strictly allopatric speciation: with that probability a
#' present ancestor leaves both daughters present, as after speciation within the community
#' (in situ or sympatric speciation) or when one daughter recolonises it before the two can be
#' told apart; otherwise one daughter is present, as above. `sympatry = 0`, the default, is the
#' model of Pigot & Etienne (2015). An absent ancestor always leaves two absent daughters.
#'
#' The states at the nodes are drawn from the exact transition probabilities of the two-state
#' process over each branch, which gives the same distribution of tip states as simulating each
#' colonisation and extinction event. The model assumes the tree holds the species of the
#' regional pool, so it is meant for ultrametric trees of extant species. Branches of length
#' zero and a root edge (`tree$root.edge`) are allowed.
#'
#' @param tree A rooted, fully dichotomous `phylo` tree with branch lengths.
#' @param gamma Rate of local colonisation per lineage, in the units of the branch lengths: one
#'   rate, or one per edge of `tree` in the order of `tree$edge`.
#' @param mu Rate of local extinction per lineage, one or one per edge. With rates per edge and
#'   a root edge, the two edges leaving the root must have the same rates, which the root edge
#'   takes.
#' @param sympatry Probability, between 0 and 1, that a locally present ancestor leaves both
#'   daughters present at speciation; 0 (the default) is strictly allopatric speciation.
#' @param n Number of communities to simulate.
#' @param root State of the ancestor at the root, just before the first speciation (or at the
#'   start of the root edge): `"both"` draws absent or present with equal probabilities for each
#'   community, as in the simulations of Pigot & Etienne (2015); `"absent"` or `"present"` fix
#'   it.
#' @return A matrix of zeros and ones with one row per community and one column per tip, named
#'   by the tips, ready to be used as `comm` in the other functions of the package.
#' @references Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for phylogenetic
#'   community structure. *Ecology Letters*, 18, 153--163. \doi{10.1111/ele.12395}
#' @seealso [damocles_loglik()] and [fit_damocles()] for the likelihood, [damocles_test()] for
#'   the null model of phylogenetic community structure.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(30)
#' comm <- simulate_damocles(tr, gamma = 0.5, mu = 0.5, n = 5)
#' rowSums(comm)
#' @export
simulate_damocles <- function(tree, gamma, mu, n = 1, root = c("both", "absent", "present"),
                              sympatry = 0) {
  root <- match.arg(root)
  check_damocles_tree(tree)
  rates <- damocles_edge_rates(tree, gamma, mu)
  check_sympatry(sympatry)
  if (length(n) != 1 || !is.finite(n) || n < 1 || n != round(n)) {
    stop("`n` must be one whole number, 1 or more", call. = FALSE)
  }
  index <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  tree <- ape::reorder.phylo(tree, "cladewise")
  gamma <- rates$gamma[index]
  mu <- rates$mu[index]
  n_tip <- ape::Ntip(tree)
  start <- switch(root, both = stats::rbinom(n, 1, 0.5), absent = rep(0L, n),
                  present = rep(1L, n))
  if (!is.null(tree$root.edge)) {
    start <- damocles_step(start, tree$root.edge, rates$root[[1]], rates$root[[2]])
  }
  # state at the end of the branch leading to each node, just before it speciates
  state <- matrix(0L, n, n_tip + tree$Nnode)
  state[, n_tip + 1] <- start
  # which of the two daughters inherits the presence of its ancestor (3: both of them)
  heir <- matrix(sample.int(2, n * tree$Nnode, replace = TRUE), n, tree$Nnode)
  if (sympatry > 0) heir[stats::runif(length(heir)) < sympatry] <- 3L
  position <- stats::ave(seq_len(nrow(tree$edge)), tree$edge[, 1], FUN = seq_along)
  for (e in seq_len(nrow(tree$edge))) {  # cladewise: ancestors before descendants
    parent <- tree$edge[e, 1]
    inherits <- heir[, parent - n_tip] == position[e] | heir[, parent - n_tip] == 3L
    begin <- state[, parent] * inherits
    state[, tree$edge[e, 2]] <- damocles_step(begin, tree$edge.length[e], gamma[e], mu[e])
  }
  out <- state[, seq_len(n_tip), drop = FALSE]
  dimnames(out) <- list(paste0("sim", seq_len(n)), tree$tip.label)
  out
}

#' Likelihood of local presences under the DAMOCLES model
#'
#' `damocles_loglik()` computes the exact log-likelihood of the presences and absences of the
#' species of a clade in a local community, given their tree and the rates of local colonisation
#' `gamma` and extinction `mu`, under the model of [simulate_damocles()]. It follows the pruning
#' algorithm of Pigot & Etienne (2015, Box 1): from the tips to the root, the probabilities of
#' the observed states below each branch, given its state, are carried along the branch by the
#' analytical solution of the backward equations of the two-state process, and combined at each
#' node by the rule of speciation (allopatric, or with probability `sympatry` leaving both
#' daughters present). `fit_damocles()` finds the maximum-likelihood rates, and optionally
#' `sympatry`.
#'
#' At the root, `"both"` averages the likelihoods given an absent and a present ancestor, as the
#' default of the DAMOCLES package (`pchoice = 0`); the article sums them, which only adds
#' `log(2)` to the log-likelihood. The optimisation works on the logarithm of the rates, with
#' [stats::optim()] from several starting points, so estimates of zero appear as very small
#' rates. Under this model the tip states hold information on `gamma` and `mu` only through
#' their pattern on the tree, so with few species, or with a community holding nearly all or
#' nearly none of them, the likelihood surface is flat and the estimates are poor. In
#' `data-raw/validation_damocles.R`, ten communities simulated with `gamma = mu = 1` on trees of
#' 60 tips gave estimates of `mu` from zero to several thousand, with the same likelihoods as
#' `DAMOCLES::DAMOCLES_ML()`; the ratio of the two rates is better estimated than each one.
#'
#' With `shared = TRUE`, `fit_damocles()` estimates one pair of rates for all communities by
#' maximising the sum of their log-likelihoods, which treats the communities as independent
#' realisations of the same process on the same tree. Pooling communities narrows the
#' estimates: in `data-raw/validation_damocles.R`, with `gamma = mu = 1` on trees of 60 tips and
#' height one (30 fits each), 90% of the estimates of `gamma` fell between 0.49 and 6.0 from one
#' community and between 0.75 and 1.27 from twenty, and those of `mu` between zero and 18.9 from
#' one and between zero and 1.8 from twenty. `mu` stays the harder rate, and some fits still put
#' it near zero. The
#' assumption of independence matters: neighbouring sites that share their colonisation history,
#' or communities nested in one another, are not independent replicates and pooling them
#' overstates the information on the rates. Twice the difference between the summed
#' log-likelihoods of separate and shared fits can be compared with a chi-squared distribution
#' on \eqn{2 (m - 1)}{2 * (m - 1)} degrees of freedom, with \eqn{m} communities, to ask whether
#' the rates differ among communities.
#'
#' With `sympatry = NULL`, `fit_damocles()` also estimates the probability that a present
#' ancestor leaves both daughters present, on the logit scale. A single community says almost
#' nothing about it: in `data-raw/validation_damocles.R`, with `sympatry = 0.5` and
#' `gamma = mu = 1` on trees of 60 tips, 90% of the estimates fell between 0 and 1 from one
#' community, between 0.26 and 0.81 from five and between 0.41 and 0.60 from twenty. Estimate it
#' with shared rates, and compare the fit with that of `sympatry = 0` by a likelihood ratio (one
#' degree of freedom, on the boundary of the parameter space, so the chi-squared p-value is
#' conservative) to ask whether strictly allopatric speciation is enough.
#'
#' Rates can change across the tree. `damocles_loglik()` and [simulate_damocles()] take one rate
#' per edge, and `fit_damocles()` estimates rates that differ among `regimes` of the edges (a
#' partition defined a priori with [clade_regimes()], or the configuration of rate shifts of a
#' method such as BAMM) or that change with an `edge_covariate` (such as branch-specific rates
#' of speciation from a method such as ClaDS), for colonisation, extinction or both (`vary`).
#' Either input is a plain vector in the order of `tree$edge`, so any method can supply it.
#' Rates of diversification describe speciation and extinction of the whole lineage, not local
#' colonisation and extinction: the regimes or the covariate are a hypothesis about the local
#' dynamics, tested by comparing the fit with the constant model by [damocles_lrt()] (a single
#' regime, or `beta = 0`, gives the constant model) and by intervals from
#' [confint.mc_damocles_fit()]. Each regime adds parameters to a likelihood that is already
#' flat, so these fits need several communities with `shared = TRUE`; the type I error and
#' power of that comparison are in `data-raw/validation_damocles.R`: with rates shared by 5 or
#' 20 communities on trees of 60 tips, the test rejected constant rates in 4--7% of 100
#' replicates without an effect and in 86--99% with a threefold difference between clades or a
#' slope of one on an autocorrelated covariate. Those communities were simulated
#' independently; real communities that share their history will reject more often.
#'
#' @param comm Presences of the species (columns, named as the tips of `tree`) in one or more
#'   communities (rows): a matrix or data frame, or a named vector for one community. Values
#'   above zero count as presences. Tips missing from `comm` count as absent, with a warning,
#'   and species of `comm` missing from the tree are dropped.
#' @inheritParams simulate_damocles
#' @param root Prior on the state of the ancestor at the root: `"both"`, `"absent"` or
#'   `"present"`.
#' @return `damocles_loglik()`: the log-likelihood of each community. `fit_damocles()`: a data
#'   frame of class `mc_damocles_fit` with one row per community and the columns `gamma`, `mu`,
#'   `sympatry` (fixed or estimated), `loglik`, `richness` and `convergence` (0 when [stats::optim()] converged). With
#'   `shared = TRUE` every row holds the same rates and convergence code, `loglik` is the
#'   contribution of each community at those rates, and the attribute `loglik_total` holds
#'   their sum.
#' @references Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for phylogenetic
#'   community structure. *Ecology Letters*, 18, 153--163. \doi{10.1111/ele.12395}
#' @seealso [simulate_damocles()], [damocles_test()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(2)
#' tr <- treeshapesim::sim_tree_shape(60)
#' comm <- simulate_damocles(tr, gamma = 1, mu = 2, n = 2)
#' damocles_loglik(comm, tr, gamma = 1, mu = 2)
#' fit_damocles(comm, tr)
#' fit_damocles(comm, tr, shared = TRUE)
#' @export
damocles_loglik <- function(comm, tree, gamma, mu, root = c("both", "absent", "present"),
                            sympatry = 0) {
  root <- match.arg(root)
  check_damocles_tree(tree)
  rates <- damocles_edge_rates(tree, gamma, mu)
  check_sympatry(sympatry)
  presence <- damocles_presence(comm, tree)
  post <- damocles_postorder(tree)
  damocles_prune(presence, post$tree, rates$gamma[post$index], rates$mu[post$index], root,
                 sympatry, rates$root)
}

#' @rdname damocles_loglik
#' @param start Optional starting rates, `c(gamma = , mu = )`; by default the optimisation
#'   starts from four combinations of 0.2 and 2 divided by the height of the tree.
#' @param shared If `TRUE`, one set of parameters shared by all communities, estimated from the
#'   sum of their log-likelihoods; if `FALSE` (the default), separate parameters for each
#'   community.
#' @param sympatry In `damocles_loglik()`, the probability that a present ancestor leaves both
#'   daughters present (see [simulate_damocles()]). In `fit_damocles()`, a fixed value, or
#'   `NULL` to estimate it with the rates (shared or not, as they are).
#' @param regimes Optional regime of each edge of `tree`, in the order of `tree$edge` (a factor
#'   or character vector), such as from [clade_regimes()] or from a configuration of rate shifts
#'   estimated by another method: the rates in `vary` are estimated separately for each regime.
#' @param edge_covariate Optional positive number for each edge of `tree`, in the order of
#'   `tree$edge`, such as branch-specific rates of speciation estimated by another method: the
#'   rates in `vary` change with it as \eqn{\mathrm{rate} \cdot (x / \bar{x}_g)^\beta}{rate *
#'   (x / geometric mean of x)^beta}, where \eqn{x} is the covariate and \eqn{\bar{x}_g}{the
#'   geometric mean} its geometric mean, with `rate` and \eqn{\beta}{beta} estimated. Give
#'   `regimes` or `edge_covariate`, not both.
#' @param vary Which rates change among `regimes` or with `edge_covariate`: `"both"`, `"gamma"`
#'   or `"mu"`; the other one is a single rate for the whole tree.
#' @export
fit_damocles <- function(comm, tree, root = c("both", "absent", "present"), start = NULL,
                         shared = FALSE, sympatry = 0, regimes = NULL, edge_covariate = NULL,
                         vary = c("both", "gamma", "mu")) {
  root <- match.arg(root)
  vary <- match.arg(vary)
  if (!isTRUE(shared) && !isFALSE(shared)) stop("`shared` must be TRUE or FALSE", call. = FALSE)
  if (!is.null(sympatry)) check_sympatry(sympatry)
  check_damocles_tree(tree)
  presence <- damocles_presence(comm, tree)
  design <- damocles_design(tree, regimes, edge_covariate, vary, sympatry)
  post <- damocles_postorder(tree)
  starts <- damocles_par_starts(damocles_starts(post$tree, start), design)
  n <- nrow(presence)
  groups <- if (shared) list(seq_len(n)) else as.list(seq_len(n))
  fits <- lapply(groups, function(rows) {
    objective <- damocles_objective(presence[rows, , drop = FALSE], tree, post, root, design)
    best <- damocles_optim(objective, starts)
    list(par = best$par, loglik = if (best$value >= 1e10) -Inf else -best$value,
         convergence = best$convergence)
  })
  which_fit <- if (shared) rep(1L, n) else seq_len(n)
  values <- do.call(rbind, lapply(fits, function(f) damocles_par_values(f$par, design)))
  values <- values[which_fit, , drop = FALSE]
  if (!design$estimate_sympatry) values <- cbind(values, sympatry = design$sympatry)
  edge_rates <- lapply(fits, function(f) {
    r <- damocles_par_rates(f$par, design)
    list(gamma = rep_len(r$gamma, nrow(tree$edge)), mu = rep_len(r$mu, nrow(tree$edge)),
         sympatry = r$sympatry)
  })
  loglik <- if (shared) {
    r <- edge_rates[[1]]
    unname(damocles_prune(presence, post$tree, r$gamma[post$index], r$mu[post$index], root,
                          r$sympatry, damocles_root_rates(tree, r$gamma, r$mu)))
  } else {
    vapply(fits, `[[`, numeric(1), "loglik")
  }
  out <- data.frame(values, loglik = loglik, richness = rowSums(presence),
                    convergence = vapply(fits, `[[`, integer(1), "convergence")[which_fit],
                    row.names = rownames(presence), check.names = FALSE)
  if (shared) attr(out, "loglik_total") <- fits[[1]]$loglik
  attr(out, "parameters") <- setdiff(colnames(values), "sympatry")
  attr(out, "n_par") <- design$n_par
  attr(out, "edge_rates") <- edge_rates[which_fit]
  attr(out, "damocles") <- list(presence = presence, tree = tree, root = root, design = design,
                                shared = shared, par = lapply(fits, `[[`, "par"),
                                loglik = vapply(fits, `[[`, numeric(1), "loglik"))
  class(out) <- c("mc_damocles_fit", "data.frame")
  out
}

#' @export
print.mc_damocles_fit <- function(x, ...) {
  cat("DAMOCLES rates by maximum likelihood (Pigot & Etienne 2015)\n")
  total <- attr(x, "loglik_total")
  if (!is.null(total)) {
    cat("Parameters shared by", nrow(attr(x, "damocles")$presence %||% x),
        "communities; total log-likelihood",
        format(total, digits = 6), "\n")
  }
  shown <- as.data.frame(unclass(x), row.names = attr(x, "row.names"))
  print(shown, digits = 3, ...)
  invisible(x)
}

#' Phylogenetic community structure against the dynamic null model of DAMOCLES
#'
#' `damocles_test()` compares the mean pairwise distance (MPD) and the mean nearest taxon
#' distance (MNTD) among the species present in each community with two null distributions: one
#' from communities simulated along the tree under the DAMOCLES model (Pigot & Etienne 2015;
#' [simulate_damocles()]), with the rates of colonisation and local extinction estimated from
#' that community (or shared by all communities) by [fit_damocles()] or given, and one from
#' random draws of the same number of species from the tips, the static null model the dynamic
#' one is compared with in that study.
#'
#' Allopatric speciation alone makes present species less related than random draws expect, so
#' communities that look overdispersed against random draws may be consistent with the dynamic
#' model (Pigot & Etienne 2015). Richness is not fixed in the communities simulated under
#' DAMOCLES; simulated communities with fewer than two species have no MPD or MNTD and are left
#' out of the null distribution, and `n_valid` counts the ones kept. Distances are the
#' cophenetic distances of `tree`, which should hold every species of the regional pool.
#'
#' Standardised effect sizes are \eqn{(x_{obs} - \bar{x}_{null}) / s_{null}}{(observed - mean
#' of the null) / sd of the null}, as in
#' [summary_stats()]: negative values mean the species present are more closely related than
#' expected (clustering) and positive values that they are less related (overdispersion). Note
#' that `DAMOCLES::DAMOCLES_bootstrap()` reports the opposite sign. The p-values are lower-tail
#' ranks, \eqn{(r + 1) / (n + 1)}{(r + 1) / (n + 1)}, with \eqn{r} the number of null values at
#' or below the observed one and \eqn{n} the number of valid null values; values near one point
#' to
#' overdispersion.
#'
#' @inheritParams damocles_loglik
#' @param n_sim Number of null communities of each model.
#' @param gamma,mu Rates of local colonisation and extinction used to simulate the null
#'   communities, one value for all communities or one per community, or lists with one rate,
#'   or one rate per edge of `tree`, per community; `NULL` (the default) estimates both with
#'   [fit_damocles()], with `shared`, `sympatry`, `regimes`, `edge_covariate` and `vary` passed
#'   on to it.
#' @param shared When the rates are estimated, `TRUE` fits one pair shared by all communities
#'   (see [fit_damocles()]) and `FALSE` (the default) one pair per community.
#' @param sympatry Probability that a present ancestor leaves both daughters present, used in
#'   the simulations (see [simulate_damocles()]); `NULL` estimates it with the rates.
#' @param root Prior on the state at the root, used both in the fit and in the simulations
#'   (`"both"` draws the state of the root of each simulated community with equal
#'   probabilities).
#' @return A data frame of class `mc_damocles_test`, one row per community, with the richness,
#'   the parameters used (`gamma`, `mu`, or the columns of the fit with `regimes` or
#'   `edge_covariate`, and `sympatry`; `NA` for rates given per edge) and their log-likelihood when estimated,
#'   and, for `mpd` and `mntd`, the observed value, the standardised effect size and p-value against the DAMOCLES
#'   null (`ses_mpd_damocles`, `p_mpd_damocles`, ...) and against random draws
#'   (`ses_mpd_random`, `p_mpd_random`, ...), plus `n_valid`.
#' @references Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for phylogenetic
#'   community structure. *Ecology Letters*, 18, 153--163. \doi{10.1111/ele.12395}
#' @seealso [simulate_damocles()], [fit_damocles()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(3)
#' tr <- treeshapesim::sim_tree_shape(40)
#' comm <- simulate_damocles(tr, gamma = 1, mu = 1, n = 2)
#' damocles_test(comm, tr, n_sim = 99)
#' @export
damocles_test <- function(comm, tree, n_sim = 999, gamma = NULL, mu = NULL,
                          root = c("both", "absent", "present"), shared = FALSE,
                          sympatry = 0, regimes = NULL, edge_covariate = NULL,
                          vary = c("both", "gamma", "mu")) {
  root <- match.arg(root)
  vary <- match.arg(vary)
  check_damocles_tree(tree)
  if (length(n_sim) != 1 || !is.finite(n_sim) || n_sim < 2 || n_sim != round(n_sim)) {
    stop("`n_sim` must be one whole number, 2 or more", call. = FALSE)
  }
  presence <- damocles_presence(comm, tree)
  n_comm <- nrow(presence)
  if (xor(is.null(gamma), is.null(mu))) {
    stop("give both `gamma` and `mu`, or neither to estimate them", call. = FALSE)
  }
  if (is.null(gamma)) {
    fit <- fit_damocles(presence, tree, root = root, shared = shared, sympatry = sympatry,
                        regimes = regimes, edge_covariate = edge_covariate, vary = vary)
    parameters <- as.data.frame(unclass(fit))[, attr(fit, "parameters"), drop = FALSE]
    edge_rates <- attr(fit, "edge_rates")
    sympatry <- fit$sympatry
    loglik <- fit$loglik
  } else {
    if (!is.null(regimes) || !is.null(edge_covariate)) {
      stop("`regimes` and `edge_covariate` are for estimated rates; give rates per edge ",
           "as lists in `gamma` and `mu` instead", call. = FALSE)
    }
    if (!length(gamma) %in% c(1, n_comm) || !length(mu) %in% c(1, n_comm)) {
      stop("`gamma` and `mu` need one value or one per community", call. = FALSE)
    }
    gamma <- rep_len(as.list(gamma), n_comm)
    mu <- rep_len(as.list(mu), n_comm)
    edge_rates <- lapply(seq_len(n_comm), function(i) damocles_edge_rates(tree, gamma[[i]],
                                                                          mu[[i]]))
    one <- function(rate) vapply(rate, function(r) if (length(r) == 1) r else NA_real_,
                                 numeric(1))
    parameters <- data.frame(gamma = one(gamma), mu = one(mu))
    if (is.null(sympatry)) {
      stop("`sympatry = NULL` (estimated) needs `gamma` and `mu` estimated too", call. = FALSE)
    }
    check_sympatry(sympatry)
    sympatry <- rep_len(sympatry, n_comm)
    loglik <- NA_real_
  }
  distances <- tip_distances(tree)[colnames(presence), colnames(presence)]
  rows <- lapply(seq_len(n_comm), function(i) {
    site <- presence[i, , drop = FALSE]
    observed <- c(mpd = unname(mean_pairwise(site, distances)),
                  mntd = unname(mean_nearest(site, distances)))
    dynamic <- simulate_damocles(tree, edge_rates[[i]]$gamma, edge_rates[[i]]$mu, n = n_sim,
                                 root = root, sympatry = sympatry[i])
    dynamic <- dynamic[, colnames(presence), drop = FALSE]
    random <- t(replicate(n_sim, sample(unname(site[1, ]))))
    null_dynamic <- cbind(mpd = mean_pairwise(dynamic, distances),
                          mntd = mean_nearest(dynamic, distances))
    null_random <- cbind(mpd = mean_pairwise(random, distances),
                         mntd = mean_nearest(random, distances))
    out <- c(richness = sum(site))
    for (metric in c("mpd", "mntd")) {
      d <- null_effect(observed[[metric]], null_dynamic[, metric])
      r <- null_effect(observed[[metric]], null_random[, metric])
      out[paste0(c("", "ses_", "p_", "ses_", "p_"), metric,
                 c("", "_damocles", "_damocles", "_random", "_random"))] <-
        c(observed[[metric]], d[["ses"]], d[["p"]], r[["ses"]], r[["p"]])
    }
    out[["n_valid"]] <- sum(stats::complete.cases(null_dynamic))
    out
  })
  stats <- as.data.frame(do.call(rbind, rows))
  out <- data.frame(richness = stats$richness, parameters, sympatry = sympatry, loglik = loglik,
                    stats[, -1], row.names = rownames(presence), check.names = FALSE)
  attr(out, "n_sim") <- n_sim
  attr(out, "parameters") <- colnames(parameters)
  class(out) <- c("mc_damocles_test", "data.frame")
  out
}

#' @export
print.mc_damocles_test <- function(x, ...) {
  cat("Phylogenetic structure against the DAMOCLES null model (Pigot & Etienne 2015)\n")
  cat(attr(x, "n_sim"), "null communities per model; SES < 0: clustering, SES > 0:",
      "overdispersion\n")
  shown <- as.data.frame(unclass(x), row.names = attr(x, "row.names"))
  columns <- c("richness", attr(x, "parameters"), if (any(x$sympatry > 0)) "sympatry", "mpd",
               "ses_mpd_damocles", "p_mpd_damocles", "ses_mpd_random", "p_mpd_random", "mntd", "ses_mntd_damocles",
               "p_mntd_damocles", "ses_mntd_random", "p_mntd_random")
  print(shown[, columns], digits = 3, ...)
  invisible(x)
}

# --- internals ---------------------------------------------------------------------------

check_damocles_tree <- function(tree) {
  if (!inherits(tree, "phylo")) stop("`tree` must be a `phylo` object", call. = FALSE)
  if (is.null(tree$edge.length) || anyNA(tree$edge.length) || any(tree$edge.length < 0)) {
    stop("`tree` needs branch lengths, zero or more", call. = FALSE)
  }
  if (!ape::is.rooted(tree) || !ape::is.binary(tree)) {
    stop("`tree` must be rooted and fully dichotomous; see ape::multi2di()", call. = FALSE)
  }
  if (is.null(tree$tip.label) || anyDuplicated(tree$tip.label)) {
    stop("`tree` needs different tip labels", call. = FALSE)
  }
  invisible(tree)
}

# rates of each edge, in the order of `tree$edge`, and of the root edge, which takes the rates of
# the two edges leaving the root
damocles_edge_rates <- function(tree, gamma, mu) {
  n_edge <- nrow(tree$edge)
  for (rate in list(gamma, mu)) {
    if (!is.numeric(rate) || !length(rate) %in% c(1, n_edge) || any(!is.finite(rate)) ||
        any(rate < 0)) {
      stop("`gamma` and `mu` must each be one rate, or one per edge of `tree` in the order of ",
           "`tree$edge`, zero or more", call. = FALSE)
    }
  }
  gamma <- rep_len(as.numeric(gamma), n_edge)
  mu <- rep_len(as.numeric(mu), n_edge)
  list(gamma = gamma, mu = mu, root = damocles_root_rates(tree, gamma, mu))
}

damocles_root_rates <- function(tree, gamma, mu) {
  at_root <- tree$edge[, 1] == ape::Ntip(tree) + 1
  rates <- c(gamma = gamma[at_root][1], mu = mu[at_root][1])
  if (!is.null(tree$root.edge) &&
      (length(unique(gamma[at_root])) > 1 || length(unique(mu[at_root])) > 1)) {
    stop("the root edge takes the rates of the two edges leaving the root, which differ; ",
         "give them the same rates (or regime) or drop `tree$root.edge`", call. = FALSE)
  }
  rates
}

# `tree` in postorder and the positions of its edges in the original `tree$edge`
damocles_postorder <- function(tree) {
  list(tree = ape::reorder.phylo(tree, "postorder"),
       index = ape::reorder.phylo(tree, "postorder", index.only = TRUE))
}

check_sympatry <- function(sympatry) {
  if (length(sympatry) != 1 || !is.numeric(sympatry) || !is.finite(sympatry) ||
      sympatry < 0 || sympatry > 1) {
    stop("`sympatry` must be one probability, between 0 and 1", call. = FALSE)
  }
  invisible(TRUE)
}

# presences (0/1) of the tips of `tree`, in their order, one row per community
damocles_presence <- function(comm, tree) {
  if (is.null(dim(comm))) {
    if (is.null(names(comm))) stop("`comm` needs species names", call. = FALSE)
    comm <- matrix(comm, 1, dimnames = list(NULL, names(comm)))
  }
  comm <- check_comm(comm)
  if (is.null(colnames(comm))) stop("`comm` needs species names in its columns", call. = FALSE)
  if (anyDuplicated(colnames(comm))) {
    stop("repeated species: ", species_list(unique(colnames(comm)[duplicated(colnames(comm))])),
         call. = FALSE)
  }
  tips <- tree$tip.label
  only_comm <- setdiff(colnames(comm), tips)
  only_tree <- setdiff(tips, colnames(comm))
  if (length(only_comm) > 0) {
    warning(length(only_comm), " species of `comm` not in the phylogeny dropped: ",
            species_list(only_comm), call. = FALSE)
  }
  if (length(only_tree) > 0) {
    warning(length(only_tree), " tips of the phylogeny not in `comm` counted as absent: ",
            species_list(only_tree), call. = FALSE)
  }
  out <- matrix(0, nrow(comm), length(tips),
                dimnames = list(rownames(comm) %||% paste0("site", seq_len(nrow(comm))), tips))
  shared <- intersect(tips, colnames(comm))
  out[, shared] <- (comm[, shared, drop = FALSE] > 0) * 1
  out
}

# probability of being present at the end of a branch of length `t`, given the state at its
# start: the transition probabilities of the two-state Markov process
damocles_transition <- function(t, gamma, mu) {
  total <- gamma + mu
  if (total == 0) return(c(from_absent = 0, from_present = 1))
  decay <- exp(-total * t)
  c(from_absent = gamma * (1 - decay) / total, from_present = (gamma + mu * decay) / total)
}

damocles_step <- function(begin, t, gamma, mu) {
  p <- damocles_transition(t, gamma, mu)
  as.integer(stats::runif(length(begin)) < ifelse(begin == 1, p[[2]], p[[1]]))
}

# log-likelihood of each row of `presence` (tips of `tree`, in order), `tree` in postorder
# `gamma` and `mu`: one rate, or one per edge in the order of `tree$edge`; `root_rates`: the
# rates of the root edge
damocles_prune <- function(presence, tree, gamma, mu, root, sympatry = 0,
                           root_rates = c(gamma[1], mu[1])) {
  n_edge <- nrow(tree$edge)
  gamma <- rep_len(gamma, n_edge)
  mu <- rep_len(mu, n_edge)
  n_tip <- ape::Ntip(tree)
  n_comm <- nrow(presence)
  # partial likelihoods given absence (column 1) and presence (column 2) at each node
  absent <- matrix(1, n_comm, n_tip + tree$Nnode)
  present <- absent
  absent[, seq_len(n_tip)] <- 1 - presence
  present[, seq_len(n_tip)] <- presence
  # the first daughter carried up to the node waits here for the second one
  first_absent <- matrix(NA_real_, n_comm, n_tip + tree$Nnode)
  first_present <- first_absent
  loglik <- numeric(n_comm)
  carry <- function(a, p, t, gamma, mu) {
    to_present <- damocles_transition(t, gamma, mu)
    list(absent = (1 - to_present[[1]]) * a + to_present[[1]] * p,
         present = (1 - to_present[[2]]) * a + to_present[[2]] * p)
  }
  for (e in seq_len(nrow(tree$edge))) {  # postorder: descendants before ancestors
    parent <- tree$edge[e, 1]
    child <- tree$edge[e, 2]
    top <- carry(absent[, child], present[, child], tree$edge.length[e], gamma[e], mu[e])
    if (all(is.na(first_absent[, parent]))) {
      first_absent[, parent] <- top$absent
      first_present[, parent] <- top$present
      next
    }
    # an absent ancestor leaves two absent daughters; a present one leaves both present with
    # probability `sympatry`, and otherwise one present daughter, either of the two
    a <- first_absent[, parent] * top$absent
    p <- (1 - sympatry) *
      (first_present[, parent] * top$absent + first_absent[, parent] * top$present) / 2 +
      sympatry * first_present[, parent] * top$present
    scale <- a + p
    loglik <- loglik + log(scale)
    absent[, parent] <- a / ifelse(scale > 0, scale, 1)
    present[, parent] <- p / ifelse(scale > 0, scale, 1)
  }
  crown <- n_tip + 1
  a <- absent[, crown]
  p <- present[, crown]
  if (!is.null(tree$root.edge)) {
    top <- carry(a, p, tree$root.edge, root_rates[[1]], root_rates[[2]])
    a <- top$absent
    p <- top$present
  }
  at_root <- switch(root, both = (a + p) / 2, absent = a, present = p)
  stats::setNames(loglik + log(at_root), rownames(presence))
}

damocles_starts <- function(tree, start) {
  if (!is.null(start)) {
    if (length(start) != 2 || any(!is.finite(start)) || any(start <= 0)) {
      stop("`start` needs two positive rates, gamma and mu", call. = FALSE)
    }
    return(matrix(start, 1))
  }
  height <- max(ape::node.depth.edgelength(tree))
  as.matrix(expand.grid(gamma = c(0.2, 2), mu = c(0.2, 2))) / height
}

# the parameters of a fit: rates on the log scale (one per regime for the rates in `vary`),
# slopes on the covariate, and sympatry on the logit scale when estimated
damocles_design <- function(tree, regimes, edge_covariate, vary, sympatry) {
  n_edge <- nrow(tree$edge)
  at_root <- tree$edge[, 1] == ape::Ntip(tree) + 1
  root_msg <- paste("the root edge takes the regime of the two edges leaving the root, which",
                    "differ; give them the same one or drop `tree$root.edge`")
  if (!is.null(regimes) && !is.null(edge_covariate)) {
    stop("give `regimes` or `edge_covariate`, not both", call. = FALSE)
  }
  type <- "constant"
  z <- NULL
  if (!is.null(regimes)) {
    if (length(regimes) != n_edge || anyNA(regimes)) {
      stop("`regimes` needs one regime for each edge of `tree`, in the order of `tree$edge`, ",
           "with no missing values", call. = FALSE)
    }
    regimes <- factor(regimes)
    if (!is.null(tree$root.edge) && length(unique(regimes[at_root])) > 1) {
      stop(root_msg, call. = FALSE)
    }
    type <- "regimes"
  }
  if (!is.null(edge_covariate)) {
    if (!is.numeric(edge_covariate) || length(edge_covariate) != n_edge ||
        any(!is.finite(edge_covariate)) || any(edge_covariate <= 0)) {
      stop("`edge_covariate` needs one positive number for each edge of `tree`, in the order ",
           "of `tree$edge`", call. = FALSE)
    }
    z <- log(edge_covariate) - mean(log(edge_covariate))
    if (!is.null(tree$root.edge) && length(unique(z[at_root])) > 1) {
      stop(sub("regime", "covariate", root_msg), call. = FALSE)
    }
    type <- "covariate"
  }
  by_regime <- function(rate) type == "regimes" && vary %in% c("both", rate) &&
    nlevels(regimes) > 1
  n_gamma <- if (by_regime("gamma")) nlevels(regimes) else 1L
  n_mu <- if (by_regime("mu")) nlevels(regimes) else 1L
  beta_gamma <- type == "covariate" && vary %in% c("both", "gamma")
  beta_mu <- type == "covariate" && vary %in% c("both", "mu")
  names <- c(if (n_gamma > 1) paste0("gamma_", levels(regimes)) else "gamma",
             if (n_mu > 1) paste0("mu_", levels(regimes)) else "mu",
             if (beta_gamma) "beta_gamma", if (beta_mu) "beta_mu",
             if (is.null(sympatry)) "sympatry")
  list(type = type, regimes = regimes, z = z, n_gamma = n_gamma, n_mu = n_mu,
       beta_gamma = beta_gamma, beta_mu = beta_mu, estimate_sympatry = is.null(sympatry),
       sympatry = sympatry, names = names, n_par = length(names))
}

# rates of each edge (or one rate for all) and sympatry from the working parameters
damocles_par_rates <- function(par, design) {
  used <- 0
  take <- function(n) {
    value <- par[used + seq_len(n)]
    used <<- used + n
    value
  }
  log_gamma <- take(design$n_gamma)
  log_mu <- take(design$n_mu)
  if (design$type == "regimes") {
    r <- as.integer(design$regimes)
    log_gamma <- if (design$n_gamma > 1) log_gamma[r] else rep(log_gamma, length(r))
    log_mu <- if (design$n_mu > 1) log_mu[r] else rep(log_mu, length(r))
  } else if (design$type == "covariate") {
    if (design$beta_gamma) log_gamma <- log_gamma + take(1) * design$z
    if (design$beta_mu) log_mu <- log_mu + take(1) * design$z
  }
  list(gamma = exp(log_gamma), mu = exp(log_mu),
       sympatry = if (design$estimate_sympatry) stats::plogis(take(1)) else design$sympatry)
}

# working parameter `j` on its natural scale
damocles_natural <- function(value, j, design) {
  if (j <= design$n_gamma + design$n_mu) return(exp(value))
  if (design$estimate_sympatry && j == design$n_par) return(stats::plogis(value))
  value
}

damocles_par_values <- function(par, design) {
  stats::setNames(vapply(seq_along(par), function(j) damocles_natural(par[j], j, design),
                         numeric(1)), design$names)
}

damocles_par_starts <- function(base, design) {
  out <- t(apply(base, 1, function(b) {
    c(rep(log(b[1]), design$n_gamma), rep(log(b[2]), design$n_mu),
      rep(0, design$beta_gamma + design$beta_mu))
  }))
  if (design$estimate_sympatry) {
    out <- cbind(out[rep(seq_len(nrow(out)), 2), , drop = FALSE],
                 rep(stats::qlogis(c(0.1, 0.5)), each = nrow(out)))
  }
  out
}

# minus the summed log-likelihood of `presence` as a function of the working parameters
damocles_objective <- function(presence, tree, post, root, design) {
  n_edge <- nrow(tree$edge)
  function(par) {
    r <- damocles_par_rates(par, design)
    gamma <- rep_len(r$gamma, n_edge)
    mu <- rep_len(r$mu, n_edge)
    value <- sum(damocles_prune(presence, post$tree, gamma[post$index], mu[post$index], root,
                                r$sympatry, damocles_root_rates(tree, gamma, mu)))
    if (is.finite(value)) -value else 1e10
  }
}

damocles_optim <- function(objective, starts) {
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    fit <- damocles_minimise(objective, starts[i, ])
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  best
}

damocles_minimise <- function(objective, start) {
  if (length(start) == 1) {
    fit <- stats::optimize(objective, start + c(-15, 15))
    return(list(par = fit$minimum, value = fit$objective, convergence = 0L))
  }
  control <- list(reltol = 1e-10, maxit = 4000)
  fit <- stats::optim(start, objective, control = control)
  # Nelder-Mead stalls more often in more dimensions: restart from where it stopped
  if (length(start) > 2) fit <- stats::optim(fit$par, objective, control = control)
  list(par = fit$par, value = fit$value, convergence = as.integer(fit$convergence))
}

# standardised effect size and lower-tail p-value of `observed` against the finite `null`
null_effect <- function(observed, null) {
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) < 2) return(c(ses = NA_real_, p = NA_real_))
  spread <- stats::sd(null)
  c(ses = if (spread > 0) (observed - mean(null)) / spread else NA_real_,
    p = (sum(null <= observed) + 1) / (length(null) + 1))
}
