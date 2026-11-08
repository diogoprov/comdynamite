#' Local presence along a tree: the dynamic null model of DAMOCLES
#'
#' `simulate_damocles()` simulates the presence and absence of each lineage of a clade in one
#' local community, as it diversifies along `tree`, under the dynamic null model of community
#' assembly of Pigot & Etienne (2015), DAMOCLES. Along each branch the lineage is either locally
#' absent (state 0) or present (state 1); an absent lineage colonises the community at rate
#' `gamma` and a present one goes locally extinct at rate `mu`, the same for every lineage. At
#' each speciation, which is allopatric, a locally present ancestor leaves one daughter present
#' and the other absent, chosen at random, and an absent ancestor leaves two absent daughters.
#' The states at the tips form the community.
#'
#' The states at the nodes are drawn from the exact transition probabilities of the two-state
#' process over each branch, which gives the same distribution of tip states as simulating each
#' colonisation and extinction event. The model assumes the tree holds the species of the
#' regional pool, so it is meant for ultrametric trees of extant species. Branches of length
#' zero and a root edge (`tree$root.edge`) are allowed.
#'
#' @param tree A rooted, fully dichotomous `phylo` tree with branch lengths.
#' @param gamma Rate of local colonisation per lineage, in the units of the branch lengths.
#' @param mu Rate of local extinction per lineage.
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
simulate_damocles <- function(tree, gamma, mu, n = 1, root = c("both", "absent", "present")) {
  root <- match.arg(root)
  check_damocles_tree(tree)
  check_damocles_rates(gamma, mu)
  if (length(n) != 1 || !is.finite(n) || n < 1 || n != round(n)) {
    stop("`n` must be one whole number, 1 or more", call. = FALSE)
  }
  tree <- ape::reorder.phylo(tree, "cladewise")
  n_tip <- ape::Ntip(tree)
  start <- switch(root, both = stats::rbinom(n, 1, 0.5), absent = rep(0L, n),
                  present = rep(1L, n))
  if (!is.null(tree$root.edge)) start <- damocles_step(start, tree$root.edge, gamma, mu)
  # state at the end of the branch leading to each node, just before it speciates
  state <- matrix(0L, n, n_tip + tree$Nnode)
  state[, n_tip + 1] <- start
  # which of the two daughters inherits the presence of its ancestor
  heir <- matrix(sample.int(2, n * tree$Nnode, replace = TRUE), n, tree$Nnode)
  position <- stats::ave(seq_len(nrow(tree$edge)), tree$edge[, 1], FUN = seq_along)
  for (e in seq_len(nrow(tree$edge))) {  # cladewise: ancestors before descendants
    parent <- tree$edge[e, 1]
    begin <- state[, parent] * (heir[, parent - n_tip] == position[e])
    state[, tree$edge[e, 2]] <- damocles_step(begin, tree$edge.length[e], gamma, mu)
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
#' node by the rule of allopatric speciation. `fit_damocles()` finds the maximum-likelihood
#' rates.
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
#' on `2 * (number of communities - 1)` degrees of freedom to ask whether the rates differ among
#' communities.
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
#'   `loglik`, `richness` and `convergence` (0 when [stats::optim()] converged). With
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
damocles_loglik <- function(comm, tree, gamma, mu, root = c("both", "absent", "present")) {
  root <- match.arg(root)
  check_damocles_tree(tree)
  check_damocles_rates(gamma, mu)
  presence <- damocles_presence(comm, tree)
  damocles_prune(presence, ape::reorder.phylo(tree, "postorder"), gamma, mu, root)
}

#' @rdname damocles_loglik
#' @param start Optional starting rates, `c(gamma = , mu = )`; by default the optimisation
#'   starts from four combinations of 0.2 and 2 divided by the height of the tree.
#' @param shared If `TRUE`, one pair of rates shared by all communities, estimated from the sum
#'   of their log-likelihoods; if `FALSE` (the default), separate rates for each community.
#' @export
fit_damocles <- function(comm, tree, root = c("both", "absent", "present"), start = NULL,
                         shared = FALSE) {
  root <- match.arg(root)
  if (!isTRUE(shared) && !isFALSE(shared)) stop("`shared` must be TRUE or FALSE", call. = FALSE)
  check_damocles_tree(tree)
  presence <- damocles_presence(comm, tree)
  tree <- ape::reorder.phylo(tree, "postorder")
  starts <- damocles_starts(tree, start)
  groups <- if (shared) list(seq_len(nrow(presence))) else as.list(seq_len(nrow(presence)))
  fits <- lapply(groups, function(rows) {
    damocles_optim(presence[rows, , drop = FALSE], tree, root, starts)
  })
  if (shared) {
    fit <- fits[[1]]
    each <- damocles_prune(presence, tree, fit$gamma, fit$mu, root)
    n <- nrow(presence)
    out <- data.frame(gamma = rep(fit$gamma, n), mu = rep(fit$mu, n), loglik = unname(each),
                      richness = rowSums(presence), convergence = rep(fit$convergence, n),
                      row.names = rownames(presence))
    attr(out, "loglik_total") <- fit$loglik
  } else {
    out <- data.frame(gamma = vapply(fits, `[[`, numeric(1), "gamma"),
                      mu = vapply(fits, `[[`, numeric(1), "mu"),
                      loglik = vapply(fits, `[[`, numeric(1), "loglik"),
                      richness = rowSums(presence),
                      convergence = vapply(fits, `[[`, integer(1), "convergence"),
                      row.names = rownames(presence))
  }
  class(out) <- c("mc_damocles_fit", "data.frame")
  out
}

#' @export
print.mc_damocles_fit <- function(x, ...) {
  cat("DAMOCLES rates by maximum likelihood (Pigot & Etienne 2015)\n")
  total <- attr(x, "loglik_total")
  if (!is.null(total)) {
    cat("Rates shared by", nrow(x), "communities; total log-likelihood",
        format(total, digits = 6), "\n")
  }
  print(as.data.frame(unclass(x), row.names = attr(x, "row.names")), digits = 3, ...)
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
#' Standardised effect sizes are `(observed - mean of the null) / sd of the null`, as in
#' [summary_stats()]: negative values mean the species present are more closely related than
#' expected (clustering) and positive values that they are less related (overdispersion). Note
#' that `DAMOCLES::DAMOCLES_bootstrap()` reports the opposite sign. The p-values are lower-tail
#' ranks, `(number of null values <= observed + 1) / (n_valid + 1)`; values near one point to
#' overdispersion.
#'
#' @inheritParams damocles_loglik
#' @param n_sim Number of null communities of each model.
#' @param gamma,mu Rates of local colonisation and extinction used to simulate the null
#'   communities, one value for all communities or one per community; `NULL` (the default)
#'   estimates both with [fit_damocles()].
#' @param shared When the rates are estimated, `TRUE` fits one pair shared by all communities
#'   (see [fit_damocles()]) and `FALSE` (the default) one pair per community.
#' @param root Prior on the state at the root, used both in the fit and in the simulations
#'   (`"both"` draws the state of the root of each simulated community with equal
#'   probabilities).
#' @return A data frame of class `mc_damocles_test`, one row per community, with the richness,
#'   the rates used (`gamma`, `mu`) and their log-likelihood when estimated, and, for `mpd` and
#'   `mntd`, the observed value, the standardised effect size and p-value against the DAMOCLES
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
                          root = c("both", "absent", "present"), shared = FALSE) {
  root <- match.arg(root)
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
    fit <- fit_damocles(presence, tree, root = root, shared = shared)
    gamma <- fit$gamma
    mu <- fit$mu
    loglik <- fit$loglik
  } else {
    if (!length(gamma) %in% c(1, n_comm) || !length(mu) %in% c(1, n_comm)) {
      stop("`gamma` and `mu` need one value or one per community", call. = FALSE)
    }
    gamma <- rep_len(gamma, n_comm)
    mu <- rep_len(mu, n_comm)
    for (i in seq_len(n_comm)) check_damocles_rates(gamma[i], mu[i])
    loglik <- NA_real_
  }
  distances <- tip_distances(tree)[colnames(presence), colnames(presence)]
  rows <- lapply(seq_len(n_comm), function(i) {
    site <- presence[i, , drop = FALSE]
    observed <- c(mpd = unname(mean_pairwise(site, distances)),
                  mntd = unname(mean_nearest(site, distances)))
    dynamic <- simulate_damocles(tree, gamma[i], mu[i], n = n_sim, root = root)
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
  out <- data.frame(richness = stats$richness, gamma = gamma, mu = mu, loglik = loglik,
                    stats[, -1], row.names = rownames(presence))
  attr(out, "n_sim") <- n_sim
  class(out) <- c("mc_damocles_test", "data.frame")
  out
}

#' @export
print.mc_damocles_test <- function(x, ...) {
  cat("Phylogenetic structure against the DAMOCLES null model (Pigot & Etienne 2015)\n")
  cat(attr(x, "n_sim"), "null communities per model; SES < 0: clustering, SES > 0:",
      "overdispersion\n")
  shown <- as.data.frame(unclass(x), row.names = attr(x, "row.names"))
  print(shown[, c("richness", "gamma", "mu", "mpd", "ses_mpd_damocles", "p_mpd_damocles",
                  "ses_mpd_random", "p_mpd_random", "mntd", "ses_mntd_damocles",
                  "p_mntd_damocles", "ses_mntd_random", "p_mntd_random")], digits = 3, ...)
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

check_damocles_rates <- function(gamma, mu) {
  for (rate in list(gamma = gamma, mu = mu)) {
    if (length(rate) != 1 || !is.numeric(rate) || !is.finite(rate) || rate < 0) {
      stop("`gamma` and `mu` must each be one number, zero or more", call. = FALSE)
    }
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
damocles_prune <- function(presence, tree, gamma, mu, root) {
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
  carry <- function(a, p, t) {
    to_present <- damocles_transition(t, gamma, mu)
    list(absent = (1 - to_present[[1]]) * a + to_present[[1]] * p,
         present = (1 - to_present[[2]]) * a + to_present[[2]] * p)
  }
  for (e in seq_len(nrow(tree$edge))) {  # postorder: descendants before ancestors
    parent <- tree$edge[e, 1]
    child <- tree$edge[e, 2]
    top <- carry(absent[, child], present[, child], tree$edge.length[e])
    if (all(is.na(first_absent[, parent]))) {
      first_absent[, parent] <- top$absent
      first_present[, parent] <- top$present
      next
    }
    # allopatric speciation: an absent ancestor leaves two absent daughters, a present one
    # leaves one present daughter, either of the two
    a <- first_absent[, parent] * top$absent
    p <- (first_present[, parent] * top$absent + first_absent[, parent] * top$present) / 2
    scale <- a + p
    loglik <- loglik + log(scale)
    absent[, parent] <- a / ifelse(scale > 0, scale, 1)
    present[, parent] <- p / ifelse(scale > 0, scale, 1)
  }
  crown <- n_tip + 1
  a <- absent[, crown]
  p <- present[, crown]
  if (!is.null(tree$root.edge)) {
    top <- carry(a, p, tree$root.edge)
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

damocles_optim <- function(presence, tree, root, starts) {
  objective <- function(log_rates) {
    value <- sum(damocles_prune(presence, tree, exp(log_rates[1]), exp(log_rates[2]), root))
    if (is.finite(value)) -value else 1e10
  }
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    fit <- stats::optim(log(starts[i, ]), objective, control = list(reltol = 1e-10,
                                                                     maxit = 2000))
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  list(gamma = exp(best$par[1]), mu = exp(best$par[2]),
       loglik = if (best$value >= 1e10) -Inf else -best$value,
       convergence = as.integer(best$convergence))
}

# standardised effect size and lower-tail p-value of `observed` against the finite `null`
null_effect <- function(observed, null) {
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) < 2) return(c(ses = NA_real_, p = NA_real_))
  spread <- stats::sd(null)
  c(ses = if (spread > 0) (observed - mean(null)) / spread else NA_real_,
    p = (sum(null <= observed) + 1) / (length(null) + 1))
}
