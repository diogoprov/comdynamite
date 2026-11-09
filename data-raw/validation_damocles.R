# simulate_damocles(), damocles_loglik() and fit_damocles() against the DAMOCLES package of
# Pigot & Etienne (2015) (GPL, only called here; install.packages("DAMOCLES"), version 2.3).
# 1. Log-likelihoods against DAMOCLES::DAMOCLES_loglik() on 20 trees, three pairs of rates and
#    the three priors at the root (pchoice 0, 1, 2), plus trees with a root edge.
# 2. Maximum-likelihood rates against DAMOCLES::DAMOCLES_ML() on 10 simulated communities.
# 3. Simulations against DAMOCLES::DAMOCLES_sim() (a Gillespie simulation of every event): the
#    frequency of presence of each tip, richness, MPD and MNTD over 4000 communities each.
# 4. Identifiability of the rates shared by several communities (fit_damocles(shared = TRUE)):
#    estimates from 1, 5 and 20 communities simulated with gamma = mu = 1 on trees of 60 tips
#    (height one), 30 replicates each. DAMOCLES fits one community at a time, so this part
#    checks only our code.
# 5. Speciation with probability `sympatry` of both daughters present (not in DAMOCLES): rates
#    and sympatry estimated together, shared, from 1, 5 and 20 communities simulated with
#    gamma = mu = 1 and sympatry = 0.5 on trees of 60 tips, 30 replicates each.
# 6. Rates that change among regimes or with a covariate of the edges: type I error and power
#    of damocles_lrt() against constant rates, 100 replicates per setting, rates shared by 5
#    or 20 communities on trees of 60 tips. The covariate is autocorrelated along the tree,
#    like branch-specific rates of speciation (its log follows Brownian motion with variance
#    0.5 per unit of time); colonisation is gamma = exp(beta * centred log covariate), with
#    beta = 0 (type I error) or 1 (power), and mu = 1. The regimes are the crown clade of one
#    of the two daughters of the root against the rest, with gamma 1 in both (type I error)
#    or 3 in the clade and 0.5 elsewhere (power).
# 7. damocles_partition_test() against the chi-squared test of damocles_lrt(), 40 replicates
#    of 10 communities on trees of 60 tips, 49 random partitions each. The partition tested
#    is the clade (with its stem) whose subtree holds about a quarter of the edges; colonisation
#    is 1 everywhere (no effect), 3 in that clade and 0.5 elsewhere (effect in the clade), or 3
#    in another clade of similar size placed at random and 0.5 elsewhere (effect elsewhere).
#    Random partitions overlap the tested clade at most 0.25 (Jaccard), the default; a first
#    run with max_overlap = 1 rejected in 0 of 40 replicates in every scenario.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(DAMOCLES)
})
cores <- as.integer(Sys.getenv("CORES", "6"))
RNGkind("L'Ecuyer-CMRG")
set.seed(2015)

pa_table <- function(tree, presence) cbind(tree$tip.label, presence[tree$tip.label])

# 1. likelihood ------------------------------------------------------------------------------
rates <- list(c(gamma = 0.8, mu = 1.2), c(gamma = 0.1, mu = 0.05), c(gamma = 3, mu = 0.5))
loglik_diff <- function(i) {
  tree <- treeshapesim::sim_tree_shape(sample(20:80, 1))
  if (i %% 4 == 0) tree$root.edge <- stats::runif(1, 0, 2)
  out <- numeric(0)
  for (r in rates) {
    presence <- simulate_damocles(tree, r[["gamma"]], r[["mu"]])[1, ]
    for (pchoice in 0:2) {
      theirs <- DAMOCLES_loglik(tree, pa_table(tree, presence), c(r[["mu"]], r[["gamma"]], 0),
                                pchoice = pchoice)
      ours <- damocles_loglik(presence, tree, r[["gamma"]], r[["mu"]],
                              root = c("both", "absent", "present")[pchoice + 1])
      out <- c(out, abs(theirs - ours))
    }
  }
  out
}
diffs <- unlist(parallel::mclapply(seq_len(20), loglik_diff, mc.cores = cores))
cat("1. |loglik difference|: max", format(max(diffs), digits = 3), "over", length(diffs),
    "comparisons\n")

# 2. maximum likelihood ----------------------------------------------------------------------
ml_one <- function(i) {
  tree <- treeshapesim::sim_tree_shape(60)
  presence <- simulate_damocles(tree, 1, 1)[1, ]
  ours <- fit_damocles(presence, tree)
  theirs <- suppressMessages(utils::capture.output(
    fit <- DAMOCLES_ML(tree, pa_table(tree, presence), initparsopt = c(0.5, 0.5))
  ))
  c(richness = sum(presence), gamma = ours$gamma, gamma_damocles = fit$gamma_0,
    mu = ours$mu, mu_damocles = fit$mu, loglik = ours$loglik, loglik_damocles = fit$loglik)
}
ml <- do.call(rbind, parallel::mclapply(seq_len(10), ml_one, mc.cores = cores))
cat("2. maximum likelihood, 10 communities on trees of 60 tips (true gamma = mu = 1)\n")
print(round(ml, 4))
cat("   max |loglik difference|:", format(max(abs(ml[, "loglik"] - ml[, "loglik_damocles"])),
                                       digits = 3),
    "; ours higher or equal in", sum(ml[, "loglik"] >= ml[, "loglik_damocles"] - 1e-6), "of 10\n")

# 3. simulation ------------------------------------------------------------------------------
tree <- treeshapesim::sim_tree_shape(40)
height <- max(ape::node.depth.edgelength(tree))
distances <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
n_sim <- 4000
sim_compare <- function(gamma, mu) {
  ours <- simulate_damocles(tree, gamma, mu, n = n_sim, root = "both")[, tree$tip.label]
  roots <- rep(c(0, 1), length.out = n_sim)
  theirs <- do.call(rbind, parallel::mclapply(roots, function(root) {
    run <- DAMOCLES_sim(tree, gamma_0 = gamma, gamma_td = 0, mu = mu, sigma = 0, psiBranch = 0,
                        psiTrait = 0, z = 0, phi = 0, traitOpt = 0, br0 = 0, br_td = 0,
                        nTdim = 1, root.state = root, root.trait.state = 0)[[1]]
    stats::setNames(run$state, tree$tip.label[run$d])[tree$tip.label]
  }, mc.cores = cores))
  summary_of <- function(m) {
    c(richness = mean(rowSums(m)), mpd = mean(mean_pairwise(m, distances), na.rm = TRUE),
      mntd = mean(mean_nearest(m, distances), na.rm = TRUE))
  }
  freq <- colMeans(ours) - colMeans(theirs)
  # standard error of the difference of two frequencies
  se <- sqrt((colMeans(ours) * (1 - colMeans(ours)) + colMeans(theirs) *
                (1 - colMeans(theirs))) / n_sim)
  list(rates = c(gamma = gamma, mu = mu),
       ours = summary_of(ours), theirs = summary_of(theirs),
       max_freq_diff = max(abs(freq)), max_z = max(abs(freq / pmax(se, 1e-12))))
}
sims <- list(sim_compare(0.5 / height, 0.5 / height), sim_compare(5 / height, 5 / height),
             sim_compare(2 / height, 0.2 / height))
cat("3. simulations on a tree of 40 tips (height", round(height, 3), "),", n_sim,
    "communities each\n")
for (s in sims) {
  cat("   gamma", format(s$rates[["gamma"]], digits = 3), "mu", format(s$rates[["mu"]], digits = 3),
      "| richness", round(s$ours[["richness"]], 2), "vs", round(s$theirs[["richness"]], 2),
      "| MPD", round(s$ours[["mpd"]], 3), "vs", round(s$theirs[["mpd"]], 3),
      "| MNTD", round(s$ours[["mntd"]], 3), "vs", round(s$theirs[["mntd"]], 3),
      "| max |freq diff|", round(s$max_freq_diff, 3), "(max |z|", round(s$max_z, 2), ")\n")
}

# 4. shared rates -----------------------------------------------------------------------------
shared_one <- function(i, k) {
  tree <- treeshapesim::sim_tree_shape(60)
  fit <- fit_damocles(simulate_damocles(tree, 1, 1, n = k), tree, shared = TRUE)
  c(communities = k, gamma = fit$gamma[1], mu = fit$mu[1])
}
shared <- as.data.frame(do.call(rbind, lapply(c(1, 5, 20), function(k) {
  do.call(rbind, parallel::mclapply(seq_len(30), shared_one, k = k, mc.cores = cores))
})))
cat("4. shared rates (true gamma = mu = 1), 5%, 50% and 95% quantiles of 30 fits\n")
print(aggregate(cbind(gamma, mu) ~ communities, shared,
                function(x) round(stats::quantile(x, c(0.05, 0.5, 0.95)), 3)))

# 5. sympatry ---------------------------------------------------------------------------------
sympatry_one <- function(i, k) {
  tree <- treeshapesim::sim_tree_shape(60)
  comm <- simulate_damocles(tree, 1, 1, n = k, sympatry = 0.5)
  fit <- fit_damocles(comm, tree, shared = TRUE, sympatry = NULL)
  c(communities = k, gamma = fit$gamma[1], mu = fit$mu[1], sympatry = fit$sympatry[1])
}
sympatry <- as.data.frame(do.call(rbind, lapply(c(1, 5, 20), function(k) {
  do.call(rbind, parallel::mclapply(seq_len(30), sympatry_one, k = k, mc.cores = cores))
})))
cat("5. sympatry estimated with shared rates (true gamma = mu = 1, sympatry = 0.5),",
    "5%, 50% and 95% quantiles of 30 fits\n")
print(aggregate(cbind(gamma, mu, sympatry) ~ communities, sympatry,
                function(x) round(stats::quantile(x, c(0.05, 0.5, 0.95)), 3)))

# 6. regimes and covariates ---------------------------------------------------------------------
bm_covariate <- function(tree, variance = 0.5) {
  value <- numeric(ape::Ntip(tree) + tree$Nnode)
  out <- numeric(nrow(tree$edge))
  for (e in ape::reorder.phylo(tree, "cladewise", index.only = TRUE)) {
    value[tree$edge[e, 2]] <- value[tree$edge[e, 1]] +
      stats::rnorm(1, 0, sqrt(variance * tree$edge.length[e]))
    out[e] <- exp((value[tree$edge[e, 1]] + value[tree$edge[e, 2]]) / 2)
  }
  out
}
lrt_one <- function(i, setting) {
  tree <- treeshapesim::sim_tree_shape(60)
  if (setting$model == "covariate") {
    x <- bm_covariate(tree)
    z <- log(x) - mean(log(x))
    comm <- simulate_damocles(tree, exp(setting$effect * z), 1, n = setting$communities)
    alternative <- fit_damocles(comm, tree, shared = TRUE, edge_covariate = x, vary = "gamma")
  } else {
    clade <- ape::extract.clade(tree, ape::Ntip(tree) + 2)$tip.label
    reg <- clade_regimes(tree, list(A = clade))
    gamma <- if (setting$effect == 0) 1 else ifelse(reg == "A", 3, 0.5)
    comm <- simulate_damocles(tree, gamma, 1, n = setting$communities)
    alternative <- fit_damocles(comm, tree, shared = TRUE, regimes = reg, vary = "gamma")
  }
  damocles_lrt(fit_damocles(comm, tree, shared = TRUE), alternative)$p_value[2]
}
lrt_settings <- expand.grid(model = c("covariate", "regimes"), effect = c(0, 1),
                            communities = c(5, 20), stringsAsFactors = FALSE)
lrt_settings$rejected <- vapply(seq_len(nrow(lrt_settings)), function(s) {
  p <- unlist(parallel::mclapply(seq_len(100), lrt_one, setting = lrt_settings[s, ],
                                 mc.cores = cores))
  mean(p < 0.05)
}, numeric(1))
cat("6. likelihood-ratio tests against constant rates, proportion of 100 with p < 0.05",
    "(effect 0: type I error; 1: power)\n")
print(lrt_settings)

# 7. partition test ----------------------------------------------------------------------------
# another clade like `tested` overlapping less than a quarter of it, or NULL if 100 draws fail
elsewhere_clade <- function(tree, tested) {
  for (attempt in seq_len(100)) {
    other <- random_regimes(tree, tested)
    if (mean(other[tested == "A"] == "A") < 0.25) return(other)
  }
  NULL
}
partition_one <- function(i, scenario) {
  repeat {  # a new tree when no clade elsewhere fits (scenario "elsewhere" only)
    tree <- treeshapesim::sim_tree_shape(60)
    size <- damocles_subtree_size(tree)
    edge <- which.min(abs(size - nrow(tree$edge) / 4))
    tested <- clade_regimes(tree, list(A = tree$edge[edge, 2]), stem = TRUE)
    true <- switch(scenario, none = NULL, clade = tested,
                   elsewhere = elsewhere_clade(tree, tested))
    if (scenario != "elsewhere" || !is.null(true)) break
  }
  gamma <- if (is.null(true)) 1 else ifelse(true == "A", 3, 0.5)
  comm <- simulate_damocles(tree, gamma, 1, n = 10)
  res <- damocles_partition_test(comm, tree, tested, n_rand = 49, vary = "gamma")
  c(chisq = res$p_chisq < 0.05, random = res$p_random < 0.05)
}
partition <- t(vapply(c("none", "clade", "elsewhere"), function(scenario) {
  out <- colMeans(do.call(rbind, parallel::mclapply(seq_len(40), partition_one,
                                                    scenario = scenario, mc.cores = cores)))
  message(format(Sys.time(), "%H:%M"), " section 7, scenario ", scenario, " done")
  out
}, numeric(2)))
cat("7. proportion of 40 replicates rejecting at 0.05: chi-squared against constant rates and",
    "test against random partitions\n")
print(partition)

saveRDS(list(loglik_diffs = diffs, ml = ml, sims = sims, shared = shared, sympatry = sympatry,
             lrt = lrt_settings, partition = partition),
        "data-raw/validation_damocles.rds")
