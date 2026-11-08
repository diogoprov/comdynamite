# simulate_damocles(), damocles_loglik() and fit_damocles() against the DAMOCLES package of
# Pigot & Etienne (2015) (GPL, only called here; install.packages("DAMOCLES"), version 2.3).
# 1. Log-likelihoods against DAMOCLES::DAMOCLES_loglik() on 20 trees, three pairs of rates and
#    the three priors at the root (pchoice 0, 1, 2), plus trees with a root edge.
# 2. Maximum-likelihood rates against DAMOCLES::DAMOCLES_ML() on 10 simulated communities.
# 3. Simulations against DAMOCLES::DAMOCLES_sim() (a Gillespie simulation of every event): the
#    frequency of presence of each tip, richness, MPD and MNTD over 4000 communities each.
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

saveRDS(list(loglik_diffs = diffs, ml = ml, sims = sims), "data-raw/validation_damocles.rds")
