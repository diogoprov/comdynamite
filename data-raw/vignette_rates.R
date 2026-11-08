# Results for the vignette "Rates of trait evolution that vary among branches"
# (vignettes/branch-rates.Rmd). Run from the root of the package with RRphylo and mvMORPH
# installed (both GPL; they are only called here). CORES sets the number of workers.
# 1. RRphylo (Castiglione et al. 2018) under the three models of trait_rates(): Spearman
#    correlation of its rate per branch with the change that happened on the branch (per unit
#    length), of its absolute value with the true rate, and the
#    ceiling given by the true change on each branch (one change per branch is all any method
#    sees).
# 2. Regimes: traits with one rate per habitat from history_rates("state") on a history from
#    trait_mk(), and the ratio of the two rates estimated by mvMORPH::mvBM(model = "BMM") on the
#    true history (O'Meara et al. 2006).
# 3. Turnover: traits with rates from turnover_rates() on a simulated metacommunity, and the
#    Spearman correlation between the absolute rate RRphylo estimates on the tip branches and
#    the turnover of partners of each species.
suppressMessages(devtools::load_all(".", quiet = TRUE))
cores <- as.integer(Sys.getenv("CORES", "6"))
run <- function(n, f) do.call(rbind, parallel::mclapply(seq_len(n), f, mc.cores = cores,
                                                         mc.set.seed = TRUE))
RNGkind("L'Ecuyer-CMRG")
spearman <- function(a, b) stats::cor(a, b, method = "spearman")

# RRphylo rate of each branch, in the order of the rows of tree$edge
rrphylo_rates <- function(tree, tips) {
  rr <- suppressMessages(suppressWarnings(RRphylo::RRphylo(tree, tips, clus = 0)))
  n <- ape::Ntip(tree)
  node <- tree$edge[, 2]
  label <- ifelse(node <= n, tree$tip.label[pmin(node, n)], as.character(node))
  rr$rates[label, 1]
}

# 1. RRphylo under the three models
set.seed(31)
models <- do.call(rbind, lapply(c("iid", "gbm", "shift"), function(model) {
  run(30, function(r) {
    tree <- treeshapesim::sim_tree_shape(60)
    x <- trait_rates(tree, model, rate_sd = 0.8, shifts = 2, shift_size = c(8, 1 / 8))
    estimate <- rrphylo_rates(tree, x$tips)
    n <- ape::Ntip(tree)
    values <- stats::setNames(c(x$tips[tree$tip.label], x$nodes), c(seq_len(n), names(x$nodes)))
    change <- values[as.character(tree$edge[, 2])] - values[as.character(tree$edge[, 1])]
    data.frame(model = model, with_change = spearman(estimate, change / tree$edge.length),
               estimated = spearman(abs(estimate), x$edge_rate),
               ceiling = spearman(abs(change) / sqrt(tree$edge.length), x$edge_rate))
  })
}))
cat("1. RRphylo, mean Spearman rho with the true rate per branch:\n")
print(aggregate(cbind(with_change, estimated, ceiling) ~ model, models,
                function(v) round(mean(v), 2)))

# 2. Regimes
set.seed(32)
regimes <- do.call(rbind, lapply(c(1, 2, 5), function(ratio) {
  run(30, function(r) {
    repeat {  # both habitats on enough branches
      tree <- treeshapesim::sim_tree_shape(60)
      habitat <- trait_mk(tree, q = 0.5, states = c("forest", "open"))
      time <- tapply(unlist(habitat$maps), names(unlist(habitat$maps)), sum)
      if (length(time) == 2 && min(time) / sum(time) > 0.2) break
    }
    rates <- history_rates(habitat, "state", state_rates = c(forest = 1, open = ratio))
    x <- trait_rates(habitat, edge_rates = rates)
    fit <- mvMORPH::mvBM(habitat, x$tips[habitat$tip.label], model = "BMM", echo = FALSE,
                         diagnostic = FALSE)
    sigma <- fit$sigma[1, 1, ]
    data.frame(ratio = ratio, estimated = unname(sigma["open"] / sigma["forest"]))
  })
}))
cat("2. Regimes, estimated ratio of rates (open / forest):\n")
print(aggregate(estimated ~ ratio, regimes, function(v) round(stats::median(v), 2)))

# 3. Turnover of partners
set.seed(33)
turnover <- run(30, function(r) {
  tree <- treeshapesim::sim_tree_shape(40)
  sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                                competition_matrix(tree, "none"), n_steps = 400,
                                record_steps = seq(200, 400, 10), extirpation = 0.02)
  rates <- turnover_rates(sim, k = 100)
  x <- trait_rates(tree, edge_rates = rates)
  estimate <- rrphylo_rates(tree, x$tips)
  tip_edge <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
  data.frame(rho = spearman(abs(estimate[tip_edge]), attr(rates, "turnover")[tree$tip.label]),
             spread = stats::sd(rates[tip_edge]) / mean(rates[tip_edge]))
})
cat("3. Turnover, Spearman rho between the RRphylo tip rate and the turnover:\n")
print(round(summary(turnover$rho), 2))

saveRDS(list(models = models, regimes = regimes, turnover = turnover),
        "vignettes/branch-rates.rds")
