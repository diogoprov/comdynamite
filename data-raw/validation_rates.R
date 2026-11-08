# Rates that vary among branches: trait_rates() and trait_mk() against the methods that estimate
# them. Run from the root of the package with RRphylo, phytools, motmot and corHMM installed
# (all GPL; they are only called here). CORES sets the number of workers.
# 1. RRphylo (Castiglione et al. 2018), under the three models of trait_rates(): Spearman
#    correlation of its rate per branch with the change that happened on the branch (per unit
#    length), and of its absolute value with the true rate. As a ceiling, the correlation of the
#    true change (per square root of length) with the true rate: one change per branch is all
#    any method sees.
# 2. phytools::multirateBM (Revell 2021): correlation between the logarithms of the estimated
#    and true rates at the nodes, under "gbm".
# 3. motmot traitMedusa (Thomas & Freckleton 2012): is the one true shift found, under "shift"?
# 4. phytools::fitgammaMk (Revell 2024): estimate of alpha under trait_mk(alpha = 0.3), and
#    under constant rates.
# 5. corHMM (Beaulieu et al. 2013): does the model with two hidden rate classes beat the model
#    with one, under trait_mk(rate_classes) and under constant rates?
suppressMessages(devtools::load_all(".", quiet = TRUE))
cores <- as.integer(Sys.getenv("CORES", "6"))
run <- function(n, f) do.call(rbind, parallel::mclapply(seq_len(n), f, mc.cores = cores,
                                                         mc.set.seed = TRUE))
RNGkind("L'Ecuyer-CMRG")

# 1. RRphylo
set.seed(21)
rrphylo <- do.call(rbind, lapply(c("iid", "gbm", "shift"), function(model) {
  run(30, function(r) {
    tree <- treeshapesim::sim_tree_shape(60)
    x <- trait_rates(tree, model, rate_sd = 0.8, shifts = 2, shift_size = c(8, 1 / 8))
    rr <- suppressMessages(suppressWarnings(RRphylo::RRphylo(tree, x$tips, clus = 0)))
    n <- ape::Ntip(tree)
    node <- tree$edge[, 2]
    label <- ifelse(node <= n, tree$tip.label[pmin(node, n)], as.character(node))
    estimate <- rr$rates[label, 1]
    values <- stats::setNames(c(x$tips[tree$tip.label], x$nodes), c(seq_len(n), names(x$nodes)))
    change <- values[as.character(tree$edge[, 2])] - values[as.character(tree$edge[, 1])]
    spearman <- function(a, b) stats::cor(a, b, method = "spearman")
    data.frame(model = model, with_change = spearman(estimate, change / tree$edge.length),
               with_rate = spearman(abs(estimate), x$edge_rate),
               ceiling = spearman(abs(change) / sqrt(tree$edge.length), x$edge_rate))
  })
}))
cat("1. RRphylo, mean Spearman rho per branch:\n")
print(aggregate(cbind(with_change, with_rate, ceiling) ~ model, rrphylo,
                function(v) round(mean(v), 2)))

# 2. multirateBM
set.seed(22)
multirate <- run(20, function(r) {
  tree <- treeshapesim::sim_tree_shape(40)
  x <- trait_rates(tree, "gbm", rate_sd = 0.8)
  fit <- tryCatch(phytools::multirateBM(tree, x$tips, lambda = 1, quiet = TRUE),
                  error = function(e) NULL)  # fails on trees with very short branches
  if (is.null(fit)) return(data.frame(r = NA_real_, slope = NA_real_))
  estimate <- log(fit$sig2)
  truth <- log(x$node_rate)
  # sig2 has the tips by label, then the internal nodes by number
  names(truth) <- c(tree$tip.label, names(truth)[-seq_len(ape::Ntip(tree))])
  truth <- truth[names(estimate)]
  data.frame(r = stats::cor(estimate, truth),
             slope = unname(stats::coef(stats::lm(estimate ~ truth))[2]))
})
cat("2. multirateBM, log rates at nodes: r", round(mean(multirate$r, na.rm = TRUE), 2), "+/-",
    round(stats::sd(multirate$r, na.rm = TRUE), 2), "; slope",
    round(mean(multirate$slope, na.rm = TRUE), 2), "; failed fits", sum(is.na(multirate$r)),
    "of", nrow(multirate), "\n")

# 3. traitMedusa
set.seed(23)
medusa <- run(30, function(r) {
  tree <- treeshapesim::sim_tree_shape(80)
  inner <- which(tree$edge[, 2] > ape::Ntip(tree))
  clade <- vapply(tree$edge[inner, 2], function(n) length(ape::extract.clade(tree, n)$tip.label),
                  integer(1))
  candidates <- inner[clade >= 8 & clade <= 40]
  edge <- candidates[sample.int(length(candidates), 1)]
  x <- trait_rates(tree, "shift", shift_edges = edge, shift_size = 8)
  y <- matrix(x$tips, dimnames = list(names(x$tips), NULL))
  invisible(utils::capture.output(m <- motmot::transformPhylo.ML(y, tree, model = "tm2",
                                                               minCladeSize = 5, nSplits = 2)))
  rates <- motmot::summary.traitMedusa(m, cutoff = 4, AICc = TRUE)$Rates
  # "Single rate" when no shift is chosen; otherwise a table of nodes and relative rates
  chosen <- if (is.data.frame(rates)) as.integer(rates$node) else integer(0)
  hit <- match(x$shifts$node, chosen)
  data.frame(found = !is.na(hit), any = length(chosen) > 0,
             multiplier = if (is.na(hit)) NA_real_ else as.numeric(rates$MLRate[hit]))
})
cat("3. traitMedusa, true shift (x8) found:", mean(medusa$found), "; any shift found:",
    mean(medusa$any), "; median estimated multiplier when found:",
    round(stats::median(medusa$multiplier, na.rm = TRUE), 1), "\n")

# 4. fitgammaMk
set.seed(24)
gamma_mk <- do.call(rbind, lapply(c(0.3, Inf), function(alpha) {
  run(12, function(r) {
    tree <- treeshapesim::sim_tree_shape(500)
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    h <- trait_mk(tree, q = 3, alpha = if (is.finite(alpha)) alpha)
    fit <- phytools::fitgammaMk(tree, h$states, nrates = 8)
    data.frame(true_alpha = alpha, alpha = fit$alpha)
  })
}))
cat("4. fitgammaMk, estimated alpha (median, quartiles):\n")
print(aggregate(alpha ~ true_alpha, gamma_mk, function(v) round(stats::quantile(v, c(.25, .5, .75)), 2)))

# 5. corHMM
set.seed(25)
hidden <- do.call(rbind, lapply(c(TRUE, FALSE), function(classes) {
  run(12, function(r) {
    tree <- treeshapesim::sim_tree_shape(300)
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    h <- if (classes) {
      trait_mk(tree, q = 2, rate_classes = c(slow = 0.05, fast = 4), class_switch = 1)
    } else {
      trait_mk(tree, q = 2)
    }
    data <- data.frame(species = tree$tip.label, state = h$states)
    fits <- lapply(1:2, function(k) {
      suppressWarnings(utils::capture.output(fit <- corHMM::corHMM(tree, data, rate.cat = k, model = "ER",
                                      node.states = "none", nstarts = 2)))
      fit
    })
    data.frame(classes = classes, delta_aicc = fits[[1]]$AICc - fits[[2]]$AICc)
  })
}))
cat("5. corHMM, two classes preferred (delta AICc > 2):\n")
print(aggregate(delta_aicc ~ classes, hidden, function(v) c(prop = mean(v > 2),
                                                             median = round(stats::median(v), 1))))

saveRDS(list(rrphylo = rrphylo, multirate = multirate, medusa = medusa, gamma_mk = gamma_mk,
             hidden = hidden), "data-raw/validation_rates.rds")
