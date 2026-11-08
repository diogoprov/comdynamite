# Herodotools against known histories. Run from the root of the package with Herodotools
# installed (GitHub: GabrielNakamura/Herodotools; it needs daee and BioGeoBEARS from GitHub).
# 1. Equivalence: calc_tip_based_trait_evo() with the true history in place of its stochastic
#    maps gives the same values as tip_history_metrics().
# 2. Recovery of the tip-based metrics from the tips alone, by its own stochastic maps (SYM model,
#    Q by MCMC, 20 maps averaged per tip), at three rates of change per tree height.
# 3. Age of arrival: calc_age_arrival() with the true ancestral areas and with areas
#    reconstructed by maximum likelihood (ape::ace, ER model, most likely state at each node),
#    against the true ages of arrival from arrival_ages(), in a metacommunity simulated with
#    region_pool().
# CORES sets the number of workers.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(Herodotools)
})
cores <- as.integer(Sys.getenv("CORES", "6"))
plain_tree <- function(h) structure(h[c("edge", "edge.length", "tip.label", "Nnode")],
                                    class = "phylo")

# 1. equivalence
with_history <- Herodotools::calc_tip_based_trait_evo
code <- paste(deparse(body(with_history)), collapse = "\n")
code <- sub("result.scm <- phytools::make.simmap[(][^)]*[)]", "result.scm <- truth", code)
stopifnot(grepl("result.scm <- truth", code))
body(with_history) <- parse(text = code)[[1]]
formals(with_history) <- c(formals(with_history), alist(truth = ))
set.seed(11)
equivalence <- do.call(rbind, lapply(1:30, function(r) {
  h <- trait_mk(treeshapesim::sim_tree_shape(25), q = stats::runif(1, 0.2, 2),
                states = c("A", "B", "C"))
  ours <- tip_history_metrics(h)
  theirs <- with_history(h, h$states, nsim = 1, truth = h)[[1]][rownames(ours), ]
  data.frame(transitions = isTRUE(all.equal(ours$transitions, c(theirs$transitions))),
             total_nodes = isTRUE(all.equal(ours$total_nodes, c(theirs$total.nodes))),
             last_transition_time = isTRUE(all.equal(ours$last_transition_time,
                                                     c(theirs$last.transition.time))),
             stasis_time = isTRUE(all.equal(ours$stasis_time, c(theirs$stasis.time))))
}))
cat("equivalence, histories where all tips agree:\n")
print(colSums(equivalence))

# 2 and 3, one replicate
rates <- c(0.3, 1, 3)
one_rep <- function(rate, seed) {
  set.seed(seed)
  tree <- treeshapesim::sim_tree_shape(40)
  repeat {
    h <- trait_mk(tree, q = rate, states = c("A", "B"))
    if (length(unique(h$states)) == 2) break
  }
  truth <- tip_history_metrics(h)
  maps <- suppressMessages(calc_tip_based_trait_evo(tree, h$states, nsim = 20))
  estimate <- Reduce(`+`, lapply(maps, function(m) {
    as.matrix(data.frame(transitions = c(m$transitions), last = c(m$last.transition.time),
                         stasis = c(m$stasis.time), row.names = rownames(m))[rownames(truth), ])
  })) / length(maps)
  tips <- data.frame(rate = rate, seed = seed, species = rownames(truth),
                     transitions = truth$transitions, changes = truth$changes,
                     last = truth$last_transition_time, stasis = truth$stasis_time,
                     transitions_hat = estimate[, "transitions"], last_hat = estimate[, "last"],
                     stasis_hat = estimate[, "stasis"])
  regions <- rep(c("A", "B"), each = 10)
  sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                                competition_matrix(tree, "none"), n_steps = 200,
                                pool = region_pool(h, regions))
  input <- suppressWarnings(as_herodotools(sim, regions, history = h))
  true_age <- arrival_ages(h, sim, regions)[rownames(input$W), colnames(input$W)]
  fit <- ape::ace(h$states[tree$tip.label], plain_tree(h), type = "discrete", model = "ER")
  reconstructed <- input
  reconstructed$ancestral.area$area <- colnames(fit$lik.anc)[max.col(fit$lik.anc, "first")]
  with_truth <- calc_age_arrival(input$W, input$tree, input$ancestral.area, input$biogeo)
  with_ml <- calc_age_arrival(reconstructed$W, reconstructed$tree, reconstructed$ancestral.area,
                              reconstructed$biogeo)
  present <- input$W == 1
  ages <- data.frame(rate = rate, seed = seed, true = true_age[present],
                     nodes_true = with_truth$age_arrival_assemblage[present],
                     nodes_ml = with_ml$age_arrival_assemblage[present],
                     nodes_correct = mean(reconstructed$ancestral.area$area ==
                                            input$ancestral.area$area))
  list(tips = tips, ages = ages)
}
jobs <- expand.grid(rate = rates, seed = 1:20)
future::plan(future::multicore, workers = cores)
started <- Sys.time()
runs <- future.apply::future_Map(one_rep, jobs$rate, jobs$seed + 1000, future.seed = NULL)
future::plan(future::sequential)
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
tips <- do.call(rbind, lapply(runs, `[[`, "tips"))
ages <- do.call(rbind, lapply(runs, `[[`, "ages"))

recovery <- do.call(rbind, lapply(split(tips, tips$rate), function(d) {
  do.call(rbind, lapply(c("transitions", "last", "stasis"), function(m) {
    data.frame(rate = d$rate[1], metric = m,
               spearman = stats::cor(d[[m]], d[[paste0(m, "_hat")]], method = "spearman"),
               bias = mean(d[[paste0(m, "_hat")]] - d[[m]]), truth_mean = mean(d[[m]]))
  }))
}))
missed <- stats::aggregate(cbind(changes, transitions) ~ rate, tips, mean)
arrival <- do.call(rbind, lapply(split(ages, ages$rate), function(d) {
  data.frame(rate = d$rate[1], nodes_correct = mean(d$nodes_correct),
             spearman_true_areas = stats::cor(d$true, d$nodes_true, method = "spearman"),
             bias_true_areas = mean(d$nodes_true - d$true),
             spearman_ml_areas = stats::cor(d$true, d$nodes_ml, method = "spearman"),
             bias_ml_areas = mean(d$nodes_ml - d$true), mean_true = mean(d$true))
}))
rownames(recovery) <- rownames(arrival) <- rownames(missed) <- NULL
saveRDS(list(equivalence = equivalence, recovery = recovery, missed = missed, arrival = arrival,
             tips = tips, ages = ages, minutes = minutes,
             versions = c(Herodotools = as.character(utils::packageVersion("Herodotools")),
                          phytools = as.character(utils::packageVersion("phytools")))),
        "data-raw/validation_herodotools.rds")
cat("minutes:", round(minutes, 1), "\n")
print(recovery, digits = 2)
print(missed, digits = 2)
print(arrival, digits = 2)
