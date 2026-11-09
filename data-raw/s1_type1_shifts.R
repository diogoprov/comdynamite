# S1 of the DAMOCLES paper (dev/plano-artigo-damocles.md): type I error of the tests of
# colonisation and local extinction among regimes when the tree has shifts in the rate of
# speciation but the communities assemble with constant rates.
# 1. Trees: a forward birth-death simulation with heritable rates, run until `n_tips` lineages
#    are alive. Speciation starts at 1 with relative extinction 0.5; each of `n_shifts` shifts
#    multiplies the speciation rate of one living lineage and all its descendants by `ratio`
#    (extinction keeps the ratio 0.5). Shift k happens when the number of living lineages first
#    reaches a number drawn between 3 and n_tips / 4, on a lineage drawn at random. Trees with a
#    regime (the background included) holding fewer than 5% of the tips are drawn again. Extinct lineages are pruned
#    and the tree rescaled to height one. An edge takes the regime of its lower end.
# 2. Communities: simulate_damocles() with gamma = mu = 1 on every edge (no effect of the
#    regimes), `n_comm` communities per tree.
# 3. Tests, likelihood-ratio against constant rates shared by the communities (damocles_lrt()):
#    - true regimes, gamma and mu free per regime (vary = "both", df = 2 n_shifts);
#    - true regimes, gamma only (df = n_shifts);
#    - DR per edge (dr_rates()) as covariate of gamma (df = 1);
#    and, in the first `n_partition` replicates of each setting, damocles_partition_test() with
#    the true regimes (vary = "both", 49 random partitions).
# Settings: n_tips 112 or 250, n_shifts 1 or 3, ratio 2 or 4, n_comm 14 or 40, plus a control
# with ratio 1 (shifts with no change in rate). REPS replicates per setting (default 500),
# PARTITION replicates with the partition test (default 50), CORES workers (default 6).
# Results in data-raw/s1_type1_shifts.rds; the summary gives the rejection rate at 5% with a
# 95% binomial interval (Clopper-Pearson).
suppressMessages(devtools::load_all(".", quiet = TRUE))
cores <- as.integer(Sys.getenv("CORES", "6"))
reps <- as.integer(Sys.getenv("REPS", "500"))
n_partition <- as.integer(Sys.getenv("PARTITION", "50"))
out <- Sys.getenv("OUT", "data-raw/s1_type1_shifts.rds")
RNGkind("L'Ecuyer-CMRG")
set.seed(101)

# 1. Trees with shifts -------------------------------------------------------------------------
sim_shift_tree <- function(n_tips, n_shifts, ratio, lambda0 = 1, epsilon = 0.5,
                           min_fraction = 0.05) {
  repeat {
    tree <- try_shift_tree(n_tips, n_shifts, ratio, lambda0, epsilon)
    if (is.null(tree)) next
    shifted <- table(factor(tree$tip_regime, levels = seq_len(n_shifts + 1) - 1))
    if (all(shifted >= min_fraction * n_tips)) return(tree)
  }
}
try_shift_tree <- function(n_tips, n_shifts, ratio, lambda0, epsilon) {
  thresholds <- sort(sample(3:max(4, floor(n_tips / 4)), n_shifts, replace = TRUE))
  # Lineages: parent node, start time, rate, regime; nodes numbered as created.
  parent <- 0L
  start <- 0
  rate <- lambda0
  regime <- 0L
  alive <- 1L
  # Edges as (parent node, child node, length, regime of the child end), filled as lineages end.
  edges <- matrix(0, 0, 4)
  tip_status <- integer(0)
  next_node <- 1L
  node_of_lineage <- integer(0)
  t <- 0
  shifts_done <- 0L
  end_lineage <- function(i, time, extinct) {
    child <- next_node
    next_node <<- next_node + 1L
    edges <<- rbind(edges, c(parent[i], child, time - start[i], regime[i]))
    tip_status[child] <<- if (extinct) 2L else 0L
    child
  }
  while (length(alive) < n_tips) {
    while (shifts_done < n_shifts && length(alive) >= thresholds[shifts_done + 1]) {
      shifts_done <- shifts_done + 1L
      i <- alive[sample.int(length(alive), 1)]
      # The shift splits the lineage: the part above keeps the old regime.
      node <- end_lineage(i, t, FALSE)
      tip_status[node] <- 1L
      parent[i] <- node
      start[i] <- t
      rate[i] <- rate[i] * ratio
      regime[i] <- shifts_done
    }
    total <- sum(rate[alive] * (1 + epsilon))
    t <- t + stats::rexp(1, total)
    i <- alive[sample.int(length(alive), 1, prob = rate[alive])]
    if (stats::runif(1) < 1 / (1 + epsilon)) {
      node <- end_lineage(i, t, FALSE)
      tip_status[node] <- 1L
      new <- length(parent) + 1:2
      parent[new] <- node
      start[new] <- t
      rate[new] <- rate[i]
      regime[new] <- regime[i]
      alive <- c(setdiff(alive, i), new)
    } else {
      end_lineage(i, t, TRUE)
      alive <- setdiff(alive, i)
      if (!length(alive)) return(NULL)
    }
  }
  if (shifts_done < n_shifts) return(NULL)
  t <- t + stats::rexp(1, sum(rate[alive] * (1 + epsilon)))
  for (i in alive) end_lineage(i, t, FALSE)
  to_phylo(edges, tip_status)
}
# From the edge list (with nodes of one child at shifts) to a phylo object with the regime of
# each edge; extinct tips and nodes of one child are removed.
to_phylo <- function(edges, tip_status) {
  children <- table(factor(edges[, 1], levels = seq_along(tip_status)))
  is_tip <- children[seq_along(tip_status)] == 0
  is_tip <- as.vector(is_tip) & seq_along(tip_status) %in% edges[, 2]
  tips <- which(is_tip)
  internal <- setdiff(c(0, edges[, 2]), tips)
  # Root is node 0 (the start of the first lineage).
  id <- integer(max(edges[, 2]) + 1)
  id[tips + 1] <- seq_along(tips)
  id[internal + 1] <- length(tips) + seq_along(internal)
  tree <- list(edge = cbind(id[edges[, 1] + 1], id[edges[, 2] + 1]), edge.length = edges[, 3],
               tip.label = paste0("t", seq_along(tips)), Nnode = length(internal))
  class(tree) <- "phylo"
  node_regime <- integer(length(tips) + length(internal))
  node_regime[tree$edge[, 2]] <- edges[, 4]
  tree$node.label <- as.character(node_regime[length(tips) + seq_along(internal)])
  tree$tip.label <- paste0("t", seq_along(tips), "_r", node_regime[seq_along(tips)])
  extinct <- tree$tip.label[tip_status[tips] == 2]
  tree <- ape::drop.tip(tree, extinct, collapse.singles = FALSE)
  tree <- ape::collapse.singles(tree)
  # Regime of each edge: of its lower end (tip label or node label).
  lower <- tree$edge[, 2]
  n_tip <- ape::Ntip(tree)
  regime <- ifelse(lower <= n_tip, sub(".*_r", "", tree$tip.label[pmin(lower, n_tip)]),
                   tree$node.label[pmax(lower - n_tip, 1)])
  height <- max(ape::node.depth.edgelength(tree))
  tree$edge.length <- tree$edge.length / height
  tip_regime <- as.integer(sub(".*_r", "", tree$tip.label))
  tree$node.label <- NULL
  tree$tip.label <- paste0("t", seq_len(n_tip))
  list(tree = tree, regimes = factor(regime), tip_regime = tip_regime)
}

# 2-3. One replicate ---------------------------------------------------------------------------
p_value <- function(null, alternative) damocles_lrt(null, alternative)$p_value[2]
one <- function(setting, partition) {
  sim <- sim_shift_tree(setting$n_tips, setting$n_shifts, setting$ratio)
  tree <- sim$tree
  comm <- simulate_damocles(tree, 1, 1, n = setting$n_comm)
  quiet <- function(x) suppressWarnings(x)
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  both <- quiet(damocles_fit(comm, tree, shared = TRUE, regimes = sim$regimes, vary = "both"))
  gamma <- quiet(damocles_fit(comm, tree, shared = TRUE, regimes = sim$regimes, vary = "gamma"))
  dr <- quiet(damocles_fit(comm, tree, shared = TRUE, edge_covariate = dr_rates(tree, "edge"),
                           vary = "gamma"))
  result <- data.frame(
    p_both = p_value(null, both), p_gamma = p_value(null, gamma), p_dr = p_value(null, dr),
    p_partition = NA_real_, shifted_tips = mean(sim$tip_regime > 0),
    shifted_edges = mean(sim$regimes != "0"),
    convergence = max(null$convergence, both$convergence, gamma$convergence, dr$convergence))
  if (partition) {
    test <- quiet(damocles_partition_test(comm, tree, sim$regimes, n_rand = 49, vary = "both"))
    result$p_partition <- test$p_random
  }
  result
}

settings <- expand.grid(n_tips = c(112, 250), n_shifts = c(1, 3), ratio = c(2, 4),
                        n_comm = c(14, 40))
settings <- rbind(settings, data.frame(n_tips = c(112, 250), n_shifts = 1, ratio = 1, n_comm = 14))
jobs <- expand.grid(setting = seq_len(nrow(settings)), rep = seq_len(reps))
started <- Sys.time()
results <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  s <- settings[jobs$setting[j], ]
  r <- tryCatch(one(s, jobs$rep[j] <= n_partition), error = function(e) {
    data.frame(p_both = NA, p_gamma = NA, p_dr = NA, p_partition = NA, shifted_tips = NA,
               shifted_edges = NA, convergence = NA, error = conditionMessage(e))
  })
  if (is.null(r$error)) r$error <- NA_character_
  cbind(s, rep = jobs$rep[j], r)
}, mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

rate <- function(p) {
  p <- p[!is.na(p)]
  k <- sum(p < 0.05)
  ci <- if (length(p)) stats::binom.test(k, length(p))$conf.int else c(NA, NA)
  sprintf("%.3f [%.3f, %.3f] (n=%d)", k / max(length(p), 1), ci[1], ci[2], length(p))
}
summary <- do.call(rbind, lapply(split(results, results[, names(settings)], drop = TRUE),
  function(d) {
    cbind(d[1, names(settings)], both = rate(d$p_both), gamma = rate(d$p_gamma),
          dr = rate(d$p_dr), partition = rate(d$p_partition),
          shifted_tips = round(mean(d$shifted_tips, na.rm = TRUE), 2),
          errors = sum(!is.na(d$error)), not_converged = sum(d$convergence != 0, na.rm = TRUE))
  }))
summary <- summary[order(summary$ratio, summary$n_tips, summary$n_shifts, summary$n_comm), ]
cat("S1:", nrow(results), "replicates in", round(elapsed, 1), "min\n")
print(summary, row.names = FALSE)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
