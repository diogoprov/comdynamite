# Spatial dependence among communities in the DAMOCLES paper (dev/plano-artigo-damocles.md, S5):
# does the test of colonisation (gamma) and local extinction (mu) among regimes of the tree stay
# calibrated when the communities are not independent realisations of one process?
# Trees as in S2 (data-raw/s2_power.R): sim_shift_tree(112, 1, 2) redrawn until the focal regime
# holds 0.25 of the tips (factor 0.8-1.2), height one, gamma = 5 and mu = 1 in the background.
# Two kinds of dependence:
# - "history": the communities fall in blocks (n = 14: two blocks of seven, as the two regions of
#   Benicio et al. 2021; n = 40: four blocks of ten). Within a block they share the states of
#   the lineages (and the outcome of each speciation) from the root down to the separation age
#   a_s, then evolve independently. a_s = Inf (older than the root) is independence. On one of
#   these trees the correlation of tip states between two communities is 0.11 between blocks;
#   within a block it is 0.14, 0.25, 0.48 and about 0.7 for a_s = 0.2, 0.1, 0.05 and 0.02
#   (scratchpad check of simulate_blocks(); the rule of speciation after the split makes it
#   lower than exp(-(gamma + mu) a_s)).
# - "regions": two groups of communities (half each) with different rates in the whole tree,
#   mu multiplied by h = 3 or 10 in the second group (lower occupancy), as the inland and
#   coastal areas of Benicio et al. (data-raw/s5_benicio_spatial.R). Besides the tests of
#   S2 layer C, a stratified test: background rates free per group, the regime effect free per
#   group (2 df for gamma or mu), the log-likelihoods summed over the groups.
# Effect: none (size) or gamma x 1.5 in the focal regime (power; history only). On each data set
# the four models of S2 layer C (constant, gamma free per regime, mu free, both free) with the
# LRT p-values against the constant model.
# REPS replicates per setting (default 200), WORKERS processes (default 16), JOBS (default
# data-raw/s5_jobs, ignored by git), OUT (default data-raw/s5_spatial.rds).
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
workers <- as.integer(Sys.getenv("WORKERS", "16"))
reps <- as.integer(Sys.getenv("REPS", "200"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s5_jobs")
out <- Sys.getenv("OUT", "data-raw/s5_spatial.rds")
RNGkind("L'Ecuyer-CMRG")
dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)

quiet <- function(x) suppressWarnings(x)
focal_tree <- function(n_tips, fraction, max_tries = 1000) {
  for (i in seq_len(max_tries)) {
    sim <- sim_shift_tree(n_tips, 1, 2)
    f <- mean(sim$tip_regime > 0)
    if (f >= 0.8 * fraction && f <= 1.2 * fraction) return(sim)
  }
  stop("no tree with the focal fraction in ", max_tries, " tries", call. = FALSE)
}

# simulate_damocles() (allopatric speciation, root "both") with communities in blocks that share
# every state and speciation outcome above the separation age `a_s` (ultrametric tree).
simulate_blocks <- function(tree, gamma, mu, block, a_s) {
  rates <- damocles_edge_rates(tree, gamma, mu)
  index <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  tree <- ape::reorder.phylo(tree, "cladewise")
  gamma <- rates$gamma[index]
  mu <- rates$mu[index]
  n <- length(block)
  n_tip <- ape::Ntip(tree)
  depth <- ape::node.depth.edgelength(tree)
  age <- max(depth) - depth
  first <- match(sort(unique(block)), block)
  expand <- function(x) x[block]
  state <- matrix(0L, n, n_tip + tree$Nnode)
  state[, n_tip + 1] <- if (age[n_tip + 1] >= a_s) expand(stats::rbinom(length(first), 1, 0.5)) else
    stats::rbinom(n, 1, 0.5)
  heir <- matrix(sample.int(2, n * tree$Nnode, replace = TRUE), n, tree$Nnode)
  shared_node <- age[n_tip + seq_len(tree$Nnode)] >= a_s
  heir[, shared_node] <- heir[first, shared_node, drop = FALSE][block, , drop = FALSE]
  position <- stats::ave(seq_len(nrow(tree$edge)), tree$edge[, 1], FUN = seq_along)
  for (e in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[e, 1]
    child <- tree$edge[e, 2]
    begin <- state[, parent] * (heir[, parent - n_tip] == position[e])
    if (age[child] >= a_s) {
      state[, child] <- expand(damocles_step(begin[first], tree$edge.length[e], gamma[e], mu[e]))
    } else if (age[parent] > a_s) {
      mid <- expand(damocles_step(begin[first], age[parent] - a_s, gamma[e], mu[e]))
      state[, child] <- damocles_step(mid, a_s - age[child], gamma[e], mu[e])
    } else {
      state[, child] <- damocles_step(begin, tree$edge.length[e], gamma[e], mu[e])
    }
  }
  out <- state[, seq_len(n_tip), drop = FALSE]
  dimnames(out) <- list(paste0("sim", seq_len(n)), tree$tip.label)
  out
}

fit_tests <- function(comm, tree, regimes) {
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  fit <- function(vary) quiet(damocles_fit(comm, tree, shared = TRUE, regimes = regimes,
                                           vary = vary, start = start))
  both <- fit("both")
  g <- fit("gamma")
  m <- fit("mu")
  ll <- function(x) attr(x, "loglik_total")
  p <- function(alt, df) stats::pchisq(pmax(2 * (ll(alt) - ll(null)), 0), df, lower.tail = FALSE)
  list(ll = c(null = ll(null), gamma = ll(g), mu = ll(m), both = ll(both)),
       row = data.frame(gamma_null = null$gamma[1], mu_null = null$mu[1],
                        gamma_ratio = both$gamma_1[1] / both$gamma_0[1],
                        mu_ratio = both$mu_1[1] / both$mu_0[1],
                        p_both = p(both, 2), p_gamma = p(g, 1), p_mu = p(m, 1),
                        convergence = max(null$convergence, both$convergence, g$convergence,
                                          m$convergence)))
}

one <- function(s) {
  sim <- focal_tree(112, 0.25)
  focal <- sim$regimes != "0"
  gamma <- rep(5, length(focal))
  mu <- rep(1, length(focal))
  if (s$effect == "gamma") gamma[focal] <- 5 * 1.5
  if (s$kind == "history") {
    block <- rep(seq_len(s$n_comm / s$block_size), each = s$block_size)
    comm <- simulate_blocks(sim$tree, gamma, mu, block, s$a_s)
    r <- fit_tests(comm, sim$tree, sim$regimes)$row
    r$p_strat_gamma <- r$p_strat_mu <- r$p_strat_both <- NA_real_
  } else {
    half <- s$n_comm / 2
    group <- rep(1:2, each = half)
    comm <- rbind(simulate_damocles(sim$tree, gamma, mu, n = half),
                  simulate_damocles(sim$tree, gamma, mu * s$h, n = half))
    rownames(comm) <- paste0("sim", seq_len(nrow(comm)))
    r <- fit_tests(comm, sim$tree, sim$regimes)$row
    lls <- lapply(1:2, function(k) fit_tests(comm[group == k, , drop = FALSE], sim$tree,
                                             sim$regimes)$ll)
    tot <- Reduce(`+`, lls)
    p <- function(alt, df) stats::pchisq(pmax(2 * (tot[[alt]] - tot[["null"]]), 0), df,
                                         lower.tail = FALSE)
    r$p_strat_both <- p("both", 4)
    r$p_strat_gamma <- p("gamma", 2)
    r$p_strat_mu <- p("mu", 2)
  }
  r$richness <- mean(rowSums(comm))
  r$similarity <- local({
    both <- tcrossprod(comm)
    rich <- rowSums(comm)
    mean((both / (outer(rich, rich, "+") - both))[lower.tri(both)])
  })
  r
}

settings <- rbind(
  data.frame(expand.grid(kind = "history", n_comm = c(14, 40), a_s = c(Inf, 0.2, 0.1, 0.05, 0.02),
                         effect = c("none", "gamma"), stringsAsFactors = FALSE), h = 1),
  data.frame(expand.grid(kind = "regions", n_comm = c(14, 40), a_s = Inf, effect = "none",
                         stringsAsFactors = FALSE), h = rep(c(3, 10), each = 2)))
settings$block_size <- ifelse(settings$n_comm == 14, 7, 10)
only <- Sys.getenv("SETTINGS", "")
run <- seq_len(nrow(settings))
if (nzchar(only)) run <- intersect(run, as.integer(strsplit(only, ",")[[1]]))
jobs <- expand.grid(setting = run, rep = seq_len(reps))
jobs$file <- file.path(jobs_dir, sprintf("s%02d_r%04d.rds", jobs$setting, jobs$rep))
todo <- which(!file.exists(jobs$file))
todo <- todo[order(-settings$n_comm[jobs$setting[todo]])]
cat("S5 spatial:", length(todo), "of", nrow(jobs), "replicates on", workers, "workers\n")
started <- Sys.time()
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  set.seed(505 + 10000 * jobs$setting[j] + jobs$rep[j])
  t0 <- Sys.time()
  r <- tryCatch(one(s), error = function(e) data.frame(error = conditionMessage(e)))
  if (is.null(r$error)) r$error <- NA_character_
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(cbind(s, rep = jobs$rep[j], r, row.names = NULL), jobs$file[j])
  NULL
}, mc.cores = workers, mc.preschedule = FALSE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
done <- file.exists(jobs$file[jobs$setting %in% run])
if (!all(done)) stop(sum(!done), " replicates missing in ", jobs_dir, call. = FALSE)
cols <- c("gamma_null", "mu_null", "gamma_ratio", "mu_ratio", "p_both", "p_gamma", "p_mu",
          "p_strat_both", "p_strat_gamma", "p_strat_mu", "convergence", "richness",
          "similarity", "error", "seconds")
results <- do.call(rbind, lapply(jobs$file[jobs$setting %in% run], function(f) {
  x <- readRDS(f)
  for (k in setdiff(cols, names(x))) x[[k]] <- NA
  x[, c(names(settings), "rep", cols)]
}))

keys <- c("kind", "n_comm", "a_s", "h", "effect")
rate <- function(p) mean(p < 0.05, na.rm = TRUE)
summary <- do.call(rbind, lapply(split(results, results[, keys], drop = TRUE), function(d)
  data.frame(d[1, keys], p_both = rate(d$p_both), p_gamma = rate(d$p_gamma),
             p_mu = rate(d$p_mu), strat_both = rate(d$p_strat_both),
             strat_gamma = rate(d$p_strat_gamma), strat_mu = rate(d$p_strat_mu),
             gamma_ratio = stats::median(d$gamma_ratio, na.rm = TRUE),
             similarity = round(mean(d$similarity, na.rm = TRUE), 3),
             richness = round(mean(d$richness, na.rm = TRUE), 1),
             errors = sum(!is.na(d$error)), not_converged = sum(d$convergence != 0, na.rm = TRUE))))
summary <- summary[order(summary$kind, summary$effect, summary$n_comm, summary$h,
                         -summary$a_s), ]
cat("S5 spatial:", nrow(results), "replicates in", round(elapsed, 1), "min\n")
options(width = 200)
print(summary, row.names = FALSE, digits = 3)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
