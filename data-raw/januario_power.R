# Power of the test of Januario et al. (2025): does a PGLS of tip speciation rate (lambda_DR)
# on metapopulation persistence (logit equilibrium occupancy, log turnover rate) detect a causal
# link between persistence and speciation?
#
#   Rscript data-raw/januario_power.R
# Results in data-raw/januario_power.rds (44 min on 6 cores).
#
# Each replicate:
# 1. a latent trait X evolves by Brownian motion while lineages speciate at lambda0 exp(b X)
#    (b = 0 is the null) and go extinct at a constant rate, until 200 lineages are alive
#    (forward simulation written here);
# 2. X of the tips, standardised, sets each species' dispersal (0.01 exp(z)) in
#    simulate_metacommunity() on a 10 x 10 grid without competition; as a negative control on a
#    subset, it sets extirpation (0.02 exp(-z)) instead, a link the pilot showed is drowned by
#    demographic turnover;
# 3. colonisation and extirpation of each species are estimated from presence over the last
#    `window` steps with a two-state Markov model, giving equilibrium occupancy
#    eps = logit(c / (c + e)) and resilience tau = log(c + e), as in Januario et al. (2025);
# 4. lambda_DR (Jetz et al. 2012) on the reconstructed tree; PGLS (nlme::gls with a Brownian
#    correlation) of log lambda_DR, and of the true log lambda, on eps and tau; p-values of each
#    term and of a likelihood-ratio test of both against the intercept.
#
# Caveat: the link from persistence to speciation is imposed (a continuous-trait SSE), not
# emergent; the experiment asks whether their design has power under a causal link, not
# whether the link exists.

devtools::load_all(quiet = TRUE)

# birth-death tree with speciation rate lambda0 * exp(b * X), X a Brownian trait; time steps of
# dt; stops when n lineages are alive; returns the reconstructed tree and X at the tips
sim_tree_trait <- function(n = 200, lambda0 = 1, mu = 0.3, b = 0, sigma = 0.5, dt = 0.01,
                           max_tries = 50) {
  for (try in seq_len(max_tries)) {
    parent <- c(1L, 1L); start <- c(0, 0); x <- c(0, 0); id <- c(1L, 2L)
    next_node <- 2L; next_lineage <- 3L
    edges <- list(); t <- 0
    while (length(id) > 0 && length(id) < n) {
      t <- t + dt
      x <- x + stats::rnorm(length(x), 0, sigma * sqrt(dt))
      u <- stats::runif(length(id))
      birth <- u < lambda0 * exp(b * x) * dt
      death <- !birth & u > 1 - mu * dt
      for (k in which(birth | death)) {
        if (birth[k]) {
          edges[[length(edges) + 1]] <- c(parent[k], -next_node, t - start[k], 0)
          parent <- c(parent, next_node, next_node); start <- c(start, t, t)
          x <- c(x, x[k], x[k]); id <- c(id, next_lineage, next_lineage + 1L)
          next_lineage <- next_lineage + 2L; next_node <- next_node + 1L
        } else {
          edges[[length(edges) + 1]] <- c(parent[k], id[k], t - start[k], 1)
        }
      }
      gone <- which(birth | death)
      if (length(gone)) {
        parent <- parent[-gone]; start <- start[-gone]; x <- x[-gone]; id <- id[-gone]
      }
    }
    if (length(id) < n) next
    t <- t + dt  # tips born in the last step get a branch of one step, not zero
    for (k in seq_along(id)) edges[[length(edges) + 1]] <- c(parent[k], id[k], t - start[k], 0)
    e <- do.call(rbind, edges)
    tips <- e[e[, 2] > 0, , drop = FALSE]
    n_tips <- nrow(tips)
    # node numbers: tips 1..n_tips, internal nodes n_tips + node label (root is node 1)
    tip_number <- stats::setNames(seq_len(n_tips), tips[, 2])
    child <- ifelse(e[, 2] > 0, tip_number[as.character(e[, 2])], n_tips - e[, 2])
    tree <- list(edge = cbind(n_tips + e[, 1], child), edge.length = e[, 3],
                 tip.label = paste0("t", tips[, 2]), Nnode = next_node - 1L)
    class(tree) <- "phylo"
    tree <- ape::reorder.phylo(tree)
    alive <- paste0("t", id)
    tree <- ape::keep.tip(tree, alive)
    return(list(tree = tree, x = stats::setNames(x, alive)[tree$tip.label]))
  }
  stop("no tree reached n tips")
}

# DR statistic (Jetz et al. 2012)
dr_stat <- function(tree) {
  n <- ape::Ntip(tree)
  len <- tree$edge.length[match(seq_len(n + tree$Nnode), tree$edge[, 2])]
  paths <- ape::nodepath(tree)
  vapply(paths, function(p) {
    l <- rev(len[p[-1]])                     # from the tip towards the root
    1 / sum(l / 2^(seq_along(l) - 1))
  }, numeric(1))
}

# two-state Markov estimates of colonisation and extirpation from presence over time
markov_rates <- function(history) {
  pres <- history > 0
  steps <- dim(pres)[3]
  t(vapply(seq_len(dim(pres)[2]), function(j) {
    a <- pres[, j, -steps, drop = FALSE]; b <- pres[, j, -1, drop = FALSE]
    n0 <- sum(!a); n1 <- sum(a)
    col <- if (n0 > 0) sum(!a & b) / n0 else NA_real_
    ext <- if (n1 > 0) sum(a & !b) / n1 else NA_real_
    c(col = col, ext = ext)
  }, numeric(2)))
}

one_replicate <- function(b, link = c("extirpation", "dispersal"), window = 30, n = 200,
                          strength = 1, n_steps = 600, base = 0.02) {
  link <- match.arg(link)
  clock <- list()
  tic <- Sys.time()
  sim_tree <- sim_tree_trait(n = n, b = b)
  tree <- sim_tree$tree
  z <- as.vector(scale(sim_tree$x))
  clock$tree <- Sys.time() - tic
  extirpation <- if (link == "extirpation") base * exp(-strength * z) else 0.01
  dispersal <- if (link == "dispersal") base * exp(strength * z) else 0.01
  tic <- Sys.time()
  sim <- simulate_metacommunity(landscape_grid(10), species_traits(tree),
                                competition_matrix(tree, "none"), dispersal = pmin(dispersal, 1),
                                extirpation = pmin(extirpation, 1), n_steps = n_steps,
                                record_steps = (n_steps - window + 1):n_steps)
  clock$metacommunity <- Sys.time() - tic
  tic <- Sys.time()
  rates <- markov_rates(sim$history)
  rownames(rates) <- dimnames(sim$history)[[2]]
  rates <- rates[tree$tip.label, , drop = FALSE]
  d <- data.frame(species = tree$tip.label, z = z, col = rates[, "col"], ext = rates[, "ext"],
                  dr = dr_stat(tree), lambda = exp(b * sim_tree$x))
  d$eps <- stats::qlogis(pmin(pmax(d$col / (d$col + d$ext), 1e-4), 1 - 1e-4))
  d$tau <- log(d$col + d$ext)
  ok <- stats::complete.cases(d) & is.finite(d$eps) & is.finite(d$tau)
  d <- d[ok, ]
  tr <- ape::keep.tip(tree, d$species)
  fit <- function(response) {
    if (stats::sd(d[[response]]) == 0) return(c(eps = NA, tau = NA, p_any = NA))
    f <- stats::as.formula(paste(response, "~ eps + tau"))
    m <- nlme::gls(f, data = d, correlation = ape::corBrownian(1, tr, form = ~species),
                   method = "ML")
    tt <- summary(m)$tTable
    c(eps = tt["eps", "p-value"], tau = tt["tau", "p-value"],
      p_any = stats::anova(m, update(m, . ~ 1))$`p-value`[2])
  }
  d$log_dr <- log(d$dr); d$log_lambda <- log(d$lambda)
  p_dr <- fit("log_dr"); p_true <- fit("log_lambda")
  clock$analysis <- Sys.time() - tic
  occupied <- apply(sim$history > 0, c(2, 3), mean)
  list(summary = data.frame(b = b, link = link, window = window, n_used = nrow(d),
                            n_alive_end = sum(colSums(sim$abundance) > 0),
                            mean_occupancy = mean(occupied),
                            cor_eps_z = stats::cor(d$eps, d$z, method = "spearman"),
                            cor_tau_z = stats::cor(d$tau, d$z, method = "spearman"),
                            cor_dr_lambda = stats::cor(d$dr, d$lambda, method = "spearman"),
                            p_dr_any = p_dr[["p_any"]], p_true_any = p_true[["p_any"]],
                            p_dr_eps = p_dr[["eps"]], p_dr_tau = p_dr[["tau"]],
                            p_true_eps = p_true[["eps"]], p_true_tau = p_true[["tau"]]),
       clock = vapply(clock, as.numeric, numeric(1), units = "secs"),
       data = d)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
runs <- rbind(
  expand.grid(b = c(0, 0.25, 0.5, 1), window = c(10, 30, 100), link = "dispersal",
              base = 0.01, rep = 1:200, stringsAsFactors = FALSE),
  expand.grid(b = c(0, 1), window = 30, link = "extirpation", base = 0.02, rep = 1:100,
              stringsAsFactors = FALSE))
start <- Sys.time()
out <- parallel::mclapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]
  tryCatch({
    x <- one_replicate(r$b, link = r$link, window = r$window, base = r$base)
    d <- x$data
    cbind(i = i, rep = r$rep, base = r$base, x$summary,
          secs = sum(x$clock),
          cor_col_z = stats::cor(d$col, d$z, method = "spearman"),
          cor_ext_z = stats::cor(d$ext, d$z, method = "spearman"))
  }, error = function(e) data.frame(i = i, error = conditionMessage(e)))
}, mc.cores = 6, mc.set.seed = TRUE)
ok <- vapply(out, function(o) is.null(o$error), logical(1))
results <- list(main = do.call(rbind, out[ok]), errors = do.call(rbind, out[!ok]), runs = runs,
                elapsed = difftime(Sys.time(), start, units = "mins"))
saveRDS(results, "data-raw/januario_power.rds")
cat("done in", format(results$elapsed), "with", sum(!ok), "errors\n")
