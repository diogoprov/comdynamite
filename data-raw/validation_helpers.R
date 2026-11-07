# Shared code of validation_offgrid.R and validation_gap.R: landscapes, layouts of areas, the
# community of validation_edges.R, the p-values of node_test() from given permutations and
# the table of rejection rates. Sourced by those scripts; not run on its own.
library(comdynamite)

square <- as.matrix(expand.grid(1:20, 1:20))
grid <- as.matrix(expand.grid(1:70, 1:30))
band_middle <- 15 + 9 * sin(grid[, 1] / 9)
band_half <- 3 + 1.5 * cos(grid[, 1] / 5)
band <- unname(grid[abs(grid[, 2] - band_middle) <= band_half, ])
strip <- as.matrix(expand.grid(1:50, 1:8))

smooth_field <- function(xy, bandwidth) {
  d <- as.matrix(dist(xy))
  drop(exp(-d^2 / (2 * bandwidth^2)) %*% rnorm(nrow(xy)))
}

# sites within a box of coordinates, which also works off a lattice
box <- function(xy, xs, ys) {
  xy[, 1] >= min(xs) - 0.5 & xy[, 1] <= max(xs) + 0.5 & xy[, 2] >= min(ys) - 0.5 &
    xy[, 2] <= max(ys) + 0.5
}

# random layouts: two blocks of about 25 sites around sites drawn at random, irregular areas
# of several patches from the tails of a smooth field, one such area against all others
blocks <- function(xy) {
  repeat {
    centres <- xy[sample(nrow(xy), 2), ]
    one <- pmax(abs(xy[, 1] - centres[1, 1]), abs(xy[, 2] - centres[1, 2])) <= 2.5
    two <- pmax(abs(xy[, 1] - centres[2, 1]), abs(xy[, 2] - centres[2, 2])) <= 2.5
    if (!any(one & two) && sum(one) >= 15 && sum(two) >= 15) return(one * 1 - two * 1)
  }
}
tails <- function(xy) {
  f <- smooth_field(xy, 2.5)
  limits <- quantile(f, c(0.1, 0.9))
  (f >= limits[2]) * 1 - (f <= limits[1]) * 1
}
others <- function(xy) {
  f <- smooth_field(xy, 2.5)
  ifelse(f >= quantile(f, 0.9), 1, -1)
}

# fixed layouts of the square of 20 by 20
fixed <- list(
  adjacent = function(xy) box(xy, 6:10, 8:12) * 1 - box(xy, 11:15, 8:12) * 1,
  corners = function(xy) box(xy, 1:5, 1:5) * 1 - box(xy, 16:20, 16:20) * 1,
  one_vs_all = function(xy) ifelse(box(xy, 3:6, 3:6) | box(xy, 12:15, 13:16), 1, -1),
  ends = function(xy) box(xy, 1:6, 2:7) * 1 - box(xy, 45:50, 2:7) * 1
)

# sparse landscapes: a fraction of the cells of a lattice, or clusters of cells of the square
sparse_sites <- function(xy, keep) xy[sort(sample(nrow(xy), round(keep * nrow(xy)))), ]
cluster_sites <- function(n_clusters = 7, radius = 2.5) {
  repeat {
    centres <- square[sample(nrow(square), n_clusters), ]
    d <- as.matrix(dist(centres))
    if (min(d[upper.tri(d)]) >= 2.5 * radius) break
  }
  d <- sqrt(outer(square[, 1], centres[, 1], "-")^2 + outer(square[, 2], centres[, 2], "-")^2)
  cluster <- apply(d, 1, which.min)
  keep <- d[cbind(seq_len(nrow(square)), cluster)] <= radius
  list(xy = square[keep, ], cluster = cluster[keep])
}
# the two clusters farthest apart of a clustered landscape
far_clusters <- function(xy, cluster) {
  centres <- t(sapply(split(seq_len(nrow(xy)), cluster), function(i) colMeans(xy[i, , drop = FALSE])))
  d <- as.matrix(dist(centres))
  pair <- which(d == max(d), arr.ind = TRUE)[1, ]
  ifelse(cluster == rownames(centres)[pair[1]], 1, ifelse(cluster == rownames(centres)[pair[2]], -1, 0))
}

# p_partial_site of the event for each node, as node_test() computes it, from given
# permutations of sites (`perms`) or free permutations within the sites of the event
event_p <- function(comm, nodes, env, event, perms = NULL, n_perm = 199) {
  sites <- event != 0
  local_comm <- comm[sites, , drop = FALSE]
  predictors <- cbind(event[sites], env[sites, , drop = FALSE])
  draws <- if (is.null(perms)) {
    lapply(seq_len(n_perm), function(r) {
      list(comm = local_comm[sample.int(nrow(local_comm)), , drop = FALSE], env = predictors)
    })
  } else {
    lapply(seq_len(nrow(perms)), function(r) {
      moved <- event[perms[r, ]]
      s <- moved != 0
      list(comm = comm[s, , drop = FALSE], env = cbind(moved[s], env[s, , drop = FALSE]))
    })
  }
  vapply(seq_len(ncol(nodes)), function(i) {
    contrast <- nodes[, i]
    members <- which(contrast != 0)
    values <- contrast[members]
    observed <- comdynamite:::node_fit(local_comm[, members, drop = FALSE], values, predictors)
    if (is.na(observed$statistic)) return(NA_real_)
    null <- vapply(draws, function(d) {
      comdynamite:::node_fit(d$comm[, members, drop = FALSE], values, d$env)$partial_f[1]
    }, numeric(1))
    comdynamite:::permutation_p(null, observed$partial_f[1])
  }, numeric(1))
}

# one replicate: the community of validation_edges.R on landscape `xy` with event `event`
simulate_rep <- function(xy, event, effect, bandwidth) {
  tree <- treeshapesim::sim_tree_shape(12)
  hidden <- smooth_field(xy, bandwidth)
  hidden <- (hidden - min(hidden)) / diff(range(hidden)) * 100
  traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
  comm <- simulate_static(traits, hidden + effect * event, incidence = TRUE)$abundance
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  if (ncol(comm) < 4) return(NULL)
  nodes <- node_contrasts(ape::keep.tip(tree, colnames(comm)))
  list(comm = comm, nodes = nodes, noise = cbind(noise = rnorm(nrow(xy))))
}

# the permutations of one method, or NULL when the method stops. `spec` has `method`
# ("patch" or "rigid") and the arguments of the function
method_perms <- function(spec, xy, events, n_perm) {
  tryCatch(suppressWarnings(
    if (spec$method == "patch") {
      do.call(patch_permutations, c(list(xy, events, n_perm), spec[setdiff(names(spec), "method")]))$event
    } else {
      do.call(rigid_permutations, c(list(xy, events, n_perm), spec[setdiff(names(spec), "method")]))$event
    }), error = function(e) NULL)
}

# one replicate of a scenario: a data frame with one row per node and the p-value of each
# method, with the number of its permutations in `n_<method>`
one_rep <- function(i, scenario, data, effect, bandwidth, methods, n_perm, free = FALSE) {
  xy <- data$xy
  event <- data$event
  if (min(sum(event == 1), sum(event == -1)) < 5) return(NULL)
  sim <- simulate_rep(xy, event, effect, bandwidth)
  if (is.null(sim)) return(NULL)
  events <- cbind(event = event)
  out <- data.frame(rep = i, scenario = scenario, effect = effect, bandwidth = bandwidth,
                    node = seq_len(ncol(sim$nodes)))
  if (free) out$free <- event_p(sim$comm, sim$nodes, sim$noise, event, NULL, n_perm)
  for (m in names(methods)) {
    t0 <- Sys.time()
    perms <- method_perms(methods[[m]], xy, events, n_perm)
    out[[paste0("t_", m)]] <- as.numeric(Sys.time() - t0, units = "secs")
    out[[m]] <- if (is.null(perms)) NA_real_ else event_p(sim$comm, sim$nodes, sim$noise, event, perms)
    out[[paste0("n_", m)]] <- if (is.null(perms)) 0L else nrow(perms)
  }
  out
}

# all replicates of a design, in parallel. `design` has one row per scenario with
# `scenario`, `effect` and `bandwidth`, and `scenario_data(scenario)` gives `xy` and `event`
run_design <- function(design, scenario_data, methods, n_rep, n_perm, cores, free = FALSE) {
  rows <- design[rep(seq_len(nrow(design)), each = n_rep), ]
  rows$i <- rep(seq_len(n_rep), nrow(design))
  t0 <- Sys.time()
  results <- parallel::mclapply(seq_len(nrow(rows)), function(k) {
    r <- rows[k, ]
    set.seed(r$i)
    tryCatch(one_rep(r$i, r$scenario, scenario_data(r$scenario), r$effect, r$bandwidth,
                     methods, n_perm, free),
             error = function(e) data.frame(err = conditionMessage(e), scenario = r$scenario))
  }, mc.cores = cores, mc.preschedule = FALSE)
  cat("time:", format(Sys.time() - t0), "\n")
  errors <- Filter(function(x) !is.null(x) && "err" %in% names(x), results)
  if (length(errors)) print(table(unlist(lapply(errors, function(e) paste(e$scenario, e$err)))))
  do.call(rbind, Filter(function(x) !is.null(x) && !"err" %in% names(x), results))
}

# rejection rate at 5% of each method per scenario, with a standard error at the level of
# the replicate (nodes of one replicate are correlated), the number of replicates where the
# method ran, the number of replicates with at least 19 permutations (the smallest number
# that can reach the 5% level), the mean number of permutations and the mean time of the
# permutations in seconds
rate_table <- function(d, methods, by = c("scenario", "effect", "bandwidth")) {
  key <- interaction(d[by], drop = TRUE)
  rows <- lapply(split(d, key), function(s) {
    per <- lapply(methods, function(m) {
      rep_rate <- tapply(s[[m]] <= 0.05, s$rep, mean, na.rm = TRUE)
      rep_rate <- rep_rate[is.finite(rep_rate)]
      per_rep <- function(v) if (is.null(v)) NA else tapply(v, s$rep, `[`, 1)
      n_perm <- per_rep(s[[paste0("n_", m)]])
      data.frame(method = m, rate = 100 * mean(rep_rate),
                 se = 100 * sd(rep_rate) / sqrt(length(rep_rate)), n_rep = length(rep_rate),
                 n_testable = if (all(is.na(n_perm))) NA else sum(n_perm >= 19),
                 n_perm = if (all(is.na(n_perm))) NA else round(mean(n_perm)),
                 seconds = if (is.null(s[[paste0("t_", m)]])) NA else
                   round(mean(per_rep(s[[paste0("t_", m)]])), 2))
    })
    data.frame(s[1, by], do.call(rbind, per), row.names = NULL)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
