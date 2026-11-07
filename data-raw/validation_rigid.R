# Evaluation of rigid_permutations() against patch_permutations() and free permutations.
# Part 1: layouts and landscapes with a kernel field (150 replicates, 199 permutations); the
# layouts used to choose `outside` and `independent` were adjacent, corners, halves, random
# blocks and tails in the square and corners with 25% of the sites. The others were held out.
# Part 2: coordinates off a lattice (60 replicates, 99 permutations).
# Part 3: node_test() itself, with a field of exponential covariance, a generator not used to
# build the method (150 replicates, 199 permutations).
# About 1.5 h on 12 cores. Run from the root of the package, with the package installed.
library(comdynamite)
n_rep <- as.integer(Sys.getenv("N_REP", "150"))
n_rep_jitter <- as.integer(Sys.getenv("N_REP_JITTER", "60"))
cores <- as.integer(Sys.getenv("CORES", "12"))
out_file <- Sys.getenv("OUT", "data-raw/validation_rigid.rds")

grid <- as.matrix(expand.grid(1:70, 1:30))
middle <- 15 + 9 * sin(grid[, 1] / 9)
half <- 3 + 1.5 * cos(grid[, 1] / 5)
square <- as.matrix(expand.grid(1:20, 1:20))
landscapes <- list(square = square,
                   strip = as.matrix(expand.grid(1:50, 1:8)),
                   band = unname(grid[abs(grid[, 2] - middle) <= half, ]))

smooth_field <- function(xy, bandwidth) {
  d <- as.matrix(dist(xy))
  drop(exp(-d^2 / (2 * bandwidth^2)) %*% rnorm(nrow(xy)))
}

# random layouts of validation_edges.R
blocks <- function(xy) {
  repeat {
    centres <- xy[sample(nrow(xy), 2), ]
    one <- pmax(abs(xy[, 1] - centres[1, 1]), abs(xy[, 2] - centres[1, 2])) <= 2
    two <- pmax(abs(xy[, 1] - centres[2, 1]), abs(xy[, 2] - centres[2, 2])) <= 2
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
# fixed layouts of the square (also evaluated on subsets of its sites)
box <- function(xy, xs, ys) xy[, 1] %in% xs & xy[, 2] %in% ys
fixed <- list(
  adjacent = function(xy) box(xy, 6:10, 8:12) * 1 - box(xy, 11:15, 8:12) * 1,
  corners = function(xy) box(xy, 1:5, 1:5) * 1 - box(xy, 16:20, 16:20) * 1,
  halves = function(xy) ifelse(xy[, 1] <= 10, 1, -1),
  strips = function(xy) ifelse(xy[, 1] <= 2, 1, ifelse(xy[, 1] >= 19, -1, 0)),
  multi = function(xy) {
    (box(xy, 2:4, 2:4) | box(xy, 9:11, 15:17) | box(xy, 16:18, 5:7)) * 1 -
      (box(xy, 2:4, 14:16) | box(xy, 10:12, 3:5) | box(xy, 17:19, 17:19)) * 1
  },
  one_vs_all = function(xy) ifelse(box(xy, 3:6, 3:6) | box(xy, 12:15, 13:16), 1, -1)
)
# the two clusters farthest apart of a clustered landscape
far_clusters <- function(xy, cluster) {
  centres <- t(sapply(split(seq_len(nrow(xy)), cluster), function(i) colMeans(xy[i, , drop = FALSE])))
  d <- as.matrix(dist(centres))
  pair <- which(d == max(d), arr.ind = TRUE)[1, ]
  ifelse(cluster == rownames(centres)[pair[1]], 1, ifelse(cluster == rownames(centres)[pair[2]], -1, 0))
}
# sparse landscapes: a fraction of the cells of the square, or clusters of cells
sparse_sites <- function(keep) square[sort(sample(nrow(square), round(keep * nrow(square)))), ]
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
simulate_rep <- function(i, xy, event, effect, bandwidth) {
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

# rejection rate at 5% of each method per scenario, with a standard error at the level of
# the replicate (nodes of one replicate are correlated), the number of replicates where the
# method could run and the mean number of distinct placements
rate_table <- function(d, methods, by = c("scenario", "effect", "bandwidth")) {
  key <- interaction(d[by], drop = TRUE)
  rows <- lapply(split(d, key), function(s) {
    per <- lapply(c("free", methods), function(m) {
      rep_rate <- tapply(s[[m]] <= 0.05, s$rep, mean, na.rm = TRUE)
      rep_rate <- rep_rate[is.finite(rep_rate)]
      n_distinct <- if (m == "free") NA else mean(tapply(s[[paste0("n_", m)]], s$rep, `[`, 1))
      data.frame(method = m, rate = 100 * mean(rep_rate),
                 se = 100 * sd(rep_rate) / sqrt(length(rep_rate)), n_rep = length(rep_rate),
                 n_distinct = round(n_distinct))
    })
    data.frame(s[1, by], do.call(rbind, per), row.names = NULL)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# landscape and event of a replicate: fixed layouts of the square, random layouts of any
# landscape, sparse and clustered sites
scenario_data <- function(land, layout) {
  if (land %in% names(landscapes)) {
    xy <- landscapes[[land]]
  } else if (land == "sparse40") {
    xy <- sparse_sites(0.4)
  } else if (land == "sparse25") {
    xy <- sparse_sites(0.25)
  } else if (land == "clusters") {
    cl <- cluster_sites()
    xy <- cl$xy
  }
  event <- if (layout %in% names(fixed)) {
    fixed[[layout]](xy)
  } else if (layout == "far_clusters") {
    far_clusters(xy, cl$cluster)
  } else {
    get(layout)(xy)
  }
  list(xy = xy, event = event)
}

one_rep <- function(i, scenario, land, layout, effect, bandwidth, methods, n_perm) {
  set.seed(i)
  data <- scenario_data(land, layout)
  xy <- data$xy
  event <- data$event
  if (min(sum(event == 1), sum(event == -1)) < 5) return(NULL)
  sim <- simulate_rep(i, xy, event, effect, bandwidth)
  if (is.null(sim)) return(NULL)
  events <- cbind(event = event)
  out <- data.frame(rep = i, scenario = scenario, land = land, layout = layout, effect = effect,
                    bandwidth = bandwidth, node = seq_len(ncol(sim$nodes)))
  out$free <- event_p(sim$comm, sim$nodes, sim$noise, event, NULL, n_perm)
  for (m in names(methods)) {
    spec <- methods[[m]]
    perms <- tryCatch(suppressWarnings(
      if (spec$method == "patch") {
        patch_permutations(xy, events, n_perm, gap = spec$gap)$event
      } else {
        rigid_permutations(xy, events, n_perm, outside = spec$outside,
                           independent = spec$independent)$event
      }), error = function(e) NULL)
    out[[m]] <- if (is.null(perms)) NA_real_ else event_p(sim$comm, sim$nodes, sim$noise, event, perms)
    out[[paste0("n_", m)]] <- if (is.null(perms)) 0L else nrow(perms)
  }
  out
}

run_design <- function(design, methods, n_rep, n_perm, cores) {
  rows <- design[rep(seq_len(nrow(design)), each = n_rep), ]
  rows$i <- rep(seq_len(n_rep), nrow(design))
  t0 <- Sys.time()
  results <- parallel::mclapply(seq_len(nrow(rows)), function(k) {
    r <- rows[k, ]
    tryCatch(one_rep(r$i, r$scenario, r$land, r$layout, r$effect, r$bandwidth, methods, n_perm),
             error = function(e) data.frame(err = conditionMessage(e), scenario = r$scenario))
  }, mc.cores = cores, mc.preschedule = FALSE)
  cat("time:", format(Sys.time() - t0), "\n")
  errors <- Filter(function(x) !is.null(x) && "err" %in% names(x), results)
  if (length(errors)) print(table(unlist(lapply(errors, function(e) paste(e$scenario, e$err)))))
  out <- do.call(rbind, Filter(function(x) !is.null(x) && !"err" %in% names(x), results))
  out
}


# part 1
n_perm <- 199
methods <- list(
  patch = list(method = "patch", gap = NULL),
  rigid = list(method = "rigid", outside = 0.1, independent = Inf),
  rigid_i10 = list(method = "rigid", outside = 0.1, independent = 10),
  rigid_i6 = list(method = "rigid", outside = 0.1, independent = 6)
)
sq <- function(scenario, layout, effect = 0, bandwidth = 3, land = "square") {
  data.frame(scenario = scenario, land = land, layout = layout, effect = effect,
             bandwidth = bandwidth)
}
design <- rbind(
  # fixed layouts of the square, two bandwidths
  sq(names(fixed), names(fixed)),
  sq(paste0(names(fixed), "_bw12"), names(fixed), bandwidth = 12),
  # random layouts of the three landscapes
  sq(paste0(rep(c("blocks", "tails", "others"), 3), "_", rep(c("square", "strip", "band"), each = 3)),
     rep(c("blocks", "tails", "others"), 3), land = rep(c("square", "strip", "band"), each = 3)),
  sq(paste0(c("blocks", "tails", "others"), "_square_bw12"), c("blocks", "tails", "others"),
     bandwidth = 12),
  sq(paste0(c("blocks", "tails"), "_band_bw12"), c("blocks", "tails"), bandwidth = 12,
     land = "band"),
  # sparse sites
  sq("corners_sparse25", "corners", land = "sparse25"),
  sq("corners_sparse40", "corners", land = "sparse40"),
  sq("tails_sparse40", "tails", land = "sparse40"),
  sq("far_clusters", "far_clusters", land = "clusters"),
  sq("tails_clusters", "tails", land = "clusters"),
  # power
  sq(paste0(c("adjacent", "corners", "multi", "one_vs_all"), "_e20"),
     c("adjacent", "corners", "multi", "one_vs_all"), effect = 20),
  sq(paste0(c("blocks", "tails", "others"), "_square_e20"), c("blocks", "tails", "others"),
     effect = 20),
  sq(paste0(c("blocks", "tails"), "_band_e20"), c("blocks", "tails"), effect = 20,
     land = "band"),
  sq("far_clusters_e20", "far_clusters", effect = 20, land = "clusters")
)
cat(nrow(design), "scenarios\n")
final_rows <- run_design(design, methods, n_rep, n_perm, cores)
final <- rate_table(final_rows, names(methods))

# part 2
jitter_rep <- function(i, layout, effect) {
  set.seed(i)
  xy <- square + matrix(runif(800, -0.3, 0.3), ncol = 2)
  event <- if (layout %in% names(fixed)) fixed[[layout]](square) else get(layout)(xy)
  sim <- simulate_rep(i, xy, event, effect, 3)
  if (is.null(sim)) return(NULL)
  events <- cbind(event = event)
  out <- data.frame(rep = i, scenario = layout, effect = effect, bandwidth = 3,
                    node = seq_len(ncol(sim$nodes)))
  out$free <- event_p(sim$comm, sim$nodes, sim$noise, event, NULL, 99)
  draws <- list(
    patch = patch_permutations(xy, events, 99)$event,
    rigid = tryCatch(suppressWarnings(rigid_permutations(xy, events, 99))$event,
                     error = function(e) NULL),
    rigid_i10 = tryCatch(suppressWarnings(rigid_permutations(xy, events, 99,
                                                                independent = 10))$event,
                         error = function(e) NULL))
  for (m in names(draws)) {
    out[[m]] <- if (is.null(draws[[m]])) NA_real_ else
      event_p(sim$comm, sim$nodes, sim$noise, event, draws[[m]])
    out[[paste0("n_", m)]] <- if (is.null(draws[[m]])) 0L else nrow(draws[[m]])
  }
  out
}
jitter_design <- expand.grid(i = seq_len(n_rep_jitter), layout = c("adjacent", "corners", "blocks", "tails"),
                      effect = 0, stringsAsFactors = FALSE)
jitter_rows <- parallel::mclapply(seq_len(nrow(jitter_design)), function(k) {
  jitter_rep(jitter_design$i[k], jitter_design$layout[k], jitter_design$effect[k])
}, mc.cores = cores, mc.preschedule = FALSE)
jitter <- rate_table(do.call(rbind, jitter_rows), c("patch", "rigid", "rigid_i10"))

# part 3
xy <- as.matrix(expand.grid(1:20, 1:20)); n <- nrow(xy); D <- as.matrix(dist(xy))
strip <- as.matrix(expand.grid(1:50, 1:8)); Ds <- as.matrix(dist(strip))
box <- function(xy, xs, ys) xy[,1] %in% xs & xy[,2] %in% ys
ch <- list(sq4 = chol(exp(-D/4) + diag(1e-8, n)), sq8 = chol(exp(-D/8) + diag(1e-8, n)), st4 = chol(exp(-Ds/4) + diag(1e-8, nrow(strip))))
lay <- list(
 adjacent = list(xy, "sq", box(xy,6:10,8:12)*1 - box(xy,11:15,8:12)*1),
 lshape   = list(xy, "sq", (box(xy,3:8,3:5)|box(xy,3:5,3:9))*1 - box(xy,12:16,12:15)*1),
 midgap   = list(xy, "sq", box(xy,2:6,8:12)*1 - box(xy,13:17,8:12)*1),
 corners  = list(xy, "sq", box(xy,1:5,1:5)*1 - box(xy,16:20,16:20)*1),
 strips   = list(xy, "sq", ifelse(xy[,1]<=2,1,ifelse(xy[,1]>=19,-1,0))),
 strip_ends = list(strip, "st", box(strip,1:6,2:7)*1 - box(strip,45:50,2:7)*1))
design <- rbind(
 data.frame(l=c("adjacent","lshape","midgap"), r="4", ind=Inf),
 data.frame(l=c("adjacent","lshape","midgap"), r="8", ind=Inf),
 data.frame(l=c("corners","strips"), r="4", ind=10),
 data.frame(l=c("corners"), r="8", ind=10),
 data.frame(l=c("strip_ends"), r="4", ind=20))
one <- function(i, d) {
  set.seed(1000 + i); L <- lay[[d$l]]; x <- L[[1]]; ev <- L[[3]]; m <- nrow(x)
  U <- if (L[[2]]=="st") ch$st4 else ch[[paste0("sq",d$r)]]
  h <- drop(crossprod(U, rnorm(m))); h <- (h-min(h))/diff(range(h))*100
  tree <- treeshapesim::sim_tree_shape(12)
  tr <- species_traits(tree, lambda=1, optima_range=c(0,100), rmax_range=c(0,30))
  comm <- simulate_static(tr, h, incidence=TRUE)$abundance; comm <- comm[,colSums(comm)>0,drop=FALSE]
  if (ncol(comm)<4) return(NULL)
  rownames(comm) <- paste0("s",seq_len(m)); nodes <- node_contrasts(ape::keep.tip(tree,colnames(comm)))
  b <- cbind(event=ev); rownames(b) <- rownames(comm); env <- cbind(noise=rnorm(m)); rownames(env) <- rownames(comm)
  pn <- suppressWarnings(rigid_permutations(x,b,n_perm,independent=d$ind))
  po <- patch_permutations(x,b,n_perm)
  f <- function(p) { r <- node_test(comm,nodes,env,b,n_perm=n_perm,area_permutations=p)$biogeography$p_partial_site; mean(r[,"event"]<=.05,na.rm=TRUE) }
  c(new=f(pn), old=f(po), layouts=nrow(pn$event))
}
res <- lapply(seq_len(nrow(design)), function(k) { d <- design[k,]
  r <- parallel::mclapply(seq_len(n_rep), function(i) tryCatch(one(i,d),error=function(e) conditionMessage(e)), mc.cores=cores, mc.preschedule=FALSE)
  bad <- vapply(r,is.character,TRUE); if(any(bad)) print(table(unlist(r[bad])))
  m <- do.call(rbind,r[!bad & !vapply(r,is.null,TRUE)])
  o <- data.frame(d, reps=nrow(m), layouts=mean(m[,"layouts"]), new=100*mean(m[,"new"]), se_new=100*sd(m[,"new"])/sqrt(nrow(m)), old=100*mean(m[,"old"]), se_old=100*sd(m[,"old"])/sqrt(nrow(m)))
  print(o, digits=3, row.names=FALSE); flush.console(); o })
exponential <- do.call(rbind, res)

saveRDS(list(final = final, jitter = jitter, exponential = exponential), out_file)
