# Block A, item 1: null models for a continuous, spatially autocorrelated predictor in
# landscapes with edges. The community follows a hidden environment h (Gaussian field with
# exponential covariance of range R); the measured predictor x is rho * h + sqrt(1 - rho^2) * z,
# with z an independent field of the same covariance. With rho = 0 there is no relation, but
# both are autocorrelated, so free permutations are liberal. We compare the rejection rate of
# the site test of pcps_test() (whole composition, p_site) under several nulls for x:
# - free: permutations of the sites;
# - torus: torus_permutations() (square grid only, which has edges here);
# - msr: Moran spectral randomization (Wagner & Dray 2015), as in trait_partition();
# - shift_rotate: an ad hoc rotation and shift of the coordinates about the centroid,
#   wrapped within the bounding box, each site taking the value of the nearest site. It is
#   inspired by the shift-&-rotate approach (Ridder et al. 2024, Methods Ecol Evol,
#   doi:10.1111/2041-210x.14443), whose full text we could not read, so it is not their
#   method;
# - param_rank: a parametric surrogate, a Gaussian field with the exponential range fitted to
#   x by maximum likelihood, given the observed values of x by rank (the alternative that a
#   referee suggested to Viladomat et al. 2014);
# - vgram: the variogram matching of Viladomat et al. (2014, Biometrics,
#   doi:10.1111/biom.12139): x permuted, smoothed over a proportion delta of its nearest
#   neighbours, with delta chosen so that the variogram of the smoothed field best matches
#   that of x up to distances of the 25% quantile of the pairs, then rescaled (and noise
#   added for the nugget) so that scale and intercept match too. We use binned variograms
#   rather than their kernel-smoothed one, and an exponential kernel over the neighbours.
# Since version 0.10.0 the parametric and variogram nulls come from env_surrogates(), so the
# rates are those of the package code.
# PILOT=1 runs a small design to time the nulls. Run from the package root, with the package
# installed. Writes data-raw/validation_continuous_nulls(_pilot).rds.
library(comdynamite)

pilot <- Sys.getenv("PILOT", "0") == "1"
cores <- as.integer(Sys.getenv("CORES", "14"))
n_rep <- as.integer(Sys.getenv("N_REP", if (pilot) "20" else "200"))
n_null <- as.integer(Sys.getenv("N_NULL", if (pilot) "99" else "199"))
scales <- c(2, 5, 10)
rhos <- c(0, 0.6)

grid20 <- as.matrix(expand.grid(1:20, 1:20))
lshape <- function(m = 280) {
  repeat {
    u <- cbind(stats::runif(3 * m, 0, 24), stats::runif(3 * m, 0, 24))
    u <- u[u[, 1] < 12 | u[, 2] < 12, ]
    if (nrow(u) >= m) return(u[seq_len(m), ])
  }
}
geometry <- function(name) switch(name,
  grid = grid20,
  lshape = lshape(),
  strip = as.matrix(expand.grid(1:50, 1:8)))
torus20 <- torus_permutations(20)

field_root <- function(xy, range) t(chol(exp(-as.matrix(stats::dist(xy)) / range) + diag(1e-8, nrow(xy))))

nulls <- function(name, x, xy, n) {
  m <- length(x)
  switch(name,
    vgram = lapply(env_surrogates(x, xy, n, "variogram"), drop),
    free = replicate(n, x[sample.int(m)], simplify = FALSE),
    torus = lapply(sample.int(nrow(torus20), n), function(i) x[torus20[i, ]]),
    msr = {
      basis <- NULL
      utils::capture.output(basis <- adespatial::dbmem(xy, MEM.autocor = "all", silent = TRUE))
      out <- adespatial::msr(matrix(x, ncol = 1), basis, nrepet = n, method = "singleton",
                             simplify = FALSE)
      lapply(out, function(e) as.numeric(as.matrix(e)))
    },
    shift_rotate = {
      lo <- apply(xy, 2, min)
      width <- apply(xy, 2, max) - lo + 1
      centre <- colMeans(xy)
      replicate(n, {
        a <- stats::runif(1, 0, 2 * pi)
        rot <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
        moved <- sweep(sweep(xy, 2, centre) %*% rot, 2, centre + stats::runif(2) * width, "+")
        moved <- sweep(sweep(sweep(moved, 2, lo), 2, width, "%%"), 2, lo, "+")
        nearest <- apply(moved, 1, function(p) which.min(colSums((t(xy) - p)^2)))
        x[nearest]
      }, simplify = FALSE)
    },
    param_rank = lapply(env_surrogates(x, xy, n, "parametric"), drop))
}
null_names <- function(geo) c("free", if (geo == "grid") "torus", "msr", "shift_rotate", "param_rank",
                                "vgram")

one <- function(i, geo, range, rho) {
  set.seed(9000 + i)
  for (try in 1:20) {
    xy <- geometry(geo)
    m <- nrow(xy)
    root <- field_root(xy, range)
    h <- drop(root %*% stats::rnorm(m))
    z <- drop(root %*% stats::rnorm(m))
    tree <- treeshapesim::sim_tree_shape(12)
    tr <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
    hs <- (h - min(h)) / diff(range(h)) * 100
    comm <- simulate_static(tr, hs, incidence = TRUE)$abundance
    comm <- comm[, colSums(comm) > 0, drop = FALSE]
    if (ncol(comm) >= 4 && all(rowSums(comm) > 0)) break
  }
  if (ncol(comm) < 4 || any(rowSums(comm) == 0)) return(NULL)
  x <- rho * h + sqrt(1 - rho^2) * z
  tree <- ape::keep.tip(tree, colnames(comm))
  weights <- comm / rowSums(comm)
  gower <- comdynamite:::pcps_gower(weights, comdynamite:::fuzzy_belonging(tree))
  scores <- comdynamite:::pcps_axes(gower)$scores[, 1:2, drop = FALSE]
  stat <- function(e) comdynamite:::pcps_fit(gower, scores, matrix(e, ncol = 1))$statistic[1]
  observed <- stat(x)
  out <- list()
  for (nm in null_names(geo)) {
    t0 <- proc.time()[["elapsed"]]
    draws <- tryCatch(nulls(nm, x, xy, n_null), error = function(e) NULL)
    p <- if (is.null(draws)) NA_real_ else {
      null <- vapply(draws, stat, numeric(1))
      (1 + sum(null >= observed - 1e-12)) / (1 + length(null))
    }
    out[[nm]] <- c(p = p, seconds = proc.time()[["elapsed"]] - t0)
  }
  data.frame(rep = i, geometry = geo, range = range, rho = rho, null = names(out),
             p = vapply(out, `[[`, numeric(1), "p"),
             seconds = vapply(out, `[[`, numeric(1), "seconds"), row.names = NULL)
}

design <- expand.grid(i = seq_len(n_rep), geo = c("grid", "lshape", "strip"), range = scales,
                      rho = rhos, stringsAsFactors = FALSE)
t_start <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(design)), function(k) {
  d <- design[k, ]
  tryCatch(one(d$i, d$geo, d$range, d$rho), error = function(e) NULL)
}, mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, res)
attr(results, "elapsed") <- difftime(Sys.time(), t_start, units = "mins")
saveRDS(results, if (pilot) "data-raw/validation_continuous_nulls_pilot.rds" else
          "data-raw/validation_continuous_nulls.rds", version = 2)

summary_table <- aggregate(cbind(reject = p <= 0.05) ~ geometry + range + rho + null, results, mean)
print(reshape(summary_table, idvar = c("geometry", "range", "rho"), timevar = "null",
              direction = "wide"), digits = 2)
print(aggregate(seconds ~ geometry + null, results, mean), digits = 2)
cat("elapsed:", format(attr(results, "elapsed")), "; failed runs:",
    nrow(design) - length(unique(paste(results$rep, results$geometry, results$range,
                                       results$rho))), "\n")
