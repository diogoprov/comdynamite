# Second check of the tests of trait_partition(), with a generator not used to build them:
# counts from Gaussian niches on a hidden field of exponential covariance, two env columns,
# a grid of 10 by 10 sites, 8 Moran's eigenvector maps, 200 data sets, 199 permutations.
# Run from the root of the package, with the package installed. About 10 min on 5 cores.
library(comdynamite)
xy <- as.matrix(expand.grid(1:10, 1:10)); n <- 100; D <- as.matrix(dist(xy))
U <- chol(exp(-D / 3) + diag(1e-8, n))
fld <- function() as.vector(scale(drop(crossprod(U, rnorm(n)))))
mem <- suppressWarnings(as.matrix(adespatial::dbmem(xy, MEM.autocor = "positive", silent = TRUE)))[, 1:8]
rownames(xy) <- paste0("s", 1:n); rownames(mem) <- rownames(xy)
gen <- function(kind) {
  e1 <- if (kind %in% c("space_only_iid")) rnorm(n) else fld()
  e2 <- if (kind %in% c("space_only_iid")) rnorm(n) else fld()
  hid <- fld()
  driver <- switch(kind, null_iid = rnorm(n), space_only_iid = hid, space_only_auto = hid, env_auto = e1, both = e1 + hid)
  tr <- setNames(rnorm(25), paste0("sp", 1:25))
  mu <- exp(1 - outer(as.vector(scale(driver)), tr, function(g, t) (g - t)^2 / 2))
  comm <- matrix(rpois(length(mu), 3 * mu), n, dimnames = list(rownames(xy), names(tr)))
  keep <- rowSums(comm) > 0; list(comm = comm, tr = tr, env = `rownames<-`(cbind(e1, e2), rownames(xy)), keep = keep)
}
one <- function(i, kind) { set.seed(9000 + i); g <- gen(kind); if (!all(g$keep)) return(NULL)
  g$comm <- g$comm[, colSums(g$comm) > 0]
  o <- trait_partition(g$comm, g$tr[colnames(g$comm)], g$env, mem, n_perm = 199, shared_null = "msr", coords = xy)
  f <- o$fractions; rownames(f) <- f$fraction; c(a = f["a", "p_site"], c = f["c", "p_site"], amax = f["a", "p_max"], cmax = f["c", "p_max"], b = f["b", "p_site"], bmax = f["b", "p_max"], tested = o$shared$tested) }
results <- list()
for (k in c("null_iid", "space_only_iid", "space_only_auto", "env_auto", "both")) {
  r <- parallel::mclapply(1:200, function(i) tryCatch(one(i, k), error = function(e) conditionMessage(e)), mc.cores = 5, mc.preschedule = FALSE)
  bad <- vapply(r, is.character, TRUE); if (any(bad)) print(table(unlist(r[bad])))
  m <- do.call(rbind, r[!bad]); rj <- function(x) round(mean(x <= .05, na.rm = TRUE), 3)
  results[[k]] <- colMeans(m <= 0.05, na.rm = TRUE)
  cat(k, "reps", nrow(m), "tested", sum(m[, "tested"]), " a", rj(m[, "a"]), " amax", rj(m[, "amax"]), " c", rj(m[, "c"]), " cmax", rj(m[, "cmax"]), " b", rj(m[, "b"]), " bmax", rj(m[, "bmax"]), "\n")
}
saveRDS(do.call(rbind, results), "data-raw/validation_fractions_check.rds")
