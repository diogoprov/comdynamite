# Power of tests of a trait-speciation link under a design that follows the empirical pipeline
# of tip-rate studies (Januario et al. 2025, Ecol. Lett. and Am. Nat.):
#
#   Rscript data-raw/januario_design.R            (REPS=200 by default; REPS=5 for a pilot)
# Results in data-raw/januario_design.rds.
#
# Each replicate:
# 1. a large tree (N_TIPS living species, default 4000) grows by an event-based (Gillespie)
#    birth-death process in which a Brownian trait X sets the speciation rate,
#    lambda = lambda0 * base * exp(b X); `base` is 1, or, in the "shifts" scenario, a
#    clade-specific background rate that jumps by a lognormal factor at rate nu along lineages,
#    independently of X (a confounder of the kind in Fig. 7 of the Am. Nat. paper);
# 2. lambda_DR (Jetz et al. 2012) on the complete tree, as tip rates are computed on large
#    trees before the species with data are selected;
# 3. ~200 species are sampled in three ways: a complete clade of 150-250 species (a
#    "family-level" study), 200 species at random, and a phylogenetically overdispersed set
#    (one random species from each of the ~200 lineages alive at the time when 200 lineages
#    exist);
# 4. the trait is measured with error: w = sqrt(R) z + sqrt(1 - R) e, z the standardised X of
#    the sampled species, e independent noise, R the reliability (1, or 0.75, about the
#    squared correlation of tau with the latent trait in data-raw/januario_power.R);
# 5. tests of log lambda_DR against w: PGLS under Brownian motion (regression of independent
#    contrasts through the origin, the same test as nlme::gls with corBrownian) and ES-sim
#    (Harvey & Rabosky 2018; Pearson correlation against 999 Brownian traits on the tree).
#
# The effect b is in biological units: species 1 SD apart in X differ in speciation rate by a
# factor exp(b sd(X)).

devtools::load_all(quiet = TRUE)

N_TIPS <- as.integer(Sys.getenv("N_TIPS", 4000))
REPS <- as.integer(Sys.getenv("REPS", 200))
CORES <- as.integer(Sys.getenv("CORES", 6))

# event-based birth-death with trait-dependent speciation and optional background shifts
sim_big_tree <- function(n = N_TIPS, b = 0, lambda0 = 1, mu = 0.3, sigma = 0.4,
                         nu = 0, shift_sd = 0.5, max_tries = 50) {
  for (try in seq_len(max_tries)) {
    cap <- 6 * n
    edges <- matrix(NA_real_, cap, 3)  # parent node, child (negative: node label), length
    n_edges <- 0L
    parent <- c(1L, 1L); start <- c(0, 0); x <- c(0, 0); base <- c(1, 1); id <- c(1L, 2L)
    next_node <- 2L; next_lineage <- 3L; t <- 0
    while (length(id) > 0 && length(id) < n) {
      lam <- lambda0 * base * exp(b * x)
      total <- sum(lam) + (mu + nu) * length(id)
      dt <- stats::rexp(1, total)
      t <- t + dt
      x <- x + stats::rnorm(length(x), 0, sigma * sqrt(dt))
      u <- stats::runif(1, 0, total)
      if (u < sum(lam)) {                       # speciation
        k <- which.max(cumsum(lam) >= u)
        if (n_edges + 2L > cap) { edges <- rbind(edges, matrix(NA_real_, cap, 3)); cap <- 2L * cap }
        n_edges <- n_edges + 1L
        edges[n_edges, ] <- c(parent[k], -next_node, t - start[k])
        parent <- c(parent[-k], next_node, next_node); start <- c(start[-k], t, t)
        x <- c(x[-k], x[k], x[k]); base <- c(base[-k], base[k], base[k])
        id <- c(id[-k], next_lineage, next_lineage + 1L)
        next_lineage <- next_lineage + 2L; next_node <- next_node + 1L
      } else if (u < sum(lam) + mu * length(id)) {  # extinction
        k <- sample.int(length(id), 1)
        n_edges <- n_edges + 1L
        edges[n_edges, ] <- c(parent[k], id[k], t - start[k])
        parent <- parent[-k]; start <- start[-k]; x <- x[-k]; base <- base[-k]; id <- id[-k]
      } else {                                   # background shift
        k <- sample.int(length(id), 1)
        base[k] <- base[k] * exp(stats::rnorm(1, 0, shift_sd))
      }
    }
    if (length(id) < n) next
    t <- t + 1e-3
    e <- rbind(edges[seq_len(n_edges), , drop = FALSE], cbind(parent, id, t - start))
    tips <- e[e[, 2] > 0, , drop = FALSE]
    n_tips <- nrow(tips)
    tip_number <- stats::setNames(seq_len(n_tips), tips[, 2])
    child <- ifelse(e[, 2] > 0, tip_number[as.character(e[, 2])], n_tips - e[, 2])
    tree <- list(edge = cbind(n_tips + e[, 1], child), edge.length = e[, 3],
                 tip.label = paste0("t", tips[, 2]), Nnode = next_node - 1L)
    class(tree) <- "phylo"
    tree <- ape::reorder.phylo(tree)
    alive <- paste0("t", id)
    tree <- ape::keep.tip(tree, alive)
    return(list(tree = tree, x = stats::setNames(x, alive)[tree$tip.label],
                base = stats::setNames(base, alive)[tree$tip.label]))
  }
  stop("no tree reached n tips")
}

# the three sampling schemes; each returns tip labels
sample_clade <- function(tree, lo = 150, hi = 250) {
  n <- ape::Ntip(tree)
  size <- tabulate(unlist(lapply(ape::nodepath(tree), function(p) p[-length(p)])),
                   nbins = n + tree$Nnode)[-(1:n)]
  nodes <- n + which(size >= lo & size <= hi)
  if (!length(nodes)) nodes <- n + which.min(abs(size - (lo + hi) / 2))
  node <- nodes[sample.int(length(nodes), 1)]
  ape::extract.clade(tree, node)$tip.label
}
sample_random <- function(tree, k = 200) sample(tree$tip.label, k)
sample_dispersed <- function(tree, k = 200) {
  depth <- ape::node.depth.edgelength(tree)
  starts <- depth[tree$edge[, 1]]; ends <- depth[tree$edge[, 2]]
  times <- sort(unique(depth[-(1:ape::Ntip(tree))]))
  crossing <- vapply(times, function(s) sum(starts <= s & ends > s), numeric(1))
  cut <- times[which.max(crossing >= k)]
  paths <- ape::nodepath(tree)
  # ancestor of each tip on the lineages alive just after the cut
  anc <- vapply(paths, function(p) p[which.max(depth[p] > cut)], numeric(1))
  picks <- tapply(seq_along(anc), anc, function(i) i[sample.int(length(i), 1)])
  tree$tip.label[unlist(picks)]
}

es_sim <- function(y, w, tree, n_null = 999) {
  root <- t(chol(ape::vcv(tree)[names(w), names(w)]))
  null <- stats::cor(y, root %*% matrix(stats::rnorm(length(w) * n_null), length(w)))
  observed <- stats::cor(y, w)
  c(r = observed, p = (1 + sum(abs(null) >= abs(observed))) / (n_null + 1))
}
pgls_pic <- function(y, w, tree) {
  cy <- ape::pic(y, tree); cw <- ape::pic(w, tree)
  m <- stats::lm(cy ~ cw - 1)
  t_hc3 <- stats::coef(m)[[1]] / sqrt(sandwich::vcovHC(m, type = "HC3")[1, 1])
  c(r_pic = sum(cy * cw) / sqrt(sum(cy^2) * sum(cw^2)),
    p = summary(m)$coefficients["cw", "Pr(>|t|)"],
    p_hc3 = 2 * stats::pt(-abs(t_hc3), df = length(cw) - 1))
}
# PGLS with Pagel's lambda estimated by maximum likelihood (phylolm)
pgls_lambda <- function(y, w, tree) {
  d <- data.frame(y = y, w = w, row.names = tree$tip.label)
  fit <- tryCatch(phylolm::phylolm(y ~ w, d, tree, model = "lambda"), error = function(e) NULL)
  if (is.null(fit)) return(c(p = NA_real_, lambda = NA_real_))
  c(p = summary(fit)$coefficients["w", "p.value"], lambda = fit$optpar)
}

one_replicate <- function(b, shifts) {
  s <- sim_big_tree(b = b, nu = if (shifts) 0.02 else 0)
  tree <- s$tree
  dr <- dr_rates(tree)[tree$tip.label]
  sets <- list(clade = sample_clade(tree), random = sample_random(tree),
               dispersed = sample_dispersed(tree))
  out <- list()
  for (scheme in names(sets)) {
    tips <- sets[[scheme]]
    tr <- ape::keep.tip(tree, tips)
    y <- log(dr[tr$tip.label])
    z <- as.vector(scale(s$x[tr$tip.label]))
    for (R in c(1, 0.75)) {
      w <- stats::setNames(sqrt(R) * z + sqrt(1 - R) * stats::rnorm(length(z)), tr$tip.label)
      pg <- pgls_pic(y, w, tr); pl <- pgls_lambda(y, w, tr); es <- es_sim(y, w, tr)
      out[[length(out) + 1]] <- data.frame(
        b = b, shifts = shifts, scheme = scheme, reliability = R, n = length(tips),
        r_raw = stats::cor(y, w), r_pic = pg[["r_pic"]], p_pgls = pg[["p"]],
        p_hc3 = pg[["p_hc3"]], p_lambda = pl[["p"]], lambda = pl[["lambda"]],
        p_es = es[["p"]], sd_x = stats::sd(s$x[tr$tip.label]),
        cor_dr_true = suppressWarnings(stats::cor(y, log(s$base[tr$tip.label]) + b * s$x[tr$tip.label])))
    }
  }
  do.call(rbind, out)
}

if (sys.nframe() == 0L) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(2028)
  runs <- rbind(expand.grid(b = c(0, 0.25, 0.5, 1), shifts = FALSE, rep = seq_len(REPS)),
                expand.grid(b = c(0, 0.5), shifts = TRUE, rep = seq_len(REPS)))
  start <- Sys.time()
  out <- parallel::mclapply(seq_len(nrow(runs)), function(i) {
    tryCatch(cbind(rep = runs$rep[i], one_replicate(runs$b[i], runs$shifts[i])),
             error = function(e) data.frame(i = i, error = conditionMessage(e)))
  }, mc.cores = CORES, mc.set.seed = TRUE)
  ok <- vapply(out, function(o) is.null(o$error), logical(1))
  res <- list(main = do.call(rbind, out[ok]), errors = do.call(rbind, out[!ok]),
              n_tips = N_TIPS, elapsed = difftime(Sys.time(), start, units = "mins"))
  saveRDS(res, Sys.getenv("OUT", "data-raw/januario_design.rds"))
  m <- res$main
  print(aggregate(cbind(pgls = p_pgls < 0.05, hc3 = p_hc3 < 0.05, lam = p_lambda < 0.05, es = p_es <= 0.05, r_raw, r_pic, n) ~
                    reliability + scheme + shifts + b, m, mean, na.action = na.pass), digits = 2)
  cat("done in", format(res$elapsed), "with", sum(!ok), "errors\n")
}
