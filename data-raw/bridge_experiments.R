# Experiments for the bridge between Leibold et al. (2019, 2022) and Januario et al. (2025):
# in which metacommunities and sampling designs do the methods recover local adaptation (h),
# character displacement (c) and the imprint of partner turnover on branch rates?
#
# Run from the root of the package with treeshapesim and RRphylo installed:
#   Rscript data-raw/bridge_experiments.R
# Results in data-raw/bridge_experiments.rds (about 40 min on 6 cores).
#
# Each replicate:
# 1. a tree of 30 species, a 10 x 10 grid with environmental autocorrelation `bandwidth`, and a
#    metacommunity with `dispersal` and patch extinction (`extirpation`);
# 2. branch rates from the turnover of partners (turnover_rates()), traits evolved with them
#    (trait_rates()) and rates estimated with RRphylo: Spearman correlation with the true rate of
#    each branch and between the tip rate and the turnover;
# 3. individuals with local adaptation h and displacement c (sample_individuals(), the same
#    draws subsampled to 20 or 100 sites and 5 or 20 individuals per population), T-statistics
#    with the species null of T_IP.IC, and population_shifts() with the width free and fixed at
#    the true value;
# 4. on a subset, the parametric bootstrap of population_shifts().
#
# Caveat of scale: h and c are imposed when individuals are sampled, not evolved in the
# dynamics, and the branch rates come from a map of turnover to rate, not from emergent
# evolution; the two layers do not share the trait. The experiment asks whether the methods
# recover planted signals under realistic structures of assembly, not whether assembly causes
# trait evolution.

devtools::load_all(quiet = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
cores <- 6

spearman <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

# RRphylo rate of each branch, in the order of the rows of tree$edge
rrphylo_rates <- function(tree, tips) {
  rr <- suppressMessages(suppressWarnings(RRphylo::RRphylo(tree, tips, clus = 0)))
  n <- ape::Ntip(tree)
  node <- tree$edge[, 2]
  label <- ifelse(node <= n, tree$tip.label[pmin(node, n)], as.character(node))
  rr$rates[label, 1]
}

cells <- expand.grid(dispersal = c(0.001, 0.01, 0.1), extirpation = c(0, 0.01, 0.05),
                     bandwidth = c(0.5, 3), rep = 1:20)
processes <- expand.grid(h = c(0, 0.5), c = c(0, 1, 3))
designs <- expand.grid(n_sites = c(20, 100), n_ind = c(5, 20))
intraspecific_sd <- 5
population_sd <- 1

# population_shifts() with a record of its warning
shifts <- function(ind, ...) {
  warned <- FALSE
  fit <- tryCatch(withCallingHandlers(population_shifts(ind, ...), warning = function(w) {
    if (grepl("poorly identified", conditionMessage(w))) warned <<- TRUE
    invokeRestart("muffleWarning")
  }), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cbind(as.data.frame(fit)[setdiff(names(fit), c("n_populations", "n_species"))],
        n_populations = fit$n_populations, n_species = fit$n_species, warned = warned)
}

subsample <- function(ind, n_sites, n_ind) {
  sites <- unique(ind$site)
  keep <- if (n_sites < length(sites)) sample(sites, n_sites) else sites
  ind[ind$site %in% keep & ind$individual <= n_ind, ]
}

structure_of <- function(ind) {
  pops <- unique(ind[c("site", "species", "env")])
  presence <- table(pops$site, pops$species) > 0
  jaccard <- as.vector(stats::dist(presence, method = "binary"))
  data.frame(species_per_site = mean(rowSums(presence)),
             beta_jaccard = mean(jaccard, na.rm = TRUE),
             env_sd_within_species = mean(tapply(pops$env, pops$species, stats::sd),
                                          na.rm = TRUE))
}

one_replicate <- function(i) {
  cell <- cells[i, ]
  tree <- treeshapesim::sim_tree_shape(30)
  traits <- species_traits(tree)
  sim <- simulate_metacommunity(landscape_grid(10, env_bandwidth = cell$bandwidth), traits,
                                competition_matrix(tree, "phylo"), dispersal = cell$dispersal,
                                extirpation = cell$extirpation, n_steps = 600,
                                record_steps = seq(300, 600, 20))
  # macro layer
  macro <- tryCatch({
    rates <- turnover_rates(sim, k = 100, window = 0.1)
    x <- trait_rates(tree, edge_rates = rates)
    estimate <- rrphylo_rates(tree, x$tips)
    tip_edge <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
    data.frame(rho_branch = spearman(abs(estimate), x$edge_rates),
               rho_turnover = spearman(abs(estimate[tip_edge]),
                                       attr(rates, "turnover")[tree$tip.label]),
               rate_cv = stats::sd(rates[tip_edge]) / mean(rates[tip_edge]))
  }, error = function(e) data.frame(rho_branch = NA, rho_turnover = NA, rate_cv = NA))
  width <- mean(sim$traits$sigma)
  # individual layer
  rows <- list()
  for (p in seq_len(nrow(processes))) {
    full <- sample_individuals(sim, n = 20, local_adaptation = processes$h[p],
                               displacement = processes$c[p], population_sd = population_sd,
                               intraspecific_sd = intraspecific_sd)
    struct <- structure_of(full)
    for (d in seq_len(nrow(designs))) {
      ind <- subsample(full, designs$n_sites[d], designs$n_ind[d])
      ts <- tryCatch(suppressWarnings(tstats(ind, n_perm = 99, null_ip = "species")),
                     error = function(e) NULL)
      ses <- if (is.null(ts)) rep(NA_real_, 3) else
        colMeans(ts[c("ses_T_IP.IC", "ses_T_IC.IR", "ses_T_PC.PR")], na.rm = TRUE)
      free <- shifts(ind)
      fixed <- shifts(ind, width = width)
      rows[[length(rows) + 1]] <- data.frame(
        i = i, cell, h = processes$h[p], c = processes$c[p], designs[d, ], struct, macro,
        true_width = width,
        ses_T_IP.IC = ses[1], ses_T_IC.IR = ses[2], ses_T_PC.PR = ses[3],
        h_free = free$local_adaptation %||% NA, c_free = free$displacement %||% NA,
        w_free = free$width %||% NA, r2_free = free$r_squared %||% NA,
        warned_free = free$warned %||% NA, n_pops_used = free$n_populations %||% NA,
        h_fixed = fixed$local_adaptation %||% NA, c_fixed = fixed$displacement %||% NA,
        row.names = NULL)
    }
  }
  # bootstrap subset
  boot <- NULL
  if (cell$rep <= 5) {
    for (cc in c(0, 3)) {
      full <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = cc,
                                 population_sd = population_sd,
                                 intraspecific_sd = intraspecific_sd)
      for (d in c(1, 4)) {
        ind <- subsample(full, designs$n_sites[d], designs$n_ind[d])
        b <- shifts(ind, nboot = 99)
        if (!is.null(b)) {
          boot <- rbind(boot, data.frame(i = i, cell, h = 0.5, c = cc, designs[d, ], b,
                                         row.names = NULL))
        }
      }
    }
  }
  list(main = do.call(rbind, rows), boot = boot)
}

start <- Sys.time()
out <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
  tryCatch(one_replicate(i), error = function(e) list(error = conditionMessage(e), i = i))
}, mc.cores = cores, mc.set.seed = TRUE)
errors <- Filter(function(o) !is.null(o$error), out)
results <- list(
  main = do.call(rbind, lapply(out, `[[`, "main")),
  boot = do.call(rbind, lapply(out, `[[`, "boot")),
  errors = errors,
  settings = list(intraspecific_sd = intraspecific_sd, population_sd = population_sd,
                  n_species = 30, grid = 10, n_steps = 600, n_perm = 99, nboot = 99,
                  elapsed = difftime(Sys.time(), start, units = "mins"))
)
saveRDS(results, "data-raw/bridge_experiments.rds")
cat("done in", format(results$settings$elapsed), "with", length(errors), "errors\n")
