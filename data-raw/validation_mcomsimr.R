# Compares comdynamite with mcomsimr, the original implementation of the model
# of Thompson et al. (2020). Both receive the same environment, traits,
# competition and dispersal, and we compare distributions across replicates,
# because the two do not share a random stream. Run from the package root;
# writes data-raw/validation_mcomsimr.rds. Needs mcomsimr
# (remotes::install_github("plthompson/mcomsimr")).
#
# Conventions: sigma here is sqrt(2) times env_niche_breadth in mcomsimr, and
# alpha[j, k] is the effect of k on j, so int_mat = t(alpha). In mcomsimr
# disp_mat[to, from]; in comdynamite dispersal[from, to]. The environment is
# held constant, and an implicit landscape is a uniform disp_mat without a
# diagonal. On a grid, mcomsimr distributes emigrants multinomially and
# comdynamite draws immigrants from a Poisson distribution, so the number of
# dispersers is conserved only in expectation.
library(comdynamite)
library(mcomsimr)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args)) as.integer(args[1]) else 200
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)

side <- 5
n_sites <- side^2
n_species <- 12
max_r <- 5
breadth <- 0.15 # env_niche_breadth of mcomsimr
n_init <- 100
n_after <- 200

scenarios <- expand.grid(dispersal = c(0.001, 0.01, 0.1, 0.5),
                         competition = c("intraspecific", "stabilising"),
                         landscape = c("implicit", "grid"),
                         stringsAsFactors = FALSE)

setup <- function(s) {
  env <- runif(n_sites)
  optima <- seq(0, 1, length.out = n_species)
  alpha <- matrix(if (s$competition == "stabilising") runif(n_species^2) else 0,
                  n_species, n_species)
  diag(alpha) <- 1
  alpha <- alpha * 0.05
  grid <- s$landscape == "grid"
  xy <- expand.grid(x = seq_len(side), y = seq_len(side))
  from_to <- if (grid) {
    dispersal_kernel(as.matrix(dist(xy)), scale = 1.5)
  } else {
    (1 - diag(n_sites)) / (n_sites - 1)
  }
  traits <- data.frame(optimum = optima, sigma = sqrt(2) * breadth, r_max = max_r,
                       row.names = paste0("sp", seq_len(n_species)))
  list(env = env, optima = optima, alpha = alpha, xy = xy, from_to = from_to,
       traits = traits, grid = grid, dispersal = s$dispersal,
       env_df = data.frame(env1 = rep(env, n_after), patch = rep(seq_len(n_sites), n_after),
                           time = rep(seq_len(n_after), each = n_sites)))
}

run_comdynamite <- function(p) {
  r <- growth_rates(p$env, p$traits)
  mc_dynamics(r, p$alpha, p$dispersal, n_steps = n_init + n_after,
              disp_matrix = if (p$grid) p$from_to)$abundance
}

run_mcomsimr <- function(p) {
  utils::capture.output(out <- suppressWarnings(suppressMessages(simulate_MC(
    patches = n_sites, species = n_species, dispersal = p$dispersal, plot = FALSE,
    timesteps = n_after, burn_in = 0, initialization = n_init, max_r = max_r,
    env_niche_breadth = breadth, landscape = p$xy, disp_mat = t(p$from_to),
    env.df = p$env_df, env_optima = p$optima, int_mat = t(p$alpha)))))
  d <- out$dynamics.df
  d <- d[d$time == max(d$time), ]
  n <- matrix(0L, n_sites, n_species)
  n[cbind(d$patch, d$species)] <- as.integer(d$N)
  n
}

statistics <- function(n) {
  pa <- (n > 0) * 1
  keep <- rowSums(pa) > 0
  c(alpha = mean(rowSums(pa)), gamma = sum(colSums(pa) > 0), total = sum(n),
    beta = if (sum(keep) > 1) mean(dist(pa[keep, ], method = "binary")) else NA,
    occupancy = mean(pa))
}

one_scenario <- function(i) {
  s <- scenarios[i, ]
  p <- setup(s)
  runs <- list(comdynamite = replicate(n_rep, run_comdynamite(p), simplify = FALSE),
               mcomsimr = replicate(n_rep, run_mcomsimr(p), simplify = FALSE))
  stats <- lapply(runs, function(x) t(vapply(x, statistics, numeric(5))))
  cells <- lapply(runs, function(x) Reduce(`+`, x) / n_rep)
  list(scenario = s, stats = stats, cells = cells)
}

results <- parallel::mclapply(seq_len(nrow(scenarios)), one_scenario,
                              mc.cores = nrow(scenarios))
saveRDS(results, "data-raw/validation_mcomsimr.rds", version = 2)

# Summary: difference in means and Kolmogorov-Smirnov test for each statistic,
# with false discovery rate control across all comparisons.
compare <- function(x) {
  a <- x$stats$comdynamite
  b <- x$stats$mcomsimr
  do.call(rbind, lapply(colnames(a), function(v) {
    data.frame(x$scenario, statistic = v, comdynamite = mean(a[, v], na.rm = TRUE),
               mcomsimr = mean(b[, v], na.rm = TRUE),
               p = suppressWarnings(ks.test(a[, v], b[, v])$p.value))
  }))
}
comparison <- do.call(rbind, lapply(results, compare))
comparison$relative <- comparison$comdynamite / comparison$mcomsimr - 1
comparison$q <- p.adjust(comparison$p, "BH")
print(comparison, digits = 3)
