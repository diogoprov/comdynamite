# Replication of the pulse-disturbance metacommunity model of Hajian-Forooshani & Chase (2026,
# Am Nat 207: 838-850, doi:10.1086/739955), from the article and the code archived with it
# (Figshare, doi:10.6084/m9.figshare.30620603.v1, file 2025_10_09_MetaDist_simluation_model.R).
#
# Their model is that of Thompson et al. (2020) with a single pulse disturbance:
# * 100 patches on a 10 x 10 lattice wrapped as a torus; dispersal kernel exp(-0.1 d), normalised
#   over the other patches;
# * 30 species with optima drawn from U(0.5, 1), r_max = 5 and a niche breadth n shared by all
#   species; growth r = r_max exp(-((z - env) / (2 n))^2) in their code;
# * competition: intraspecific 0.5 x 0.05, interspecific U(0, 1) x 0.05, an asymmetric matrix;
# * environment of each patch drawn from U(0.5, 1), spatially uncorrelated and constant except
#   for the pulse; at step 500 a random set of 1 or 5% to 100% of patches is set to 0.1 (below
#   every optimum) for steps 500 to 600, after which the patches return to their values;
# * Poisson Beverton-Holt growth, binomial emigration with probability d, and redistribution of
#   all emigrants of a species, conserving them; Poisson(1) individuals of every species in
#   every patch at the start, no later seeding; 1000 steps;
# * alpha = mean patch richness over all patches, gamma = richness of the landscape, changes as
#   log ratios relative to step 499, and Delta beta = Delta gamma - Delta alpha (the log ratio of
#   Whittaker's beta, as in their Fig. 2 script).
#
# Translation used here:
# * landscape_grid(side = 10, disp_scale = 10, torus = TRUE), since dispersal_kernel() is
#   exp(-d / scale); the environment of the landscape is replaced by env_pulse();
# * sigma = sqrt(2) n, since growth_rates() uses exp(-(z - env)^2 / (2 sigma^2)) and their code
#   exp(-(z - env)^2 / (4 n^2));
# * their competition matrix passed as is (asymmetric, not from competition_matrix());
# * the pulse with env_pulse(start = 500, duration = 101, value = 0.1), steps 500 to 600 as in
#   their code;
# * immigration = "multinomial": the emigrants of each patch are distributed among the other
#   patches by the kernel. Their code pools the emigrants of a species and sends each to a patch
#   with probability proportional to the expected inflow, so an emigrant can return to its own
#   patch. With a kernel this flat (exp(-0.1 d) on a 10 x 10 torus, at most 7.1 cells apart) the
#   two differ little;
# * seed_steps = numeric(0) and seed_mean = 1 for the Poisson(1) start.
#
# Their Fig. 2 uses n = 0.25, d = 1e-5, 1e-4 and 1e-3, step 700 (100 steps after the pulse);
# their Fig. 3 uses n = 0.2, 0.3 and 0.4 and steps 601 to 1000. Their archived code runs 50
# replicates per setting; the legends say 100. We run N_REP (default 100) per setting.
# Run from the root of the package. N_REP, CORES and OUTPUT change the size and the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_rep <- as.integer(Sys.getenv("N_REP", "100"))
cores <- as.integer(Sys.getenv("CORES", "6"))
output <- Sys.getenv("OUTPUT", "data-raw/replication_hajian2026.rds")

n_steps <- 1000
record <- c(499, 601, seq(610, n_steps, 10))

one_run <- function(dispersal, n_disturbed, niche, seed) {
  set.seed(seed)
  k <- 30
  sp <- paste0("sp", seq_len(k))
  landscape <- landscape_grid(side = 10, disp_scale = 10, torus = TRUE)
  env <- stats::runif(100, 0.5, 1)
  traits <- data.frame(optimum = stats::runif(k, 0.5, 1), sigma = sqrt(2) * niche, r_max = 5,
                       row.names = sp)
  alpha <- matrix(stats::runif(k * k, 0, 1), k, k, dimnames = list(sp, sp)) * 0.05
  diag(alpha) <- 0.5 * 0.05
  env_time <- env_pulse(env, n_steps, sites = sample(100, n_disturbed), start = 500,
                        duration = 101, value = 0.1)
  sim <- simulate_metacommunity(landscape, traits, alpha, dispersal = dispersal,
                                n_steps = n_steps, env_time = env_time,
                                seed_steps = numeric(0), seed_mean = 1,
                                immigration = "multinomial", record_steps = record)
  present <- sim$history > 0
  a <- apply(present, 3, function(x) mean(rowSums(x)))
  g <- apply(present, 3, function(x) sum(colSums(x) > 0))
  data.frame(step = record[-1], alpha_pre = a[1], gamma_pre = g[1], alpha = a[-1],
             gamma = g[-1], d_alpha = log(a[-1] / a[1]), d_gamma = log(g[-1] / g[1]),
             row.names = NULL)
}

grid <- expand.grid(rep = seq_len(n_rep), n_disturbed = c(1, seq(5, 100, 5)),
                    dispersal = c(1e-5, 1e-4, 1e-3), niche = c(0.2, 0.25, 0.3, 0.4))
grid$seed <- 2026 + seq_len(nrow(grid))
started <- Sys.time()
runs <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  out <- one_run(grid$dispersal[i], grid$n_disturbed[i], grid$niche[i], grid$seed[i])
  cbind(grid[i, c("niche", "dispersal", "n_disturbed", "rep")], out, row.names = NULL)
}, mc.cores = cores)
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
sims <- do.call(rbind, runs)
sims$d_beta <- sims$d_gamma - sims$d_alpha

# control without disturbance (not in the article, whose smallest extent is one patch), with
# its own seeds so that the grid above is unchanged: changes in it are the background drift of
# a metacommunity that is still losing species at step 499
control_grid <- expand.grid(rep = seq_len(n_rep), dispersal = c(1e-5, 1e-4, 1e-3),
                            niche = c(0.2, 0.25, 0.3, 0.4))
control_grid$seed <- 9026 + seq_len(nrow(control_grid))
control <- do.call(rbind, parallel::mclapply(seq_len(nrow(control_grid)), function(i) {
  out <- one_run(control_grid$dispersal[i], 0, control_grid$niche[i], control_grid$seed[i])
  cbind(control_grid[i, c("niche", "dispersal", "rep")], n_disturbed = 0, out,
        row.names = NULL)
}, mc.cores = cores))
control$d_beta <- control$d_gamma - control$d_alpha
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

# replicates in which the whole metacommunity went extinct (gamma = 0, log change -Inf) are
# counted in `n_extinct` and left out of the means, as in the script of their Fig. 3
summarise_runs <- function(d) {
  keys <- c("niche", "dispersal", "n_disturbed", "step")
  d$extinct <- d$gamma == 0
  out <- stats::aggregate(cbind(n = 1, n_extinct = d$extinct), d[keys], sum)
  ok <- d[!d$extinct, ]
  for (v in c("d_alpha", "d_gamma", "d_beta")) {
    m <- stats::aggregate(ok[v], ok[keys], mean)
    se <- stats::aggregate(ok[v], ok[keys], function(x) stats::sd(x) / sqrt(length(x)))
    names(se)[5] <- paste0(v, "_se")
    out <- merge(merge(out, m, all.x = TRUE), se, all.x = TRUE)
  }
  out[order(out$niche, out$dispersal, out$n_disturbed, out$step), ]
}
summary_table <- summarise_runs(sims)
control_table <- summarise_runs(control)
saveRDS(list(summary = summary_table, control = control_table, n_rep = n_rep, minutes = minutes,
             session = utils::sessionInfo()), output)
cat(sprintf("%d runs in %.1f min; saved to %s\n", nrow(grid) + nrow(control_grid), minutes,
            output))
