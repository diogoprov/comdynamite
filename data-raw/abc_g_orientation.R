# Can the orientation of the additive genetic covariance matrix G be recovered by ABC with
# random forests from data that a field study could collect, and is the phenotypic matrix P a
# good proxy for it (Cheverud 1988)? Six species with two traits, one per environmental axis, on
# a 10 x 10 grid whose two axes are correlated (r = 0.6), so that the major axis of
# environmental variation lies at 45 degrees. G = g R(theta) diag(1, e) R(theta)', where theta
# is the angle between the leading eigenvector of G (g_max) and the major environmental axis
# (0: aligned, pi/2: orthogonal) and e the ratio of its eigenvalues. The non-genetic matrix V_E
# has its own size, angle and an eccentricity of 0.5, so that P = G + V_E departs from G.
# Parameters and uniform priors: log10 g (0 to 1.5), e (0.05 to 1), theta (0 to pi/2),
# log10 V_E (-0.5 to 1.5), theta_E (0 to pi/2), log10 dispersal (-3 to -0.5), niche breadth
# sigma (10 to 30, both axes). Each run is summarised at step 400 in two views:
#   community    abundances with both environmental axes and coordinates;
#   populations  plus 10 individuals measured per population: the pooled within-population
#                covariance (an estimate of P), the covariance among population means within
#                species (D) and the regressions of traits on the environment. Angles enter as
#                cos(2a) and sin(2a), since an axis has no sign.
# The orientation is judged by the angular error of the estimate of theta among runs with
# e <= 0.5 (where g_max is defined), against the angle of the leading eigenvector of P, the
# proxy of the Cheverud conjecture, and of D. Criterion for the other parameters as in
# data-raw/abc_recovery.R (R2 >= 0.5, coverage of the 90% interval 85% to 95%).
# Run from the root of the package, with abcrf installed. N_REF, N_TEST, CORES and NTREE change
# the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "10000"))
n_test <- as.integer(Sys.getenv("N_TEST", "300"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_g_orientation.rds")

set.seed(20261119)
side <- 10
land <- landscape_grid(side, env_bandwidth = 3, disp_scale = 1.5)
z1 <- as.numeric(scale(land$env))
z2 <- as.numeric(scale(landscape_grid(side, env_bandwidth = 3)$env))
z2 <- as.numeric(scale(0.6 * z1 + sqrt(1 - 0.6^2) * z2))
env <- cbind(a = 50 + 15 * z1, b = 50 + 15 * z2)
major <- pi / 4                                   # major axis of environmental variation
n_sp <- 6
species <- paste0("sp", seq_len(n_sp))
optima <- env[sample(nrow(env), n_sp), ]          # ancestral optima where sites exist
alpha <- matrix(0.01, n_sp, n_sp, dimnames = list(species, species))
diag(alpha) <- 0.02
n_ind <- 10
prior <- list(log_g = c(0, 1.5), e = c(0.05, 1), theta = c(0, pi / 2), log_VE = c(-0.5, 1.5),
              theta_E = c(0, pi / 2), log_dispersal = c(-3, -0.5), sigma = c(10, 30))

rotated <- function(size, ecc, angle) {
  r <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
  m <- size * r %*% diag(c(1, ecc)) %*% t(r)
  dimnames(m) <- list(c("a", "b"), c("a", "b"))
  m
}
# angle of the leading eigenvector relative to the major environmental axis, in [0, pi/2]
axis_angle <- function(m) {
  v <- eigen(m, symmetric = TRUE)$vectors[, 1]
  a <- (atan2(v[2], v[1]) - major) %% pi
  if (a > pi / 2) pi - a else a
}
shape <- function(m, prefix) {
  ev <- eigen(m, symmetric = TRUE)
  v <- ev$vectors[, 1]
  a <- atan2(v[2], v[1]) - major
  stats::setNames(c(cos(2 * a), sin(2 * a), ev$values[2] / ev$values[1], log10(sum(ev$values))),
                  paste0(prefix, c("_cos2", "_sin2", "_ecc", "_logtrace")))
}

one_run <- function(scenario, shared) {
  G <- rotated(10^scenario$log_g, scenario$e, major + scenario$theta)
  VE <- rotated(10^scenario$log_VE, 0.5, major + scenario$theta_E)
  traits <- lapply(c(a = 1, b = 2), function(k) {
    data.frame(optimum = optima[, k], sigma = scenario$sigma, r_max = sqrt(5), row.names = species)
  })
  simulate_metacommunity(land, traits, alpha, env = env, dispersal = 10^scenario$log_dispersal,
                         n_steps = 400, record_steps = 400,
                         evolution = evolution_control(G, ne_ratio = 0.5, env_var = VE))
}

summarise_views <- function(sim, shared) {
  comm <- sim$abundance
  coords <- land$coords
  ca <- summary_stats(comm, env = env[, "a"], coords = coords, n_null = 19)
  cb <- summary_stats(comm, env = env[, "b"], coords = coords, n_null = 19)
  cb <- cb[grepl("trait|env", names(cb))]
  community <- c(stats::setNames(ca, paste0("c_", names(ca))),
                 stats::setNames(cb, paste0("cb_", names(cb))))
  P <- sim$evolution$settings$genetic_var + sim$evolution$settings$env_var
  pheno <- sim$evolution$phenotype                # sites x species x traits
  within <- matrix(0, 2, 2)
  df_within <- 0
  among <- matrix(0, 2, 2)
  w_among <- 0
  slopes <- matrix(0, 2, 2)
  w_slopes <- 0
  for (j in seq_len(n_sp)) {
    occ <- which(comm[, j] > 0)
    if (length(occ) == 0) next
    means <- t(vapply(occ, function(i) {
      ind <- MASS::mvrnorm(n_ind, pheno[i, j, ], P)
      within <<- within + stats::cov(ind) * (n_ind - 1)
      colMeans(ind)
    }, numeric(2)))
    df_within <- df_within + length(occ) * (n_ind - 1)
    if (length(occ) >= 4) {
      w <- sum(comm[occ, j])
      among <- among + stats::cov(means) * w
      fit <- stats::lm(means ~ env[occ, ])
      slopes <- slopes + stats::coef(fit)[2:3, ] * w
      w_among <- w_among + w
    }
  }
  if (df_within == 0 || w_among == 0) return(c(community, p_ok = 0))
  P_hat <- within / df_within
  D <- among / w_among
  B <- slopes / w_among
  pop <- c(shape(P_hat, "P"), shape(D, "D"),
           B = as.numeric(B), D_corr = stats::cov2cor(D)[1, 2], P_corr = stats::cov2cor(P_hat)[1, 2],
           p_ok = 1)
  c(community, pop)
}

future::plan(future::multicore, workers = cores)
started <- Sys.time()
table <- abc_table(prior, one_run, n = n_ref + n_test, setup = function() list(),
                   summarise = summarise_views, seed = 2031)
future::plan(future::sequential)
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
table <- table[!is.na(table$p_ok) & table$p_ok == 1, ]
attr(table, "prior") <- prior
attr(table, "parameters") <- names(prior)
stats <- setdiff(names(table), c(names(prior), "p_ok"))
attr(table, "statistics") <- stats
views <- list(community = stats[startsWith(stats, "c")],
              populations = stats)
recovery <- lapply(views, function(v) {
  abc_recovery(table, statistics = v, n_test = n_test, ntree = ntree, ncores = cores, seed = 7)
})

# angular error of theta (degrees) where g_max is defined, against the proxies P and D
angles <- lapply(names(recovery), function(v) {
  est <- attr(recovery[[v]], "estimates")
  th <- est[est$parameter == "theta", ]
  ec <- est[est$parameter == "e", ]
  data.frame(view = v, truth = th$truth, estimate = th$estimate, e = ec$truth,
             low = th$low, high = th$high)
})
angles <- do.call(rbind, angles)
test_rows <- table[match(round(angles$truth[angles$view == "populations"], 10),
                         round(table$theta, 10)), ]
proxy <- function(cos2, sin2) {
  a <- (atan2(sin2, cos2) / 2) %% pi
  ifelse(a > pi / 2, pi - a, a)
}
pp <- data.frame(view = "P proxy", truth = test_rows$theta,
                 estimate = proxy(test_rows$P_cos2, test_rows$P_sin2), e = test_rows$e,
                 low = NA, high = NA)
dd <- data.frame(view = "D proxy", truth = test_rows$theta,
                 estimate = proxy(test_rows$D_cos2, test_rows$D_sin2), e = test_rows$e,
                 low = NA, high = NA)
angles <- rbind(angles, pp, dd)
angles$error_deg <- abs(angles$estimate - angles$truth) * 180 / pi
defined <- angles[angles$e <= 0.5, ]
angle_summary <- aggregate(error_deg ~ view, defined, function(x) c(mean = mean(x),
                                                                    median = stats::median(x)))
cat("\nSimulation:", round(minutes, 1), "min;", nrow(table), "runs kept\n")
for (v in names(recovery)) {
  cat("--", v, "\n")
  print(recovery[[v]], digits = 3, row.names = FALSE)
}
cat("\nAngular error of theta (degrees) with e <= 0.5; uniform guessing gives 30 on average\n")
print(angle_summary)
saveRDS(list(recovery = recovery, angles = angles, angle_summary = angle_summary,
             table = table, prior = prior, minutes = minutes, n_ref = n_ref, n_test = n_test,
             date = Sys.Date()), output, compress = "xz")
