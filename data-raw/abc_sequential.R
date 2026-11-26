# Sequential ABC with random forests (abc_sequential(), the ABC-SMC-RF of Dinh et al. 2025,
# doi:10.1007/s11222-025-10748-x) against ABC-RF in one round with the same number of
# simulations. Item (xii) of block C of dev/ROADMAP.md. Dinh et al. report no coverage of the
# credible intervals, so the procedure is repeated on pseudo-observations with known parameters.
# Two cases:
#   toy        a ~ U(-10, 10) and b ~ U(0, 5), with s = a + N(0, 0.1) and t = b + N(0, 1): a wide
#              prior for a, a prior about as wide as the posterior for b, and posteriors known
#              (truncated normals). N_TOY pseudo-observations from the prior; ROUNDS rounds of
#              N_ROUND simulations against one round of ROUNDS * N_ROUND.
#   community  the design of data-raw/abc_recovery.R (tree of 30 species, grid of 25 sites,
#              series of five surveys; log_dispersal, inter_max, sigma_mean and lambda), with
#              N_COMMUNITY pseudo-observations drawn from its table. The first round is
#              N_FIRST rows of that table, and rounds 2 to ROUNDS simulate N_FIRST each; one
#              round of ROUNDS * N_FIRST rows of the table is the comparison.
# Criterion set before running: at the 90% level, intervals at least 20% narrower than in one
# round (toy: or nearer the width of the exact posterior), with coverage of 85% or more. With
# few pseudo-observations the coverage of the community case is only indicative.
# Run from the root of the package, with abcrf installed, after data-raw/abc_recovery.R:
# Rscript data-raw/abc_sequential.R. CASES (toy,community), N_TOY, N_ROUND, N_COMMUNITY,
# N_FIRST, ROUNDS, CORES and NTREE change the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
cases <- strsplit(Sys.getenv("CASES", "toy,community"), ",")[[1]]
n_toy <- as.integer(Sys.getenv("N_TOY", "200"))
n_round <- as.integer(Sys.getenv("N_ROUND", "1000"))
n_community <- as.integer(Sys.getenv("N_COMMUNITY", "20"))
n_first <- as.integer(Sys.getenv("N_FIRST", "2000"))
rounds <- as.integer(Sys.getenv("ROUNDS", "4"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_sequential.rds")
levels <- c(0.5, 0.8, 0.9, 0.95)

# one row per pseudo-observation, parameter, level and method: interval and truth
intervals_sequential <- function(fit, truth, obs) {
  do.call(rbind, lapply(levels, function(l) {
    do.call(rbind, lapply(seq_along(fit$rounds), function(r) {
      e <- predict(fit, level = l, round = r)
      data.frame(observation = obs, parameter = e$parameter, level = l,
                 method = paste0("sequential_", r), estimate = e$estimate, low = e$low,
                 high = e$high, truth = unlist(truth[e$parameter]), ess = e$ess)
    }))
  }))
}
intervals_single <- function(fit, pseudo, method) {
  do.call(rbind, lapply(levels, function(l) {
    e <- predict(fit, pseudo, level = l, ncores = cores)
    data.frame(observation = e$observation, parameter = e$parameter, level = l,
               method = method, estimate = e$estimate, low = e$low, high = e$high,
               truth = mapply(function(i, p) pseudo[[p]][i], e$observation, e$parameter),
               ess = NA_real_)
  }))
}
score <- function(d) {
  d$covered <- d$truth >= d$low & d$truth <= d$high
  d$width <- d$high - d$low
  d$abs_error <- abs(d$estimate - d$truth)
  out <- stats::aggregate(cbind(covered, width, abs_error, ess) ~ case + parameter + level +
                            method, d, mean, na.action = stats::na.pass)
  names(out)[names(out) == "covered"] <- "coverage"
  n <- stats::aggregate(observation ~ case + parameter + level + method, d, length)
  out$se <- sqrt(out$level * (1 - out$level) / n$observation)
  out
}

toy_prior <- list(a = c(-10, 10), b = c(0, 5))
toy_sd <- c(a = 0.1, b = 1)
toy_sim <- function(scenario, shared) {
  c(s = scenario$a + stats::rnorm(1, 0, 0.1), t = scenario$b + stats::rnorm(1, 0, 1))
}
toy_summary <- function(sim, shared) sim
# width of the central interval of the exact posterior, a normal around the statistic
# truncated to the prior
exact_width <- function(stat, sd, limits, level) {
  lo <- stats::pnorm(limits[1], stat, sd)
  hi <- stats::pnorm(limits[2], stat, sd)
  q <- stats::qnorm(lo + (hi - lo) * c((1 - level) / 2, (1 + level) / 2), stat, sd)
  q[2] - q[1]
}

run_toy <- function() {
  set.seed(2061)
  truth <- data.frame(a = stats::runif(n_toy, -10, 10), b = stats::runif(n_toy, 0, 5))
  pseudo <- cbind(truth, s = truth$a + stats::rnorm(n_toy, 0, 0.1),
                  t = truth$b + stats::rnorm(n_toy, 0, 1))
  sequential <- parallel::mclapply(seq_len(n_toy), function(i) {
    fit <- abc_sequential(toy_prior, toy_sim, observed = pseudo[i, c("s", "t")], n = n_round,
                          rounds = rounds, summarise = toy_summary, ntree = ntree,
                          seed = 20000 + 10 * i)
    intervals_sequential(fit, truth[i, ], i)
  }, mc.cores = cores)
  single_table <- abc_table(toy_prior, toy_sim, n = rounds * n_round, summarise = toy_summary,
                            seed = 2062)
  single <- abc_fit(single_table, ntree = ntree, ncores = cores)
  d <- rbind(do.call(rbind, sequential), intervals_single(single, pseudo, "single"))
  stat <- c(a = "s", b = "t")
  d$exact_width <- mapply(function(o, p, l) {
    exact_width(pseudo[[stat[[p]]]][o], toy_sd[[p]], toy_prior[[p]], l)
  }, d$observation, d$parameter, d$level)
  cbind(case = "toy", d)
}

run_community <- function() {
  recovery <- readRDS("data-raw/abc_recovery.rds")$table
  set.seed(20261107)
  fixed <- list(tree = treeshapesim::sim_tree_shape(30), landscape = landscape_grid(side = 5))
  prior <- attr(recovery, "prior")
  one_run <- function(scenario, shared) {
    simulate_metacommunity(shared$landscape,
                           species_traits(shared$tree, scenario$lambda,
                                          sigma_mean = scenario$sigma_mean, sigma_sd = 2),
                           competition_matrix(shared$tree, "phylo",
                                              inter_max = scenario$inter_max),
                           dispersal = 10^scenario$log_dispersal, n_steps = 1000,
                           record_steps = seq(800, 1000, by = 50))
  }
  summarise <- function(sim, shared) summary_stats(sim, shared$tree, n_null = 49)
  parameters <- names(prior)
  statistics <- attr(recovery, "statistics")
  data <- as.data.frame(recovery)
  set.seed(2063)
  rows <- sample.int(nrow(data), n_community + rounds * n_first)
  pseudo <- data[rows[seq_len(n_community)], ]
  pool <- rows[-seq_len(n_community)]
  as_table <- function(d) {
    structure(d, class = c("mc_abc_table", "data.frame"), prior = prior,
              parameters = parameters, statistics = statistics)
  }
  first <- as_table(data[pool[seq_len(n_first)], ])
  single <- abc_fit(as_table(data[pool, ]), ntree = ntree, ncores = cores)
  future::plan(future::multicore, workers = cores)
  on.exit(future::plan(future::sequential))
  sequential <- lapply(seq_len(n_community), function(i) {
    started <- Sys.time()
    fit <- abc_sequential(prior, one_run, observed = pseudo[i, statistics], n = n_first,
                          rounds = rounds, setup = function() fixed, summarise = summarise,
                          table = first, ntree = ntree, ncores = cores, seed = 30000 + 10 * i)
    cat("pseudo-observation", i, "of", n_community, round(as.numeric(difftime(
      Sys.time(), started, units = "mins")), 1), "min\n")
    intervals_sequential(fit, pseudo[i, parameters], i)
  })
  d <- rbind(do.call(rbind, sequential), intervals_single(single, pseudo, "single"))
  cbind(case = "community", d, exact_width = NA_real_)
}

if (sys.nframe() == 0) {
  started <- Sys.time()
  d <- do.call(rbind, lapply(cases, function(case) {
    switch(case, toy = run_toy(), community = run_community())
  }))
  results <- score(d)
  exact <- stats::aggregate(exact_width ~ case + parameter + level + method, d, mean)
  results <- merge(results, exact, all.x = TRUE)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(list(results = results, intervals = d, minutes = minutes,
               settings = list(n_toy = n_toy, n_round = n_round, n_community = n_community,
                               n_first = n_first, rounds = rounds, ntree = ntree),
               date = Sys.Date()), output, compress = "xz")
  shown <- results[results$level == 0.9 & results$method %in%
                     c("single", "sequential_1", paste0("sequential_", rounds)), ]
  print(shown[order(shown$case, shown$parameter, shown$method), ], digits = 3,
        row.names = FALSE)
  cat(round(minutes, 1), "min\n")
}
