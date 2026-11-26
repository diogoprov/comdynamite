# Why the coverage of abc_sequential() falls over the rounds in the toy case of
# data-raw/abc_sequential.R (90% intervals of a: 0.99, 0.905, 0.90, 0.85 in rounds 1 to 4).
# Item (xii) of block C of dev/ROADMAP.md. The same 200 pseudo-observations and seeds, four
# rounds of 1000 simulations, and these variants:
#   forest       abc_sequential() as it is (min_node_size 5, kernel of twice the variance);
#   oracle       the same rounds with exact weights, the likelihood of the observation, in place
#                of the weights of the forests: checks the proposal and the importance weights;
#   no_is        forests without the importance correction (prior / proposal);
#   one_param    a alone, with s and a statistic of pure noise, without b and t;
#   mns50        forests with min_node_size 50;
#   kernel_4, kernel_8, kernel_16   kernels of 4, 8 and 16 times the weighted variance;
#   local_cal    forest, with the last round recalibrated locally: 200 of its simulations held
#                out, the PIT of their true values read from the forests and the importance
#                weights, and the interval taken at the quantiles of the PIT weighted by
#                prior / proposal (after Prangle et al. 2014).
# The wider kernels replace the internal kernel_sd() of the package for the run.
# Run from the root of the package, with abcrf installed: Rscript
# data-raw/abc_sequential_diagnostics.R. CORES changes the threads; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
cores <- as.integer(Sys.getenv("CORES", "16"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_sequential_diagnostics.rds")
n_toy <- 200
n_round <- 1000
rounds <- 4
n_cal <- 200
ntree <- 500
prior <- list(a = c(-10, 10), b = c(0, 5))
sim <- function(scenario, shared) {
  c(s = scenario$a + stats::rnorm(1, 0, 0.1), t = scenario$b + stats::rnorm(1, 0, 1))
}
summ <- function(sim, shared) sim
set.seed(2061)
truth <- data.frame(a = stats::runif(n_toy, -10, 10), b = stats::runif(n_toy, 0, 5))
pseudo <- cbind(truth, s = truth$a + stats::rnorm(n_toy, 0, 0.1),
                t = truth$b + stats::rnorm(n_toy, 0, 1))
seed_of <- function(i) 20000 + 10 * i

# a kernel of `m` times the weighted variance, in place of the package's twice
set_kernel <- function(m) {
  ns <- asNamespace("comdynamite")
  f <- function(x, w) {
    mu <- sum(w * x)
    sqrt(m * sum(w * (x - mu)^2))
  }
  environment(f) <- list2env(list(m = m), parent = ns)
  unlockBinding("kernel_sd", ns)
  assign("kernel_sd", f, envir = ns)
  lockBinding("kernel_sd", ns)
}

# one row per round and parameter: PIT of the truth, coverage, width and effective size
score_rounds <- function(weights, i) {
  do.call(rbind, lapply(seq_along(weights), function(r) {
    w <- weights[[r]]
    do.call(rbind, lapply(intersect(c("a", "b"), names(w)), function(p) {
      x <- w[[p]]
      ww <- w[[paste0("w_", p)]]
      q <- weighted_quantile(x, ww, c(0.05, 0.1, 0.9, 0.95))
      tr <- truth[i, p]
      pit <- sum(ww[x < tr])
      data.frame(observation = i, round = r, parameter = p, pit = pit,
                 lost = pit < 0.005 || pit > 0.995,
                 coverage80 = tr >= q[2] && tr <= q[3], coverage90 = tr >= q[1] && tr <= q[4],
                 width90 = q[4] - q[1], ess = 1 / sum(ww^2))
    }))
  }))
}
forest_rounds <- function(i, min_node_size = 5, model = c("toy", "one_param")) {
  model <- match.arg(model)
  fit <- if (model == "toy") {
    suppressWarnings(abc_sequential(prior, sim, pseudo[i, c("s", "t")], n_round, rounds,
                                    summarise = summ, ntree = ntree,
                                    min_node_size = min_node_size, seed = seed_of(i)))
  } else {
    sim1 <- function(scenario, shared) {
      c(s = scenario$a + stats::rnorm(1, 0, 0.1), u = stats::rnorm(1))
    }
    suppressWarnings(abc_sequential(prior["a"], sim1, data.frame(s = pseudo$s[i], u = 0),
                                    n_round, rounds, summarise = summ, ntree = ntree,
                                    seed = seed_of(i)))
  }
  lapply(fit$rounds, `[[`, "weights")
}
oracle_rounds <- function(i) {
  cell_seed(seed_of(i))
  previous <- NULL
  out <- vector("list", rounds)
  for (r in seq_len(rounds)) {
    theta <- if (r == 1) {
      data.frame(a = stats::runif(n_round, -10, 10), b = stats::runif(n_round, 0, 5))
    } else {
      sequential_proposal(previous, prior, n_round)
    }
    lik <- stats::dnorm(pseudo$s[i], theta$a, 0.1) * stats::dnorm(pseudo$t[i], theta$b, 1)
    for (p in c("a", "b")) {
      w <- lik
      if (!is.null(previous)) {
        w <- w / sequential_density(theta[[p]], previous[[p]], previous[[paste0("w_", p)]])
      }
      theta[[paste0("w_", p)]] <- w / sum(w)
    }
    out[[r]] <- theta
    previous <- theta
  }
  out
}
# the rounds after the first, each fitted on `table` and weighted with or without correction
next_round <- function(i, previous, r, correct = TRUE, hold_out = 0) {
  cell_seed(seed_of(i) + 2 + r)
  theta <- sequential_proposal(previous, prior, n_round)
  table <- run_scenarios(theta, sim, NULL, summ, seed = seed_of(i) + r * 1e7)
  table$replicate <- NULL
  held <- table[seq_len(hold_out), , drop = FALSE]
  if (hold_out > 0) table <- table[-seq_len(hold_out), ]
  fit <- abc_fit(table, parameters = c("a", "b"), statistics = c("s", "t"), ntree = ntree)
  weigh <- function(p, obs) {
    w <- parameter_weights(fit, p, obs)
    if (correct) {
      w <- w / sequential_density(fit$training[[p]], previous[[p]],
                                  previous[[paste0("w_", p)]])
    }
    sweep(w, 2, colSums(w), "/")
  }
  out <- fit$training[c("a", "b")]
  for (p in c("a", "b")) out[[paste0("w_", p)]] <- weigh(p, pseudo[i, c("s", "t")])[, 1]
  attr(out, "weigh") <- weigh
  attr(out, "held") <- held
  out
}
no_is_rounds <- function(i) {
  first <- forest_rounds(i)[[1]]
  out <- list(first)
  for (r in 2:rounds) out[[r]] <- next_round(i, out[[r - 1]], r, correct = FALSE)
  out
}
# local recalibration of the last round: the interval of each parameter at the quantiles of
# the PIT of the held-out simulations, weighted by prior / proposal
local_cal <- function(i) {
  weights <- forest_rounds(i)[seq_len(rounds - 1)]
  previous <- weights[[rounds - 1]]
  last <- next_round(i, previous, rounds, hold_out = n_cal)
  held <- attr(last, "held")
  do.call(rbind, lapply(c("a", "b"), function(p) {
    x <- last[[p]]
    wc <- attr(last, "weigh")(p, held[c("s", "t")])
    pit <- colSums(wc * outer(x, held[[p]], "<"))
    iw <- 1 / sequential_density(held[[p]], previous[[p]], previous[[paste0("w_", p)]])
    probs <- pmin(pmax(weighted_quantile(pit, iw, c(0.05, 0.1, 0.9, 0.95)), 0.001), 0.999)
    q <- weighted_quantile(x, last[[paste0("w_", p)]], probs)
    tr <- truth[i, p]
    data.frame(observation = i, round = rounds, parameter = p, pit = NA_real_, lost = NA,
               coverage80 = tr >= q[2] && tr <= q[3], coverage90 = tr >= q[1] && tr <= q[4],
               width90 = q[4] - q[1], ess = 1 / sum(last[[paste0("w_", p)]]^2))
  }))
}

run_variant <- function(name, fun, kernel = 2) {
  set_kernel(kernel)
  started <- Sys.time()
  rows <- parallel::mclapply(seq_len(n_toy), function(i) {
    out <- fun(i)
    if (is.data.frame(out)) out else score_rounds(out, i)
  }, mc.cores = cores)
  cat(name, round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "min\n")
  cbind(variant = name, do.call(rbind, rows))
}

if (sys.nframe() == 0) {
  started <- Sys.time()
  d <- rbind(run_variant("forest", forest_rounds),
             run_variant("oracle", oracle_rounds),
             run_variant("no_is", no_is_rounds),
             run_variant("one_param", function(i) forest_rounds(i, model = "one_param")),
             run_variant("mns50", function(i) forest_rounds(i, min_node_size = 50)),
             run_variant("kernel_4", forest_rounds, kernel = 4),
             run_variant("kernel_8", forest_rounds, kernel = 8),
             run_variant("kernel_16", forest_rounds, kernel = 16),
             run_variant("local_cal", local_cal))
  set_kernel(2)
  results <- stats::aggregate(cbind(coverage80, coverage90, width90, ess, lost) ~ variant +
                                parameter + round, d, mean, na.action = stats::na.pass)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(list(results = results, rows = d, minutes = minutes, date = Sys.Date()), output,
          compress = "xz")
  print(results[results$round == rounds, ], digits = 3, row.names = FALSE)
  cat(round(minutes, 1), "min\n")
}
