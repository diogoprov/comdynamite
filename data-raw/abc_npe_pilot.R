# Pilot of neural posterior estimation (NPE) against ABC with random forests (ABC-RF), on reference
# tables already simulated. Item (xiii) of block C of dev/ROADMAP.md. NPE comes from the package
# neuralsbi (CRAN, MIT), with a masked autoregressive flow (`density_estimator = "maf"`) on the
# torch back end, trained once and then evaluated on every pseudo-observation (amortised). Unlike
# abcrf, it gives a joint posterior. Three cases, each with the same training simulations and the
# same pseudo-observations for both methods:
#   toy        the toy model of data-raw/abc_sequential.R: a ~ U(-10, 10), b ~ U(0, 5),
#              s = a + N(0, 0.1), t = b + N(0, 1); 4000 simulations and 200 pseudo-observations,
#              with the exact posterior known;
#   community  the reference table of data-raw/abc_recovery.R (log_dispersal, inter_max,
#              sigma_mean, lambda; series of five surveys), 200 pseudo-observations, training on
#              N_SMALL and on all the other rows;
#   g_ne       the reference tables of data-raw/abc_temporal_drift.R with constant G (static and
#              warming environment), the sampling design `surveys2` of item (viii), all eight
#              parameters in the joint posterior, 200 pseudo-observations; log10(G / (Ne/N)) is
#              taken from the joint draws for NPE and from a forest fitted to it for ABC-RF.
# Statistics with more than 5% missing or infinite values are dropped, then incomplete
# simulations, as in abc_fit(). Criterion set before running: NPE is worth a place in the package
# if (1) in `community` its R2 is at most 0.05 below that of ABC-RF, or better, with coverage of
# the 90% interval between 0.85 and 0.95 for every parameter; and (2) in `g_ne`, either the R2 of
# log_ne_ratio rises by at least 0.10 over ABC-RF, or log10(G / (Ne/N)) from the joint draws is
# within 0.05 of the R2 of the forest fitted to it directly.
# Run from the root of the package, with abcrf, torch (with torch::install_torch()) and neuralsbi
# installed, after data-raw/abc_recovery.R and data-raw/abc_temporal_drift.R run with TABLES_DIR
# set: TABLES_DIR=<that folder> Rscript data-raw/abc_npe_pilot.R. CASES (toy,community,g_ne),
# N_SMALL, N_DRAWS, CORES, NTREE and OUTPUT change the run.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
cases <- strsplit(Sys.getenv("CASES", "toy,community,g_ne"), ",")[[1]]
n_small <- as.integer(Sys.getenv("N_SMALL", "10000"))
n_draws <- as.integer(Sys.getenv("N_DRAWS", "1000"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
tables_dir <- Sys.getenv("TABLES_DIR")
output <- Sys.getenv("OUTPUT", "data-raw/abc_npe_pilot.rds")
levels <- c(0.5, 0.8, 0.9, 0.95)
n_test <- 200

# statistics usable by both methods and rows complete in them
usable <- function(data, statistics) {
  values <- as.matrix(data[statistics])
  keep <- colMeans(is.finite(values)) > 0.95 &
    apply(values, 2, function(v) isTRUE(stats::sd(v[is.finite(v)]) > 0))
  statistics[keep]
}
complete_rows <- function(data, columns) {
  rowSums(!is.finite(as.matrix(data[columns]))) == 0
}

# NPE: train once, draw from the posterior of every pseudo-observation; returns the draws as a
# list of matrices (one per pseudo-observation, columns named by parameter) and the minutes
npe_draws <- function(train, test, parameters, statistics, prior, seed) {
  started <- Sys.time()
  fit <- neuralsbi::npe(
    neuralsbi::prior_uniform(vapply(prior[parameters], `[`, numeric(1), 1),
                             vapply(prior[parameters], `[`, numeric(1), 2)),
    theta = as.matrix(train[parameters]), x = as.matrix(train[statistics]),
    density_estimator = "maf", seed = seed
  )
  minutes_train <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  draws <- lapply(seq_len(nrow(test)), function(i) {
    post <- neuralsbi::posterior(fit, x_obs = unlist(test[i, statistics]))
    d <- as.matrix(neuralsbi::sample(post, n_draws))
    colnames(d) <- parameters
    d
  })
  list(draws = draws, minutes = minutes_train,
       minutes_total = as.numeric(difftime(Sys.time(), started, units = "mins")))
}

# one row per pseudo-observation, parameter, level and method
intervals_npe <- function(draws, test, parameters, derived = NULL) {
  do.call(rbind, lapply(seq_along(draws), function(i) {
    d <- draws[[i]]
    if (!is.null(derived)) d <- cbind(d, derived(d))
    do.call(rbind, lapply(colnames(d), function(p) {
      q <- stats::quantile(d[, p], c(0.5, (1 - levels) / 2, (1 + levels) / 2), names = FALSE)
      data.frame(observation = i, parameter = p, level = levels, method = "npe",
                 estimate = q[1], low = q[1 + seq_along(levels)],
                 high = q[1 + length(levels) + seq_along(levels)], truth = test[[p]][i],
                 cor_G_ne = if (all(c("log_G", "log_ne_ratio") %in% colnames(d))) {
                   stats::cor(d[, "log_G"], d[, "log_ne_ratio"])
                 } else NA_real_)
    }))
  }))
}
intervals_rf <- function(train, test, parameters, statistics) {
  started <- Sys.time()
  fit <- abc_fit(train, parameters = parameters, statistics = statistics, ntree = ntree,
                 ncores = cores)
  out <- do.call(rbind, lapply(levels, function(l) {
    e <- predict(fit, test, level = l, ncores = cores)
    data.frame(observation = e$observation, parameter = e$parameter, level = l, method = "rf",
               estimate = e$estimate, low = e$low, high = e$high,
               truth = mapply(function(i, p) test[[p]][i], e$observation, e$parameter),
               cor_G_ne = NA_real_)
  }))
  attr(out, "minutes") <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  out
}
score <- function(d) {
  d$covered <- d$truth >= d$low & d$truth <= d$high
  d$width <- d$high - d$low
  by <- split(d, list(d$case, d$parameter, d$level, d$method), drop = TRUE)
  do.call(rbind, lapply(by, function(e) {
    data.frame(case = e$case[1], parameter = e$parameter[1], level = e$level[1],
               method = e$method[1], n = nrow(e),
               r2 = 1 - sum((e$truth - e$estimate)^2) / sum((e$truth - mean(e$truth))^2),
               coverage = mean(e$covered), width = mean(e$width),
               abs_error = mean(abs(e$estimate - e$truth)),
               cor_G_ne = mean(e$cor_G_ne))
  })) |> `rownames<-`(NULL)
}
run_case <- function(case, train, test, parameters, statistics, prior, seed, rf_parameters =
                       parameters, derived = NULL) {
  statistics <- usable(train, statistics)
  train <- train[complete_rows(train, c(parameters, statistics)), ]
  test <- test[complete_rows(test, c(parameters, statistics)), ]
  npe <- npe_draws(train, test, parameters, statistics, prior, seed)
  rf <- intervals_rf(train, test, rf_parameters, statistics)
  d <- rbind(intervals_npe(npe$draws, test, parameters, derived), rf)
  cat(case, ": NPE", round(npe$minutes, 1), "min to train,", round(npe$minutes_total, 1),
      "in all; ABC-RF", round(attr(rf, "minutes"), 1), "min\n")
  list(intervals = cbind(case = case, d),
       timing = data.frame(case = case, n_train = nrow(train), n_test = nrow(test),
                           n_statistics = length(statistics), npe_train = npe$minutes,
                           npe_total = npe$minutes_total, rf = attr(rf, "minutes")))
}

run_toy <- function() {
  prior <- list(a = c(-10, 10), b = c(0, 5))
  set.seed(2071)
  n <- 4000 + n_test
  d <- data.frame(a = stats::runif(n, -10, 10), b = stats::runif(n, 0, 5))
  d$s <- d$a + stats::rnorm(n, 0, 0.1)
  d$t <- d$b + stats::rnorm(n, 0, 1)
  run_case("toy", d[seq_len(4000), ], d[4000 + seq_len(n_test), ], c("a", "b"), c("s", "t"),
           prior, seed = 2072)
}
run_community <- function() {
  table <- readRDS("data-raw/abc_recovery.rds")$table
  prior <- attr(table, "prior")
  parameters <- names(prior)
  statistics <- attr(table, "statistics")
  data <- as.data.frame(table)
  set.seed(2073)
  rows <- sample.int(nrow(data))
  test <- data[rows[seq_len(n_test)], ]
  rest <- rows[-seq_len(n_test)]
  small <- run_case(paste0("community_", n_small), data[rest[seq_len(n_small)], ], test,
                    parameters, statistics, prior, seed = 2074)
  all <- run_case(paste0("community_", length(rest)), data[rest, ], test, parameters,
                  statistics, prior, seed = 2075)
  list(intervals = rbind(small$intervals, all$intervals),
       timing = rbind(small$timing, all$timing))
}
run_g_ne <- function() {
  if (!nzchar(tables_dir)) stop("Set TABLES_DIR to the tables of abc_temporal_drift.R.")
  out <- lapply(c("static", "warming"), function(environment) {
    table <- readRDS(file.path(tables_dir, paste0(environment, "_constant_unknown.rds")))
    prior <- attr(table, "prior")
    parameters <- names(prior)
    data <- as.data.frame(table)
    data$log_G_ne <- data$log_G - data$log_ne_ratio
    statistics <- names(data)[startsWith(names(data), "c_") |
                                startsWith(names(data), "surveys2_")]
    set.seed(2076)
    rows <- sample.int(nrow(data))
    test <- data[rows[seq_len(n_test)], ]
    train <- data[rows[-seq_len(n_test)], ]
    run_case(paste0("g_ne_", environment), train, test, parameters, statistics, prior,
             seed = 2077, rf_parameters = c("log_G", "log_ne_ratio", "log_G_ne"),
             derived = function(d) cbind(log_G_ne = d[, "log_G"] - d[, "log_ne_ratio"]))
  })
  list(intervals = do.call(rbind, lapply(out, `[[`, "intervals")),
       timing = do.call(rbind, lapply(out, `[[`, "timing")))
}

if (sys.nframe() == 0) {
  if (!requireNamespace("neuralsbi", quietly = TRUE) || !torch::torch_is_installed()) {
    stop("Needs neuralsbi and torch with torch::install_torch().")
  }
  torch::torch_set_num_threads(cores)
  started <- Sys.time()
  runs <- lapply(cases, function(case) {
    switch(case, toy = run_toy(), community = run_community(), g_ne = run_g_ne())
  })
  intervals <- do.call(rbind, lapply(runs, `[[`, "intervals"))
  timing <- do.call(rbind, lapply(runs, `[[`, "timing"))
  results <- score(intervals)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(list(results = results, intervals = intervals, timing = timing, minutes = minutes,
               settings = list(n_small = n_small, n_draws = n_draws, ntree = ntree,
                               neuralsbi = as.character(utils::packageVersion("neuralsbi")),
                               torch = as.character(utils::packageVersion("torch"))),
               date = Sys.Date()), output, compress = "xz")
  shown <- results[results$level == 0.9, ]
  print(shown[order(shown$case, shown$parameter, shown$method),
              c("case", "parameter", "method", "r2", "coverage", "width", "cor_G_ne")],
        digits = 3, row.names = FALSE)
  print(timing, digits = 3, row.names = FALSE)
  cat(round(minutes, 1), "min\n")
}
