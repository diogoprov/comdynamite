# Follow-up of data-raw/abc_npe_pilot.R: can the overconfidence of neural posterior estimation (NPE)
# be corrected without losing its narrower intervals? Item (xiii) of block C of dev/ROADMAP.md.
# In the pilot, NPE (masked autoregressive flow, neuralsbi) gave the same or better R2 than ABC-RF
# with intervals 25% to 47% narrower, but 90% intervals of G and Ne/N held the truth in only 76% to
# 83% of pseudo-observations with 10 000 simulations and eight parameters. Variants, on the same
# tables, training rows and pseudo-observations as the pilot:
#   maf          the pilot's estimator;
#   maf_recal    global recalibration: 500 training simulations held out, the PIT of their true
#                values read from their NPE posteriors, and intervals taken at the quantiles of
#                the PIT instead of the nominal probabilities (as abc_calibrate() does for
#                ABC-RF); the network is trained on the other simulations;
#   nsf, mdn     other density estimators (neural spline flow, mixture density network);
#   ens5         five MAF networks with different seeds, their draws pooled;
#   ens5_recal   the ensemble, recalibrated as maf_recal;
#   maf_2500, maf_5000   MAF trained on the first 2500 or 5000 training rows (case g_ne).
# Cases: toy (maf, maf_recal), community with 10 000 training simulations and g_ne (both
# environments) with all variants. Criterion set before running: a variant wins if its 90%
# intervals hold the truth in 85% to 95% of pseudo-observations for every parameter of interest
# and are narrower than those of ABC-RF in the pilot.
# Run from the root of the package as data-raw/abc_npe_pilot.R, with the same TABLES_DIR:
# TABLES_DIR=<folder> Rscript data-raw/abc_npe_tests.R. CASES, N_CAL, N_DRAWS, CORES and OUTPUT
# change the run; CACHE_DIR, if set, keeps the result of each case and of each variant, so that
# an interrupted run starts again where it stopped.
source("data-raw/abc_npe_pilot.R")
cases <- strsplit(Sys.getenv("CASES", "toy,community,g_ne"), ",")[[1]]
n_cal <- as.integer(Sys.getenv("N_CAL", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_npe_tests.rds")
cache_dir <- Sys.getenv("CACHE_DIR")
cached <- function(key, expr) {
  if (!nzchar(cache_dir)) return(expr)
  file <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  out <- expr
  saveRDS(out, file)
  out
}

fit_npe <- function(train, parameters, statistics, prior, seed, estimator = "maf") {
  neuralsbi::npe(
    neuralsbi::prior_uniform(vapply(prior[parameters], `[`, numeric(1), 1),
                             vapply(prior[parameters], `[`, numeric(1), 2)),
    theta = as.matrix(train[parameters]), x = as.matrix(train[statistics]),
    density_estimator = estimator, seed = seed
  )
}
# draws of the posterior of each row of `data` from the networks `fits`, pooled; counts the rows
# for which the estimator leaked mass outside the prior
draws_of <- function(fits, data, parameters, statistics, derived) {
  leaks <- 0
  draws <- lapply(seq_len(nrow(data)), function(i) {
    d <- do.call(rbind, lapply(fits, function(fit) {
      withCallingHandlers({
        post <- neuralsbi::posterior(fit, x_obs = unlist(data[i, statistics]))
        as.matrix(neuralsbi::sample(post, ceiling(n_draws / length(fits))))
      }, warning = function(w) {
        if (grepl("prior support", conditionMessage(w))) {
          leaks <<- leaks + 1
          invokeRestart("muffleWarning")
        }
      })
    }))
    colnames(d) <- parameters
    if (!is.null(derived)) d <- cbind(d, derived(d))
    d
  })
  attr(draws, "leaks") <- leaks
  attr(draws, "nan") <- sum(vapply(draws, function(d) sum(!stats::complete.cases(d)), numeric(1)))
  draws
}
# PIT of each column of the draws at the truth, over the calibration simulations
pit_of <- function(draws, data) {
  ps <- colnames(draws[[1]])
  stats::setNames(lapply(ps, function(p) {
    vapply(seq_along(draws), function(i) mean(draws[[i]][, p] < data[[p]][i], na.rm = TRUE),
           numeric(1))
  }), ps)
}
intervals_of <- function(draws, test, variant, pit = NULL) {
  do.call(rbind, lapply(seq_along(draws), function(i) {
    d <- draws[[i]]
    do.call(rbind, lapply(colnames(d), function(p) {
      probs <- c((1 - levels) / 2, (1 + levels) / 2)
      if (!is.null(pit)) {
        probs <- pmin(pmax(stats::quantile(pit[[p]], probs, type = 8, names = FALSE,
                                           na.rm = TRUE), 0.001), 0.999)
      }
      q <- stats::quantile(d[, p], c(0.5, probs), names = FALSE, na.rm = TRUE)
      data.frame(observation = i, parameter = p, level = levels, method = variant,
                 estimate = q[1], low = q[1 + seq_along(levels)],
                 high = q[1 + length(levels) + seq_along(levels)], truth = test[[p]][i],
                 cor_G_ne = if (all(c("log_G", "log_ne_ratio") %in% colnames(d))) {
                   stats::cor(d[, "log_G"], d[, "log_ne_ratio"], use = "complete.obs")
                 } else NA_real_)
    }))
  }))
}

variants_all <- c("maf", "maf_recal", "nsf", "mdn", "ens5", "ens5_recal", "maf_2500", "maf_5000")
run_variants <- function(case, train, test, parameters, statistics, prior, seed, variants,
                         derived = NULL) {
  statistics <- usable(train, statistics)
  train <- train[complete_rows(train, c(parameters, statistics)), ]
  test <- test[complete_rows(test, c(parameters, statistics)), ]
  cal <- train[seq_len(n_cal), ]
  rest <- train[-seq_len(n_cal), ]
  one <- function(variant) {
    started <- Sys.time()
    recal <- grepl("_recal$", variant)
    base <- sub("_recal$", "", variant)
    data <- if (recal) rest else train
    if (grepl("^maf_[0-9]+$", base)) {
      data <- train[seq_len(as.integer(sub("maf_", "", base))), ]
      base <- "maf"
    }
    estimator <- if (base == "ens5") "maf" else base
    n_nets <- if (base == "ens5") 5 else 1
    fits <- lapply(seq_len(n_nets), function(k) {
      fit_npe(data, parameters, statistics, prior, seed + 100 * k, estimator)
    })
    draws <- draws_of(fits, test, parameters, statistics, derived)
    cal_draws <- if (recal) draws_of(fits, cal, parameters, statistics, derived)
    pit <- if (recal) pit_of(cal_draws, cal)
    out <- intervals_of(draws, test, variant, pit)
    minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    nan <- attr(draws, "nan") + if (recal) attr(cal_draws, "nan") else 0
    cat(case, variant, round(minutes, 1), "min,", attr(draws, "leaks"), "leaks,", nan,
        "draws with NaN\n")
    list(intervals = cbind(case = case, out),
         timing = data.frame(case = case, variant = variant, n_train = nrow(data),
                             minutes = minutes, leaks = attr(draws, "leaks"), nan = nan))
  }
  # the recalibrated variants reuse nothing, so each trains its own networks
  runs <- lapply(variants, function(v) cached(paste(case, v, sep = "_"), one(v)))
  list(intervals = do.call(rbind, lapply(runs, `[[`, "intervals")),
       timing = do.call(rbind, lapply(runs, `[[`, "timing")))
}

tests_toy <- function() {
  prior <- list(a = c(-10, 10), b = c(0, 5))
  set.seed(2071)
  n <- 4000 + n_test
  d <- data.frame(a = stats::runif(n, -10, 10), b = stats::runif(n, 0, 5))
  d$s <- d$a + stats::rnorm(n, 0, 0.1)
  d$t <- d$b + stats::rnorm(n, 0, 1)
  run_variants("toy", d[seq_len(4000), ], d[4000 + seq_len(n_test), ], c("a", "b"),
               c("s", "t"), prior, seed = 2072, variants = c("maf", "maf_recal"))
}
tests_community <- function() {
  table <- readRDS("data-raw/abc_recovery.rds")$table
  prior <- attr(table, "prior")
  data <- as.data.frame(table)
  set.seed(2073)
  rows <- sample.int(nrow(data))
  test <- data[rows[seq_len(n_test)], ]
  rest <- rows[-seq_len(n_test)]
  run_variants(paste0("community_", n_small), data[rest[seq_len(n_small)], ], test, names(prior),
               attr(table, "statistics"), prior, seed = 2074,
               variants = setdiff(variants_all, c("maf_2500", "maf_5000")))
}
tests_g_ne <- function() {
  if (!nzchar(tables_dir)) stop("Set TABLES_DIR to the tables of abc_temporal_drift.R.")
  out <- lapply(c("static", "warming"), function(environment) {
    table <- readRDS(file.path(tables_dir, paste0(environment, "_constant_unknown.rds")))
    prior <- attr(table, "prior")
    data <- as.data.frame(table)
    data$log_G_ne <- data$log_G - data$log_ne_ratio
    statistics <- names(data)[startsWith(names(data), "c_") |
                                startsWith(names(data), "surveys2_")]
    set.seed(2076)
    rows <- sample.int(nrow(data))
    run_variants(paste0("g_ne_", environment), data[rows[-seq_len(n_test)], ],
                 data[rows[seq_len(n_test)], ], names(prior), statistics, prior, seed = 2077,
                 variants = variants_all,
                 derived = function(d) cbind(log_G_ne = d[, "log_G"] - d[, "log_ne_ratio"]))
  })
  list(intervals = do.call(rbind, lapply(out, `[[`, "intervals")),
       timing = do.call(rbind, lapply(out, `[[`, "timing")))
}

if (sys.nframe() == 0) {
  torch::torch_set_num_threads(cores)
  started <- Sys.time()
  runs <- lapply(cases, function(case) {
    switch(case, toy = tests_toy(), community = tests_community(), g_ne = tests_g_ne())
  })
  intervals <- do.call(rbind, lapply(runs, `[[`, "intervals"))
  timing <- do.call(rbind, lapply(runs, `[[`, "timing"))
  results <- score(intervals)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(list(results = results, intervals = intervals, timing = timing, minutes = minutes,
               settings = list(n_small = n_small, n_cal = n_cal, n_draws = n_draws,
                               neuralsbi = as.character(utils::packageVersion("neuralsbi"))),
               date = Sys.Date()), output, compress = "xz")
  shown <- results[results$level == 0.9, ]
  print(shown[order(shown$case, shown$parameter, shown$method),
              c("case", "parameter", "method", "r2", "coverage", "width")],
        digits = 3, row.names = FALSE)
  print(timing, digits = 3, row.names = FALSE)
  cat(round(minutes, 1), "min\n")
}
