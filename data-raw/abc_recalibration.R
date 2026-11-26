# Global recalibration of the credible intervals of ABC with random forests (abc_calibrate()).
# Item (xi) of block C of dev/ROADMAP.md. The intervals of predict() are approximate: with many
# statistics they are conservative (90% intervals held the truth in 92% to 99% of
# pseudo-observations), and with one statistic and small leaves too narrow (about 65%).
# abc_calibrate() reads, on simulations not used to fit, the posterior distribution function at
# the true value (the PIT), and predict(calibration = ) asks the forest for the quantiles of the
# PIT instead of the nominal probabilities. Two cases, without new simulations of the community
# model:
#   conservative  the reference table of data-raw/abc_recovery.R (tree of 30 species, grid of 25
#                 sites, series of five surveys; parameters log_dispersal, inter_max, sigma_mean
#                 and lambda), split at random into N_TRAIN simulations to fit, N_CAL to
#                 calibrate and N_TEST to test;
#   narrow        a toy table with a known posterior, as in the extended tests: a ~ U(0, 1) and
#                 one statistic s = a + N(0, 0.05), with the same split.
# Criterion set before running: after recalibration, the coverage of the 50%, 80%, 90% and 95%
# intervals of every parameter within two standard errors of the level, in N_TEST
# pseudo-observations; in the conservative case the intervals get narrower.
# Run from the root of the package, with abcrf installed, after data-raw/abc_recovery.R:
# Rscript data-raw/abc_recalibration.R. N_TRAIN, N_CAL, N_TEST, CORES and NTREE change the size;
# OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
n_train <- as.integer(Sys.getenv("N_TRAIN", "10000"))
n_cal <- as.integer(Sys.getenv("N_CAL", "500"))
n_test <- as.integer(Sys.getenv("N_TEST", "1000"))
cores <- as.integer(Sys.getenv("CORES", "16"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_recalibration.rds")
levels <- c(0.5, 0.8, 0.9, 0.95)

# coverage and mean width of the intervals of `fit` on `test`, with and without `calibration`
score <- function(fit, test, calibration, case) {
  do.call(rbind, lapply(levels, function(l) {
    do.call(rbind, lapply(list(nominal = NULL, recalibrated = calibration), function(cal) {
      est <- predict(fit, test, level = l, ncores = cores, calibration = cal)
      est$truth <- mapply(function(i, p) test[[p]][i], est$observation, est$parameter)
      do.call(rbind, lapply(split(est, est$parameter), function(e) {
        data.frame(case = case, parameter = e$parameter[1], level = l,
                   intervals = if (is.null(cal)) "nominal" else "recalibrated",
                   coverage = mean(e$truth >= e$low & e$truth <= e$high),
                   width = mean(e$high - e$low))
      }))
    }))
  }))
}
split_rows <- function(n_all, seed) {
  set.seed(seed)
  rows <- sample.int(n_all, n_train + n_cal + n_test)
  list(train = rows[seq_len(n_train)], cal = rows[n_train + seq_len(n_cal)],
       test = rows[n_train + n_cal + seq_len(n_test)])
}

if (sys.nframe() == 0) {
  started <- Sys.time()
  recovery <- readRDS("data-raw/abc_recovery.rds")$table
  rows <- split_rows(nrow(recovery), 2051)
  data <- as.data.frame(recovery)
  train <- recovery[rows$train, ]
  attributes(train)[c("prior", "parameters", "statistics")] <-
    attributes(recovery)[c("prior", "parameters", "statistics")]
  fit_a <- abc_fit(train, ntree = ntree, ncores = cores)
  cal_a <- abc_calibrate(fit_a, data[rows$cal, ], ncores = cores)
  test_a <- data[rows$test, ]
  test_a <- test_a[rowSums(!is.finite(as.matrix(test_a[fit_a$statistics]))) == 0, ]
  result_a <- score(fit_a, test_a, cal_a, "conservative")

  set.seed(2052)
  n_all <- n_train + n_cal + n_test
  toy <- data.frame(a = stats::runif(n_all))
  toy$s <- toy$a + stats::rnorm(n_all, 0, 0.05)
  rows <- split_rows(n_all, 2053)
  fit_b <- abc_fit(toy[rows$train, ], parameters = "a", statistics = "s", ntree = ntree,
                   ncores = cores)
  cal_b <- abc_calibrate(fit_b, toy[rows$cal, ], ncores = cores)
  result_b <- score(fit_b, toy[rows$test, ], cal_b, "narrow")

  results <- rbind(result_a, result_b)
  n <- c(conservative = nrow(test_a), narrow = n_test)[results$case]
  results$se <- sqrt(results$level * (1 - results$level) / n)
  results$within_2se <- abs(results$coverage - results$level) <= 2 * results$se
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(list(results = results, calibration_coverage = rbind(
    cbind(case = "conservative", cal_a$coverage), cbind(case = "narrow", cal_b$coverage)),
    minutes = minutes, settings = list(n_train = n_train, n_cal = n_cal, n_test = n_test,
                                       ntree = ntree), date = Sys.Date()),
    output, compress = "xz")
  print(results, digits = 3, row.names = FALSE)
  cat(round(minutes, 1), "min\n")
}
