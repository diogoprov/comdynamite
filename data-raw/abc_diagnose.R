# Does the Pareto k-hat or the effective sample size of the forest weights (abc_diagnose())
# warn about bad ABC estimates? Uses the reference table of data-raw/abc_recovery.R (50 000
# simulations, series view) and the same 200 pseudo-observations left out (seed 2028), plus the
# 100 pseudo-observations without competition of data-raw/abc_model_choice.R, which come from
# another model. Question set before running: are k-hat above its threshold, a low effective
# sample size, or statistics outside the table associated with larger errors or intervals that
# miss the truth?
# Run from the root of the package after abc_recovery.R and abc_model_choice.R. CORES and NTREE
# change the size; OUTPUT the file.
suppressMessages(devtools::load_all(".", quiet = TRUE))
cores <- as.integer(Sys.getenv("CORES", "6"))
ntree <- as.integer(Sys.getenv("NTREE", "500"))
output <- Sys.getenv("OUTPUT", "data-raw/abc_diagnose.rds")

table <- readRDS("data-raw/abc_recovery.rds")$table
statistics <- attr(table, "statistics")
cell_seed(2028) # as in abc_recovery()
test <- sort(sample.int(nrow(table), 200))
train <- table[-test, ]
attributes(train)[c("prior", "parameters", "statistics")] <-
  attributes(table)[c("prior", "parameters", "statistics")]
started <- Sys.time()
fit <- abc_fit(train, statistics = statistics, ntree = ntree, ncores = cores)

pseudo <- as.data.frame(table)[test, ]
pseudo <- pseudo[rowSums(!is.finite(as.matrix(pseudo[fit$statistics]))) == 0, ]
neutral <- readRDS("data-raw/abc_model_choice.rds")$neutral
neutral <- as.data.frame(neutral)[seq_len(nrow(neutral)) > nrow(neutral) - 100, ]
neutral$inter_max <- 0
neutral <- neutral[rowSums(!is.finite(as.matrix(neutral[fit$statistics]))) == 0, ]
sets <- list(right = pseudo, other_model = neutral)
results <- do.call(rbind, lapply(names(sets), function(set) {
  d <- sets[[set]]
  diag <- abc_diagnose(fit, d)
  est <- predict(fit, d, ncores = cores)
  est$truth <- unlist(lapply(seq_len(nrow(d)), function(i) d[i, fit$parameters]))
  cbind(set = set, merge(est, diag, by = c("observation", "parameter")))
}))
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
prior <- attr(table, "prior")
results$range <- vapply(results$parameter, function(p) diff(prior[[p]]), numeric(1))
results$error <- abs(results$estimate - results$truth) / results$range
results$covered <- results$truth >= results$low & results$truth <= results$high
results$high_k <- results$khat > results$khat_threshold

saveRDS(list(results = results, minutes = minutes, settings = list(ntree = ntree)), output)
print(aggregate(cbind(ess, khat, high_k, n_outside, error, covered) ~ set + parameter, results,
                mean), digits = 3)
cat("Spearman correlation with the absolute error, by parameter (right model):\n")
right <- results[results$set == "right", ]
print(do.call(rbind, lapply(split(right, right$parameter), function(d) {
  data.frame(parameter = d$parameter[1],
             ess = stats::cor(d$ess, d$error, method = "spearman"),
             khat = stats::cor(d$khat, d$error, method = "spearman"),
             n_outside = suppressWarnings(stats::cor(d$n_outside, d$error, method = "spearman")))
})), digits = 3)
cat("Coverage by k-hat above or below its threshold (right model):\n")
print(aggregate(covered ~ parameter + high_k, right, function(x) c(mean = mean(x), n = length(x))),
      digits = 3)
cat(round(minutes, 1), "min\n")
