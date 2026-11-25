# Reparametrising the genetic variance G and the effective size by ABC with random forests.
# Item (ix) of block C of dev/ROADMAP.md. By Lande's equations, the drift of population means has
# a variance proportional to G / Ne and the response to selection is proportional to
# G / (sigma^2 + G + V_E), so abundances and trait means should inform G and G / Ne, and Ne alone
# only as far as G is known. In data-raw/abc_identifiability.R (item viii), G was recovered with
# R2 0.91 to 0.93 and Ne/N with 0.29 to 0.61. Its reference tables were not kept, so this script
# runs them again with the same simulator, statistics and seeds (sourced from that script) and
# fits forests to functions of the parameters, which needs no new simulation design:
#   log_G_ne      log10(G / (Ne/N)) = log_G - log_ne_ratio, with constant G;
#   log_response  log10(G / (sigma^2 + G + V_E)), with constant G;
#   log_Vm_ne     log10(V_m * Ne/N) = log_Vm + log_ne_ratio, with G set by mutation, since at the
#                 balance of mutation and drift G grows with Ne * V_m.
# Questions and criterion set before running (R2 >= 0.5 and coverage of the 90% interval of at
# least 85%, in the sampling designs `base` and `best` of item viii):
# (i) are the derived parameters recovered better than Ne/N alone?
# (ii) is log_ne_ratio implied by two forests (the estimate of log_G minus that of log_G_ne, or
#     the estimate of log_Vm_ne minus that of log_Vm) better or worse than the forest fitted to
#     log_ne_ratio directly? Only R2 and bias can be given for it, since abcrf gives marginal
#     posteriors, one parameter at a time.
# Run from the root of the package, with abcrf installed: Rscript data-raw/abc_g_ne.R. N_REF,
# N_TEST, CORES, NTREE and PARTS as in data-raw/abc_identifiability.R; OUTPUT the file; TABLES_DIR,
# if set, keeps the reference tables (about 5 MB each).
source("data-raw/abc_identifiability.R")
output <- Sys.getenv("OUTPUT", "data-raw/abc_g_ne.rds")
tables_dir <- Sys.getenv("TABLES_DIR")

add_derived <- function(table, genetics) {
  if (genetics == "constant") {
    table$log_G_ne <- table$log_G - table$log_ne_ratio
    g <- 10^table$log_G
    table$log_response <- log10(g / (table$sigma^2 + g + 10^table$log_VE))
  } else {
    table$log_Vm_ne <- table$log_Vm + table$log_ne_ratio
  }
  table
}
derived_of <- function(genetics) {
  if (genetics == "constant") c("log_G_ne", "log_response") else "log_Vm_ne"
}
# log_ne_ratio implied by two forests, on the same pseudo-observations
implied_ne <- function(rec, genetics) {
  est <- attr(rec, "estimates")
  wide <- function(p) {
    e <- est[est$parameter == p, ]
    stats::setNames(e$estimate, e$name)
  }
  truth <- est[est$parameter == "log_ne_ratio", ]
  truth <- stats::setNames(truth$truth, truth$name)
  implied <- if (genetics == "constant") wide("log_G") - wide("log_G_ne") else
    wide("log_Vm_ne") - wide("log_Vm")
  implied <- implied[names(truth)]
  direct <- wide("log_ne_ratio")[names(truth)]
  score <- function(x) {
    c(r2 = 1 - sum((truth - x)^2) / sum((truth - mean(truth))^2),
      bias = mean(x - truth) / 1.3)
  }
  rbind(direct = score(direct), implied = score(implied))
}

if (sys.nframe() == 0) {
  future::plan(future::multicore, workers = cores)
  results <- list()
  for (i in parts) {
    cs <- cases[i, ]
    key <- paste(cs, collapse = "_")
    started <- Sys.time()
    shared <- list(env_time = environments[[cs$environment]])
    table <- abc_table(make_prior(cs$optima, cs$genetics), one_run, n = n_ref + n_test,
                       setup = function() shared, summarise = summarise_run, seed = 3000 + i)
    minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    if (nzchar(tables_dir)) {
      dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(table, file.path(tables_dir, paste0(key, ".rds")), compress = "xz")
    }
    statistics <- views_of(attr(table, "statistics"))[c("base", "best")]
    parameters <- c(names(make_prior(cs$optima, cs$genetics)), derived_of(cs$genetics))
    table <- add_derived(table, cs$genetics)
    recovery <- lapply(statistics, function(v) {
      abc_recovery(table, parameters = parameters, statistics = v, n_test = n_test,
                   ntree = ntree, ncores = cores, seed = 7)
    })
    implied <- lapply(recovery, implied_ne, genetics = cs$genetics)
    results[[key]] <- list(case = cs, recovery = recovery, implied = implied, minutes = minutes)
    cat("\n==", key, "(", round(minutes, 1), "min of simulation )\n")
    for (v in names(recovery)) {
      r <- recovery[[v]]
      cat("--", v, ":", paste(sprintf("%s %.2f/%.2f", r$parameter, r$r2, r$coverage),
                              collapse = "; "), "\n")
      print(implied[[v]], digits = 3)
    }
  }
  future::plan(future::sequential)
  summary_table <- do.call(rbind, lapply(results, function(x) {
    do.call(rbind, lapply(names(x$recovery), function(v) {
      data.frame(x$case, view = v,
                 x$recovery[[v]][, c("parameter", "r2", "bias", "coverage", "width")])
    }))
  }))
  summary_table$recoverable <- summary_table$r2 >= 0.5 & summary_table$coverage >= 0.85
  rownames(summary_table) <- NULL
  implied_table <- do.call(rbind, lapply(results, function(x) {
    do.call(rbind, lapply(names(x$implied), function(v) {
      data.frame(x$case, view = v, estimate = rownames(x$implied[[v]]), x$implied[[v]],
                 row.names = NULL)
    }))
  }))
  saveRDS(list(summary = summary_table, implied = implied_table,
               minutes = vapply(results, `[[`, 0, "minutes"), cases = cases,
               n_ref = n_ref, n_test = n_test, date = Sys.Date()),
          output, compress = "xz")
  message("saved ", output)
}
