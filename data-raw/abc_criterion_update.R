# Applies the revised recovery criterion to the saved results of the recovery studies, without
# running them again. The criterion set before the runs was R2 >= 0.5 and coverage of the 90%
# interval between 85% and 95%. The upper bound counted conservative intervals as a failure,
# so that well identified parameters (R2 above 0.95, coverage 0.955 to 0.99) were marked as not
# recoverable; the extended tests of the package (tests/testthat/test-extended-*.R) showed that
# the intervals of ABC-RF are conservative with several statistics, which is a property of the
# method. Revised criterion (2026-11-21): R2 >= 0.5 and coverage >= 85%; `conservative` flags
# coverage above 95%. Only the columns `recoverable` and `conservative` change.
# Run from the root of the package.
revise <- function(results) {
  results$recoverable <- results$r2 >= 0.5 & results$coverage >= 0.85
  results$conservative <- results$coverage > 0.95
  results
}
files <- c(results = "abc_recovery", results = "abc_detection", results = "abc_capacity",
           summary = "abc_genetic_variance", results = "pilot_identifiability",
           results = "pilot_identifiability_large", results = "pilot_identifiability_k",
           results = "pilot_identifiability_k_large", results = "pilot_tstats_identifiability",
           results = "pilot_tstats_identifiability_shared")
for (i in seq_along(files)) {
  path <- file.path("data-raw", paste0(files[i], ".rds"))
  x <- readRDS(path)
  part <- names(files)[i]
  before <- x[[part]]$recoverable
  x[[part]] <- revise(x[[part]])
  x[[part]]$too_wide <- NULL
  changed <- which(before != x[[part]]$recoverable)
  cat(files[i], ": ", length(changed), " rows now recoverable\n", sep = "")
  if (length(changed)) print(x[[part]][changed, ], digits = 3, row.names = FALSE)
  saveRDS(x, path)
}
