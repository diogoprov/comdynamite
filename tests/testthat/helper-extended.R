# Extended tests (rOpenSci standards G5.10-G5.12) check statistical properties over many
# simulated data sets: type I error and power of the permutation tests, recovery of known
# parameters over several seeds, sensitivity to noise and convergence as data grow. They take
# minutes, so they run only when the environment variable COMDYNAMITE_EXTENDED_TESTS is "true";
# see tests/README.md.
skip_if_not_extended <- function() {
  testthat::skip_if_not(identical(Sys.getenv("COMDYNAMITE_EXTENDED_TESTS"), "true"),
                        "extended tests run with COMDYNAMITE_EXTENDED_TESTS=true")
}

# TRUE when the rejection rate `rejections / n` is compatible with the nominal `level`, that is,
# when `level` falls inside the exact binomial interval of the observed rate at `confidence`
rejection_compatible <- function(rejections, n, level = 0.05, confidence = 0.99) {
  ci <- stats::binom.test(rejections, n, conf.level = confidence)$conf.int
  level >= ci[1] && level <= ci[2]
}
