#' srr_stats: General standards
#'
#' rOpenSci statistical standards (`srr`) of this category, with the parts of the
#' package that comply with each one or the reason it does not apply.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} The Description field of DESCRIPTION cites, with DOI, the papers of the main
#'   methods (Leibold et al. 2010, Peres-Neto et al. 2012, Peres-Neto & Kembel 2015, Thompson et
#'   al. 2020, Kraft et al. 2007). Each function that implements a published method lists it
#'   under the roxygen tag references, `method_references()` (R/cite.R) gives the papers by
#'   function, and inst/CITATION holds them next to the entry of the package.
#' @srrstats {G1.1} The Description field and the help pages of `node_test()`,
#'   `node_contrasts()`, `phylo_vector()` and `trait_partition()` state that these are
#'   translations to R of the MATLAB code in the supplements of Leibold et al. (2010) and
#'   Peres-Neto et al. (2012), and that the simulation models were written from the published
#'   equations. The section "Related packages" of README.md says that mcomsimr is the original R
#'   implementation of the model of Thompson et al. (2020), how the conventions of the two
#'   differ and what comdynamite adds, and places the package next to gen3sis.
#' @srrstats {G1.2} The section "Development stage" of .github/CONTRIBUTING.md is the life cycle
#'   statement: active development before version 1.0.0, interfaces that can still change
#'   between minor versions, and renamed functions and arguments kept with a deprecation warning
#'   until 1.0.0 (R/deprecated.R, `renamed_arg()`). README.md carries the repostatus badge
#'   "Active" and NEWS.md lists each change.
#' @srrstats {G1.3} Statistical terms are defined where they are used: the Details and the
#'   argument descriptions of each help page give the statistic, its null model and the meaning
#'   of each parameter (for example the T-statistics in `tstats()` and the fractions of
#'   `trait_partition()`), and the articles in vignettes/ and vignettes/articles/ go through
#'   them with examples. README.md states where conventions differ from those of the original
#'   papers (`sigma` and `alpha` in the model of Thompson et al. 2020). There is no separate
#'   glossary.
#' @srrstats {G1.4} All 131 exported functions and the S3 methods are documented with roxygen2
#'   (markdown enabled in DESCRIPTION), and each export has an alias in man/.
#' @srrstats {G1.4a} Internal functions are documented by roxygen2 blocks ending in the noRd
#'   tag, placed above their definition in every file of R/.
#' @srrstats {G1.5} The scripts that produce the results of the preprint and of the articles
#'   are in the repository with their saved outputs (`.rds`): paper/figures.R draws the figures
#'   of the manuscript from the replications in data-raw/ (`replication_leibold2019.R`,
#'   `replication_urban2009.R`, `replication_hajian2026.R`, `figures_*.R`), and the validation
#'   table of the manuscript points to the `validation_*.R` scripts of data-raw/. Both folders
#'   are in the source repository and excluded from the built package by .Rbuildignore.
#' @srrstats {G1.6} Scripts in data-raw/ compare the package with other implementations, which
#'   are only called there and are not dependencies: `validation_mcomsimr.R` and
#'   `validation_mcomsimr_time.R` (mcomsimr), `validation_damocles.R` (DAMOCLES),
#'   `validation_pcps.R` (PCPS), `validation_herodotools.R` (Herodotools),
#'   `validation_mvtraits.R` (mvMORPH), `validation_rates.R` (RRphylo, phytools, motmot,
#'   corHMM) and `validation_bamm.R` (BAMM). `benchmark_dynamics.R` times
#'   `simulate_metacommunity()`. The article vignettes/articles/validation.Rmd reports the
#'   comparison with mcomsimr.
#' @srrstats {G2.0} Scalar arguments are checked for length one by `check_count()`,
#'   `check_number()`, `check_flag()` and `check_string()` (R/checks.R), called at the top of
#'   the exported functions; data inputs are checked for matching dimensions by `check_env()`,
#'   `check_sites()` and `check_permutations()`. Tests in tests/testthat/test-checks.R and
#'   tests/testthat/test-validation-sim.R, -evol.R, -traits.R, -inference.R and -analysis.R.
#' @srrstats {G2.0a} The description of each argument states the expected length, for example
#'   "One whole number, 1 or more", "`TRUE` or `FALSE`", "one row per site" or "one value per
#'   trait", and the section "Arguments" of the package help page (R/comdynamite-package.R)
#'   says that every scalar argument is checked for its length, type and range.
#' @srrstats {G2.1} The same helpers assert types: numeric and whole for counts, numeric for
#'   numbers, logical for flags, character for strings (R/checks.R); `check_tree()` asks for a
#'   `phylo` object, and `check_comm()`, `check_env()` and `check_numeric_columns()` ask for
#'   numeric or logical columns. Tests in tests/testthat/test-checks.R and
#'   tests/testthat/test-validation-*.R.
#' @srrstats {G2.1a} The argument descriptions give the type of each vector input (numeric
#'   vector, character vector of names, matrix of sites by species, `phylo` object), and the
#'   article "Preparing your data" (vignettes/articles/preparing-data.Rmd) shows the formats
#'   the analyses expect.
#' @srrstats {G2.2} `check_count()`, `check_number()`, `check_flag()` and `check_string()` stop
#'   when an argument expected to be a single value has length other than one, with a message
#'   that names the argument; tested with vectors of length two in
#'   tests/testthat/test-checks.R and in tests/testthat/test-validation-*.R.
#' @srrstats {G2.3} See G2.3a and G2.3b.
#' @srrstats {G2.3a} Arguments that choose among options (such as `method`, `null`, `model`,
#'   `rule` and `distance`) are matched with `match.arg()`, used in 29 files of R/, so that an
#'   invalid value stops with the list of valid options.
#' @srrstats {G2.3b} The section "Arguments" of the package help page
#'   (R/comdynamite-package.R) states that these arguments are case-sensitive, with an example
#'   (`null = "variogram"` works and `null = "Variogram"` stops).
#' @srrstats {G2.4} See G2.4a to G2.4e. Conversions are explicit and happen once, in the
#'   functions that check the inputs.
#' @srrstats {G2.4a} Whole numbers given as doubles are accepted by `check_count()` and
#'   indices, counts and fold labels are converted with `as.integer()` where integers are
#'   needed (for example `fold_summarise()` in R/abc_cv.R and the range log in R/evolving.R).
#' @srrstats {G2.4b} Community and predictor tables are converted to numeric matrices by
#'   `as.matrix(x) * 1` in `check_comm()` and `check_env()`, so logical and integer input
#'   becomes double; `as.numeric()` is used for single columns.
#' @srrstats {G2.4c} Labels of sites, species and groups are converted with `as.character()`
#'   (for example in `tg_populations()`, R/trait_gradient.R, and for `kernel_group` in
#'   R/dynamics.R); `paste()` is used only to build messages and compound keys.
#' @srrstats {G2.4d} Grouping variables are converted with `factor()` and explicit levels where
#'   a factor is needed (for example `regimes` in R/damocles.R and `group` in R/drift_test.R).
#' @srrstats {G2.4e} Factor input is converted with `as.integer()` (`folds` in R/abc_cv.R) or
#'   `as.character()` (`kernel_group` in R/dynamics.R; site and species columns in
#'   R/trait_gradient.R) before use, never through its internal codes implicitly.
#' @srrstats {G2.5} Factors used as labels (`regimes` of the DAMOCLES functions, `group` of
#'   `drift_test()`, `strata` of `env_surrogates()`, and the site and species columns of the
#'   functions for individuals) are taken as unordered; the section "Arguments" of the package
#'   help says so and that the order of an ordered factor is ignored.
#' @srrstats {G2.6} One-dimensional input is accepted as a vector, a one-column matrix or a
#'   one-column data frame (`check_env()`, `metric_input()` in R/neighbour_metrics.R).
#'   `metric_input()` keeps only the `dim`, `dimnames` and `names` attributes, so a vector of
#'   class `units` gives the same neighbour metrics and `trait_range()` as a plain vector (test
#'   "vectors of class units give the same metrics as plain vectors" in tests/testthat/test-
#'   neighbour_metrics.R).
#' @srrstats {G2.7} Community, predictor, coordinate and trait tables can be a matrix, a
#'   `data.frame` or a tibble, converted by `check_comm()` and `check_env()`; data on
#'   individuals are a `data.frame` or tibble in long format. Domain forms are accepted too:
#'   `phylo` trees of ape, `dist` objects in the neighbour metrics, `simmap` histories, and
#'   `mc_data` objects built by `mc_data()` and `as_mc_data()`. The test "tibbles enter and
#'   results leave as tables" is in tests/testthat/test-parallel.R.
#' @srrstats {G2.8} Exported functions convert their inputs once, at the top: `check_comm()`
#'   and `check_env()` return numeric matrices with sites matched by name, `check_tree()`
#'   passes only `phylo` objects, and `mc_data()` and `as_mc_data()` (R/mc_data.R) gather
#'   communities, predictors, coordinates, traits and tree in one aligned object. Internal
#'   functions receive these matrices.
#' @srrstats {G2.9} Conversions that drop or add information give a warning or message:
#'   species of `comm` that are not in the tree and tips that are not in `comm`
#'   (R/phylo_beta.R, R/damocles.R), species that are not simulated (R/dynamics.R), sites
#'   without species (R/herodotools.R), sites with missing predictors (R/trait_variation.R),
#'   species left out of the Bayesian model (R/trait_gradient_bayes.R) and species names
#'   inserted when `traits` has none (R/assembly.R). Tests in
#'   tests/testthat/test-species-match.R and tests/testthat/test-assembly.R.
#' @srrstats {G2.10} Tabular input is converted to a matrix before any extraction, columns of
#'   data frames are taken with `[[`, and subsets of matrices use `drop = FALSE`, so results do
#'   not depend on the class of the table; the test with tibbles in
#'   tests/testthat/test-parallel.R compares classes of input.
#' @srrstats {G2.11} Columns of class `units` run in `check_env()`, `tstats()` and
#'   `trait_gradient()`, and vectors of that class in the neighbour metrics (tested in
#'   tests/testthat/test-neighbour_metrics.R). The units are dropped without a message, by
#'   decision: the values are used in the units given and results are in those units, as stated
#'   in the section "Arguments" of the package help.
#' @srrstats {G2.12} List columns stop with a message that names them:
#'   `check_numeric_columns()` (R/checks.R) for `comm`, `env` and coordinates, and the
#'   checks of the trait column in `tstats()` and `trait_gradient()`. Tested with list
#'   columns in `tests/testthat/test-column-classes.R`.
#' @srrstats {G2.13} Missing values are checked before any computation: `check_comm()`,
#'   `check_env()` and the checks of coordinates, traits and regimes stop with a message that
#'   names the argument (R/checks.R, R/mc_data.R, R/assembly.R, R/phylo_gradient.R,
#'   R/surrogates.R). Tests in tests/testthat/test-validation-analysis.R,
#'   tests/testthat/test-surrogates.R and tests/testthat/test-checks.R.
#' @srrstats {G2.14} The functions for individuals and the regressions on long tables
#'   (`tstats()`, `community_index()`, `trait_flex()`, `trait_rarefaction()`, `drift_test()`,
#'   `population_shifts()`) take `na_action = c("omit", "fail")`, handled by `handle_missing()`
#'   in R/checks.R; the neighbour metrics take `na.rm`. In `trait_gradient()` a missing trait
#'   value marks a population that was not measured, which the method uses by design (its help
#'   page), so it has no such argument. Tests in tests/testthat/test-checks.R and
#'   tests/testthat/test-population_shifts.R.
#' @srrstats {G2.14a} Missing values in `comm`, `env`, coordinates and species traits stop the
#'   function with an error that names the argument, as stated in the section "Arguments" of
#'   the package help page; the neighbour metrics return `NA` with the default `na.rm = FALSE`.
#' @srrstats {G2.14b} With `na_action = "omit"`, the functions of G2.14 leave out the
#'   incomplete rows with a message giving their number, and `trait_flex()` and
#'   `trait_partition()` warn when they leave out sites with missing predictors. The neighbour
#'   metrics with `na.rm = TRUE` drop missing values without a message, since they are called
#'   once per site inside `community_index()`, whose own message already counts them.
#' @srrstats {G2.15} Data reach `mean()`, `var()`, `cor()` and similar routines after the
#'   checks of G2.13, so without missing values; where missing values are allowed (trait values
#'   of individuals, statistics that cannot be computed for a site) the code selects complete
#'   cases or calls the routine with `na.rm = TRUE` (R/individuals.R, R/trait_variation.R,
#'   R/trait_gradient.R, R/neighbour_metrics.R).
#' @srrstats {G2.16} `Inf` and `-Inf` in `comm`, `env`, coordinates and traits stop the function
#'   with a message that names the argument, and `NaN` is treated as a missing value.
#'   `check_number()` and `check_count()` refuse infinite scalars unless `Inf` has a documented
#'   meaning (for example `ne_ratio = Inf`, `threshold = Inf`, `env_pulse(duration = Inf)`,
#'   `n = Inf`). The ABC functions drop simulations with non-finite statistics and say so in
#'   their help (R/abc.R). Tests in tests/testthat/test-checks.R and
#'   tests/testthat/test-validation-*.R.
#' @srrstats {G3.0} Zero variance is tested with a tolerance by `no_variance()` (R/checks.R:
#'   range at most `sqrt(.Machine$double.eps)` times the scale), used by `check_env()` and
#'   tested in tests/testthat/test-checks.R ("no_variance() ignores floating-point noise").
#'   The remaining `==` comparisons are between integers (counts, lengths, indices, codes of
#'   events), test whether a number is whole (`x == round(x)` in `check_count()`), test a
#'   value set by the user as a switch (for example `sigma == 0`, `nugget == 0`), or guard a
#'   division against an exact zero (`trait_kurtosis()`, `mst_ratio()`).
#' @srrstats {G3.1} The two functions that use a covariance matrix of the data take one from
#'   the user in place of the default estimate: `W` in `drift_test()` (default: pooled
#'   covariance within groups) and `trait_cov` in `assemble_community()` for Mahalanobis
#'   distances (default: `stats::cov()` of the pool). Any estimator can therefore be used.
#' @srrstats {G3.1a} The arguments `W` and `trait_cov` are described in the help pages of
#'   `drift_test()` and `assemble_community()`, and tests/testthat/test-drift-test.R runs
#'   `drift_test()` with a covariance matrix given by the user. The examples do not show it.
#' @srrstats {G4.0} `animate_food_web()` (R/food_web_plot.R) is the only function that writes a
#'   file; `file` is checked as one character string and the suffix `.gif` is added when the name
#'   has none (test in tests/testthat/test-food-web-plot.R).
#' @srrstats {G5.0} Tests use the data sets of the package, `daphniids` (lakes of the EMAP
#'   survey of the US EPA, as analysed by Leibold et al. 2010) in tests/testthat/test-node.R,
#'   and the worked example of Appendix A of Peres-Neto et al. (2012), typed in
#'   tests/testthat/test-phylo-gradient.R. Other tests use small data typed in the test or
#'   simulated with a fixed seed, for which the expected result is known.
#' @srrstats {G5.1} The data sets `daphniids` and `anurans` are exported (data/, documented in
#'   R/data.R, built by data-raw/daphniids.R and data-raw/anurans.R). All other test data are
#'   created inside the tests with exported functions and fixed seeds, or typed in the test
#'   files, so users can run them again.
#' @srrstats {G5.2} Errors and warnings are tested with `expect_error()`, `expect_warning()`
#'   and `expect_message()` and a pattern for the text, in the test file of each function and
#'   in tests/testthat/test-checks.R and tests/testthat/test-validation-*.R.
#' @srrstats {G5.2a} Every call to `stop()`, `warning()` and `message()` in R/ has a text
#'   that no other call shares; where the same check is made in several functions, the text
#'   names the function. Checked with a script that collects the literal text of each call.
#' @srrstats {G5.2b} Conditions are triggered by tests that match their text with
#'   `expect_error()`, `expect_warning()` or `expect_message()`. That each one is reached is
#'   checked through line coverage: codecov reports every `stop()`, `warning()` and
#'   `message()` line that no test runs, and the patch check fails on uncovered new lines.
#' @srrstats {G5.3} Tests assert that results have no missing or non-finite values with
#'   `expect_false(anyNA())` and `expect_true(all(is.finite()))`, for example in
#'   tests/testthat/test-evolution.R, test-individuals.R, test-abc.R, test-evolving.R,
#'   test-food-web.R, test-phylo-gradient.R, test-trait-gradient.R and test-ranges.R.
#' @srrstats {G5.4} Correctness tests compare results with fixed expected values: analytic
#'   equilibria of the simulation model (tests/testthat/test-analytic.R), values computed by
#'   hand (tests/testthat/test-trait-gradient.R, test-node.R) and values of other
#'   implementations (G5.4b). The extended tests compare with theoretical expectations over
#'   many simulated data sets (tests/testthat/test-extended-*.R).
#' @srrstats {G5.4a} The simulators, which have no reference with a shared random stream, are
#'   tested against trivial cases and analytic results: single-species and multispecies
#'   equilibria, competitive exclusion, source and sink, flow through the dispersal kernel
#'   (tests/testthat/test-analytic.R), and derivatives of the selection gradients
#'   (tests/testthat/test-food-web.R); further analytic results are in
#'   tests/testthat/test-extended-sims.R.
#' @srrstats {G5.4b} Tests call previous implementations when installed: picante
#'   (tests/testthat/test-summary-stats.R, test-phylo-beta.R), vegan (test-phylo-beta.R), cati
#'   (test-individuals.R, test-trait_variation.R) and nlme (test-trait-gradient.R). Stored
#'   outputs are used for the original MATLAB code of Leibold et al. (2010) and Peres-Neto et
#'   al. (2012), run in Octave (test-node.R, test-phylo-gradient.R), and for the package PCPS,
#'   version 1.0.9 (test-pcps.R).
#' @srrstats {G5.4c} The test "standardised values match the appendix" in
#'   tests/testthat/test-phylo-gradient.R uses the values printed in Appendix A of Peres-Neto
#'   et al. (2012).
#' @srrstats {G5.5} Tests with random numbers call `set.seed()` or pass `seed` first. The test
#'   files without a seed check deterministic code (`food_web_dynamics()`, the citation
#'   messages) or errors that do not depend on the values drawn.
#' @srrstats {G5.6} Parameter recovery from data simulated with known parameters is in the
#'   extended tests: rates of DAMOCLES, Pagel's lambda, parameters of the Ornstein-Uhlenbeck,
#'   rate-shift and Mk models, colonisation and extinction rates
#'   (tests/testthat/test-extended-evol.R); slopes of trait gradients and variance components
#'   (test-extended-traits.R); parameters of the simulators through ABC (test-extended-sims.R).
#' @srrstats {G5.6a} These tests pass within tolerances stated in a comment before each one
#'   and derived from theory or from Monte Carlo error: likelihood or binomial intervals, with
#'   `rejection_compatible()` of tests/testthat/helper-extended.R, or a stated relative error.
#' @srrstats {G5.6b} Each recovery test runs over several seeds or replicates (from 5 seeds to
#'   300 replicates, as its comment states) and sits in the extended suite because of its run
#'   time.
#' @srrstats {G5.7} Tests marked G5.7 in tests/testthat/test-extended-evol.R, -traits.R,
#'   -sims.R and -perm.R check that the error of the estimates falls as data grow (more
#'   communities, tips, sites, individuals or simulations), at the expected rate where theory
#'   gives one, and that the power of the permutation tests rises with the effect.
#' @srrstats {G5.8} See G5.8a to G5.8d.
#' @srrstats {G5.8a} Zero-length input stops with a message that names the argument, tested
#'   for example in tests/testthat/test-validation-traits.R (`traits = character(0)`,
#'   `rate_matrix(numeric(0))`), test-validation-inference.R (a table of scenarios without
#'   rows), test-maladaptation.R and test-env-pulse.R.
#' @srrstats {G5.8b} Input of unsupported type stops with a message: character vectors and
#'   columns in tests/testthat/test-neighbour_metrics.R, test-mc-data.R, test-surrogates.R,
#'   test-node.R and test-traits-external.R, and wrong types of each scalar argument in
#'   tests/testthat/test-validation-*.R and test-checks.R.
#' @srrstats {G5.8c} Missing values stop the function or give `NA`
#'   (tests/testthat/test-validation-analysis.R, test-surrogates.R, test-neighbour_metrics.R),
#'   and predictors or traits without variance stop with a message (test-checks.R,
#'   test-trait_plots.R, test-trait-gradient.R).
#' @srrstats {G5.8d} Data outside the scope of a method stop with a message: collinear
#'   predictors (tests/testthat/test-phylo-beta.R), fewer species than traits for
#'   Mahalanobis distances (test-assembly.R), too few individuals, species or sites
#'   (test-drift-test.R, test-trait-gradient.R, test-surrogates.R) and communities without
#'   individuals (test-pcc.R, test-pcps.R, test-phylo-beta.R).
#' @srrstats {G5.9} See G5.9a and G5.9b; the files tests/testthat/test-noise-perm.R,
#'   test-noise-traits.R, test-noise-evol.R and test-noise-sims.R run with the regular tests.
#' @srrstats {G5.9a} In these files a relative jitter of 1e-10 on predictors, branch lengths,
#'   trait values and abundances, with the same seed, leaves statistics, p-values and
#'   likelihoods unchanged within a stated tolerance.
#' @srrstats {G5.9b} In the same files, permutation and null-model tests run under several
#'   seeds give p-values that differ by no more than their Monte Carlo error.
#' @srrstats {G5.10} The extended tests (tests/testthat/test-extended-*.R) use testthat like
#'   the other tests and start with `skip_if_not_extended()` (tests/testthat/helper-extended.R),
#'   which skips them unless the environment variable `COMDYNAMITE_EXTENDED_TESTS` is `"true"`.
#'   The workflow .github/workflows/extended-tests.yaml sets it and runs them once a month and
#'   on demand.
#' @srrstats {G5.12} tests/README.md describes the extended tests: what each file checks, how
#'   to switch them on, their run time (several minutes), the need for cmdstanr and CmdStan in
#'   the Bayesian tests, and the workflow that runs them.
#' @noRd
NULL

#' NA_standards
#'
#' General standards that do not apply to this package, with the reason.
#'
#' @srrstatsNA {G2.14c} The package does not impute. Imputation of trait or community data
#'   depends on the phylogeny, the sampling design and the question, and belongs to dedicated
#'   software; the package help page tells users to remove or impute missing values before the
#'   analysis.
#' @srrstatsNA {G5.11} The extended tests simulate their data and need no large data sets or
#'   downloads.
#' @srrstatsNA {G5.11a} The extended tests download nothing.
#' @noRd
NULL
