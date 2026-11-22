#' srr_stats: Bayesian and Monte Carlo standards
#'
#' rOpenSci statistical standards (`srr`) of this category, with the parts of the
#' package that comply with each one or the reason it does not apply.
#'
#' The Bayesian parts of the package are (1) `trait_gradient_bayes()`
#' (`R/trait_gradient_bayes.R`, model in `inst/stan/trait_gradient.stan`), fitted by
#' MCMC in Stan through cmdstanr, called "the Stan model" below, and (2) approximate
#' Bayesian computation with random forests (`abc_table()`, `abc_fit()`,
#' `abc_recovery()`, `abc_choose()`, `abc_ppc()`, `abc_diagnose()` in `R/abc.R` and
#' `abc_cv()` in `R/abc_cv.R`, through abcrf), called "ABC" below. ABC has no chains,
#' starting values or convergence: standards on these apply to the Stan model only.
#'
#' @srrstats {BS1.1} Stan model: the data are the result of `trait_gradient()` (a table
#'   of individuals or populations with site, species and trait) and a `phylo` tree;
#'   the help page of `trait_gradient_bayes()` describes both and gives an example, and
#'   the article `vignettes/articles/trait-gradient.Rmd` shows the same steps. ABC: the
#'   help pages of `abc_table()` and `abc_fit()` describe the reference table (one row
#'   per simulation, parameters and statistics) and the observed statistics (a named
#'   vector, or a data frame with one data set per row), with examples of both, and
#'   `vignettes/abc.Rmd` works through a real data set.
#' @srrstats {BS1.2} Priors are described in text and code at the places listed in
#'   BS1.2b and BS1.2c. ABC takes uniform priors as a named list of lower and upper
#'   limits (`prior` of `abc_table()`). The Stan model has fixed, explicit priors on the
#'   standardised scale (normal(0, 1), exponential(1), lkj(2)) and a beta prior of the
#'   phylogenetic share that the user sets with `lambda_prior`.
#' @srrstatsTODO {BS1.2a} The README names `trait_gradient_bayes()` only in the
#'   installation notes and says nothing about ABC or about how priors are specified.
#'   Missing: a short paragraph or example in `README.Rmd` on the `prior` list of
#'   `abc_table()` and on `lambda_prior` of `trait_gradient_bayes()`.
#' @srrstats {BS1.2b} `vignettes/abc.Rmd` explains that parameters are drawn uniformly
#'   from the prior, defines the `prior` list of six parameters in code (with parameters
#'   that span orders of magnitude drawn on the log scale), passes it to `abc_table()`
#'   and sets the posterior draws against the prior ranges. For the Stan model the
#'   article `vignettes/articles/trait-gradient.Rmd` says only that the phylogenetic
#'   share has its own prior and refers to the help page; it has no example with
#'   `lambda_prior`.
#' @srrstats {BS1.2c} `abc_table()` documents `prior` and its description advises the
#'   log scale for parameters that span orders of magnitude; its example builds a prior
#'   list. `trait_gradient_bayes()` lists every prior of the model under "Priors" in
#'   Details, explains why the prior is put on the phylogenetic share and not on the
#'   two standard deviations, and documents `lambda_prior`; its example uses the
#'   default prior.
#' @srrstats {BS1.3} Stan model: `chains`, `parallel_chains`, `iter_warmup`,
#'   `iter_sampling`, `adapt_delta` and `seed` are documented with their valid ranges
#'   in `trait_gradient_bayes()`, and `...` passes any other argument of the
#'   `$sample()` method of cmdstanr (`max_treedepth`, `thin`, `refresh`, `init`). The
#'   convergence criteria are given in Details. ABC: `n` and `seed` of `abc_table()`,
#'   `ntree`, `ncores` and `min_node_size` of `abc_fit()`, `n_test`, `level` and `seed`
#'   of `abc_recovery()`, and `n`, `seed` of `abc_ppc()` and `abc_cv()` are documented
#'   in the respective help pages.
#' @srrstatsTODO {BS1.3a} Not documented. For the Stan model it works: `init` passed
#'   through `...` accepts the `fit` element of a previous result (checked by hand with
#'   cmdstanr 0.9.0.9002), but neither the help page nor an example says so. Missing: a
#'   sentence and an example line in `trait_gradient_bayes()`, such as
#'   `trait_gradient_bayes(tg, tree, init = previous$fit)`. Does not apply to ABC, which
#'   has no starting point.
#' @srrstats {BS1.4} The convergence checks of the Stan model (R-hat, bulk and tail
#'   effective sample sizes, divergent transitions) are computed after sampling in
#'   every call and returned in `diagnostics`; they do not change or stop the sampling,
#'   so there is no separate mode without them. Details of `trait_gradient_bayes()`
#'   gives the thresholds and the parameters they cover, and the example and
#'   `vignettes/articles/trait-gradient.Rmd` show `bayes$diagnostics` before the
#'   estimates are read. Does not apply to ABC.
#' @srrstats {BS1.5} `tests/testthat/test-trait-gradient-bayes.R` ("diagnostics flag
#'   poor convergence and divergences") changes one checker at a time (R-hat, bulk
#'   effective sample size, divergent transitions) and checks that each alone turns
#'   `ok` to `FALSE`, and that chains in different places raise R-hat above 1.1.
#' @srrstats {BS2.1} Stan model: `tg_bayes_data()` matches the species of the data
#'   with the tips of the tree, leaves out with a message the species that are not in
#'   the tree, prunes the tree to the species kept, stops with fewer than 4 species
#'   and builds the indices that link each observation to its species. ABC: `abc_fit()`
#'   stops when parameters or statistics are not columns of the table and drops
#'   incomplete simulations; `observed_statistics()` stops when the observed data lack
#'   a statistic of the fit; `abc_ppc()` stops when the table is not the one of the
#'   fit; `fold_summarise()` stops unless there is one fold per site.
#' @srrstats {BS2.1a} Tested in `tests/testthat/test-trait-gradient-bayes.R` ("the data
#'   of the model are standardised and indexed by species": numbers of observations,
#'   species index, message for species left out, error with fewer than 4 species), in
#'   `tests/testthat/test-abc.R` ("abc_fit checks its input", errors of `predict()` on
#'   missing statistics) and in `tests/testthat/test-abc-cv.R`.
#' @srrstats {BS2.2} `trait_gradient_bayes()` validates `lambda_prior` and `ou_alpha`
#'   at the top of the function, before the data are prepared and before the chains
#'   start. `abc_table()` validates `prior` with `check_prior()` before any parameter
#'   is drawn or simulation run.
#' @srrstats {BS2.3} `lambda_prior` must have exactly two values and each element of
#'   `prior` exactly two limits; any other length is an error and nothing is recycled
#'   or discarded. The error for a prior of length one is tested in
#'   `tests/testthat/test-abc.R` ("abc_table checks the prior").
#' @srrstats {BS2.4} The two values of `lambda_prior` are the two shape parameters of
#'   one beta prior, which the Stan model declares as a vector of length 2; the names
#'   of the `prior` list define the parameters of ABC, and `check_prior()` requires
#'   unique, non-empty names with one pair of limits each.
#' @srrstats {BS2.5} `lambda_prior` must be finite and greater than 0 (shape
#'   parameters of a beta distribution) and `ou_alpha` finite and greater than 0; the
#'   limits of each uniform prior of ABC must be finite and increasing. Tested in
#'   `tests/testthat/test-trait-gradient-bayes.R`, `tests/testthat/test-validation-analysis.R`
#'   and `tests/testthat/test-abc.R`.
#' @srrstats {BS2.6} `trait_gradient_bayes()` checks with `check_count()` and
#'   `check_number()` that `chains`, `parallel_chains` and `iter_sampling` are whole
#'   numbers of 1 or more, `iter_warmup` 0 or more, `adapt_delta` strictly between 0
#'   and 1 and `seed` a whole number of 0 or more. The ABC functions check `n`, `ntree`,
#'   `ncores`, `min_node_size`, `n_test`, `level` and `seed` in the same way. Tested in
#'   `tests/testthat/test-validation-analysis.R` and
#'   `tests/testthat/test-validation-inference.R`.
#' @srrstats {BS2.7} Stan model: starting values are controlled with `init`, passed
#'   through `...` to the `$sample()` method of cmdstanr, which takes one value or
#'   function for all chains or a list with one set of values per chain. The help page
#'   does not name `init` among the examples of `...` (see BS1.3a). Does not apply to
#'   ABC.
#' @srrstats {BS2.8} Stan model: `init` (through `...`) accepts the `fit` element of a
#'   previous result, so that the chains start from draws of the earlier posterior;
#'   checked by hand with cmdstanr 0.9.0.9002, not covered by a test and not documented
#'   (see BS1.3a). ABC: `abc_fit()` takes any data frame, so a reference table from an
#'   earlier run of `abc_table()` can be reused or extended.
#' @srrstats {BS2.9} Stan model: cmdstanr passes one seed and the chain identifiers 1
#'   to `chains` to CmdStan, which combines the two so that each chain has its own
#'   stream of random numbers; with `seed = NULL` cmdstanr draws the seed. ABC:
#'   `run_scenarios()`, called by `abc_table()`, sets a different seed (`seed` plus an
#'   offset) for every simulation, and `abc_cv()` for every fold and observation.
#' @srrstats {BS2.12} `trait_gradient_bayes()` has no verbosity argument of its own:
#'   `refresh`, `show_messages` and `show_exceptions` of the `$sample()` method of
#'   cmdstanr are passed through `...` (`refresh` is named in the help page), and by
#'   default the sampler prints progress and all messages. The ABC functions print
#'   nothing but a progress bar in interactive sessions.
#' @srrstatsTODO {BS2.13} Enabled but not tested. `refresh = 0` and
#'   `show_messages = FALSE` silence the progress and messages of the sampler while the
#'   warnings of `trait_gradient_bayes()` are kept (used in
#'   `tests/testthat/test-extended-traits.R`), but no test asserts that the output is
#'   silent and that the warnings remain. Missing: an extended test with
#'   `expect_silent()` or `capture.output()` on a converged fit, and one with
#'   `expect_warning()` on a fit that fails the diagnostics with the same arguments.
#' @srrstats {BS2.14} All warnings of `trait_gradient_bayes()` (weights, `ou_alpha`,
#'   failed summary, failed diagnostics) are raised with `warning()` and can be
#'   suppressed with `suppressWarnings()` or handled one by one;
#'   `tests/testthat/test-trait-gradient-bayes.R` does both ("trait_gradient_bayes()
#'   samples and builds the result"). The warnings of the sampler itself are turned off
#'   with `show_exceptions = FALSE` through `...`.
#' @srrstats {BS2.15} When the posterior cannot be summarised, `trait_gradient_bayes()`
#'   catches the error, turns it into a warning with the original message and returns
#'   the object with `signal`, `variances` and `species` set to `NULL` and the fit kept
#'   in `fit`; `print()` reports it. Tested in
#'   `tests/testthat/test-trait-gradient-bayes.R` ("could not summarise").
#' @srrstats {BS3.0} Stan model: `trait_gradient()`, which prepares the data, documents
#'   that trait values must be finite or missing and how missing values are treated;
#'   `trait_gradient_bayes()` documents (section Limits) that populations without a
#'   site mean enter only the trait response and that a species alone in all its sites
#'   gets `NA` for beta and alpha. ABC: `abc_fit()` documents that statistics with more
#'   than 5 percent of missing or infinite values are dropped and then the simulations
#'   with any missing statistic; observed statistics that are missing or infinite are
#'   an error (`observed_statistics()`).
#' @srrstats {BS3.1} Neither model has user-supplied predictors. The analogue in the
#'   Stan model is a singular phylogenetic correlation matrix (tips at distance 0),
#'   which `tg_bayes_data()` detects through the Cholesky factor and reports with an
#'   error before sampling, and an `ou_alpha` so large that all correlations are below
#'   0.01, which gives a warning that the phylogenetic and species effects cannot be
#'   told apart. Both are tested in `tests/testthat/test-trait-gradient-bayes.R`. In
#'   ABC, constant statistics are dropped by `abc_fit()` and `abc_choose()`, and
#'   correlated statistics are handled by the random forests.
#' @srrstats {BS4.0} `trait_gradient_bayes()` documents that the model is fitted in
#'   Stan with the `$sample()` method of cmdstanr, and where the Stan file is. The ABC
#'   functions cite the algorithms of abcrf: Raynal et al. (2019) and Meinshausen
#'   (2006) in `abc_fit()`, Pudlo et al. (2016) in `abc_choose()`.
#' @srrstats {BS4.2} Stan model: the `pp_check()` method for `mc_trait_gradient_bayes`
#'   (posterior predictive check of either response, with bayesplot) and the pointwise
#'   log-likelihoods for cross-validation with loo. ABC: `abc_recovery()` (recovery and
#'   coverage of the intervals on pseudo-observations), `abc_ppc()` (posterior
#'   predictive check), `abc_cv()` (cross-validation over folds of sites) and
#'   `abc_diagnose()` (effective sample size of the weights, statistics outside the
#'   reference table).
#' @srrstatsTODO {BS4.3} The checkers exist (R-hat, bulk and tail effective sample
#'   sizes from `posterior::summarise_draws()`, divergent transitions from the sampler;
#'   `tg_diagnostics_from()`), but the help page of `trait_gradient_bayes()` gives only
#'   the thresholds, with no reference for the statistics or their implementation.
#'   Missing: a `@references` entry for R-hat and the effective sample sizes as
#'   computed by the posterior package, and a mention of that package in Details.
#' @srrstats {BS4.5} The diagnostics count transitions that hit the maximum tree depth
#'   (`max_treedepth` in `diagnostics`) besides divergences, R-hat and effective sample
#'   sizes. When the diagnostics fail, `trait_gradient_bayes()` warns, returns
#'   `diagnostics$ok = FALSE` and `print()` marks the result "CHECK BEFORE READING";
#'   the fit is still returned for inspection. The help page says not to read such a
#'   fit and names longer chains (`iter_sampling`) as the first remedy; `adapt_delta`
#'   and `max_treedepth` are available for divergent transitions.
#' @srrstatsTODO {BS5.0} Partly met. The object of `trait_gradient_bayes()` holds the
#'   `CmdStanMCMC` fit, and `fit$metadata()$seed` gives the seed used, also when
#'   `seed = NULL`; this is not documented, and starting values are not kept unless
#'   given. The ABC objects do not keep their seeds: `abc_table()` returns the
#'   attributes `prior`, `parameters` and `statistics` only, and `abc_recovery()`,
#'   `abc_ppc()` and `abc_cv()` return no seed. Missing: a `seed` element in
#'   `mc_trait_gradient_bayes`, a `seed` attribute in `mc_abc_table` and in the results
#'   of the other three functions, and their documentation.
#' @srrstats {BS5.1} `mc_trait_gradient_bayes` returns `data` (the table of populations
#'   used, with species, site, trait, site mean and weight), `A` (the phylogenetic
#'   correlation matrix, with one row per species) and `ou_alpha`; `print()` gives the
#'   number of species. `mc_abc_fit` returns `training` (the simulations used),
#'   `parameters` and `statistics`, and `print()` gives their numbers; `mc_abc_table`
#'   carries the names of parameters and statistics as attributes.
#' @srrstats {BS5.2} `mc_trait_gradient_bayes` returns `lambda_prior`, the only prior
#'   that can be changed, and the same values inside the data of the fit; the fixed
#'   priors are in the help page and in `inst/stan/trait_gradient.stan`. `mc_abc_table`
#'   has the attribute `prior`, which `abc_fit()` copies to the element `prior` of
#'   `mc_abc_fit` (it is `NULL` when the table is a plain data frame).
#' @srrstats {BS5.3} `trait_gradient_bayes()` returns `diagnostics`: the largest R-hat,
#'   the smallest bulk and tail effective sample sizes, the number of divergent
#'   transitions and `ok`. `print()` shows them. For ABC the equivalent is
#'   `abc_diagnose()`, which gives the effective sample size and the Pareto k-hat of
#'   the weights behind a posterior.
#' @srrstats {BS5.4} All the checkers are always applied, and `diagnostics` has one
#'   named column for each (`max_rhat`, `min_ess_bulk`, `min_ess_tail`, `divergent`).
#' @srrstats {BS5.5} `diagnostics` is an element of the returned object and `ok` is
#'   `FALSE` when any check fails; the full diagnostics of the sampler are available
#'   from `fit` (for instance `fit$diagnostic_summary()`).
#' @srrstats {BS6.0} `print.mc_trait_gradient_bayes()`, `print.mc_abc_fit()` and
#'   `print.mc_abc_choice()`; `mc_abc_table` and `mc_abc_cv` are data frames and print
#'   as such. Tested in `tests/testthat/test-trait-gradient-bayes.R` and
#'   `tests/testthat/test-abc.R`.
#' @srrstatsTODO {BS6.1} There is no `plot()` method for `mc_trait_gradient_bayes`,
#'   `mc_abc_fit`, `mc_abc_choice` or `mc_abc_cv`. The only plot is the `pp_check()`
#'   method for `mc_trait_gradient_bayes`. Missing: `plot()` methods, for instance the
#'   posterior medians and intervals of `signal` and `species` for the Stan model and
#'   posterior against prior for `mc_abc_fit` (as drawn by hand in
#'   `vignettes/abc.Rmd`).
#' @srrstatsTODO {BS6.2} Applies to the Stan model only. Trace plots are neither
#'   provided nor documented; they can be drawn from the returned fit with bayesplot
#'   (`bayesplot::mcmc_trace(x$fit$draws("lambda"))`), and the warmup draws are not
#'   saved unless `save_warmup = TRUE` is passed through `...`. Missing: a documented
#'   example, or an option of a `plot()` method (BS6.1), that shows the chains with the
#'   warmup marked.
#' @srrstatsTODO {BS6.3} Partly met. `vignettes/abc.Rmd` plots histograms of the
#'   posterior draws of each ABC parameter against its prior range with ggplot2 code,
#'   but the package has no function for it, and nothing is documented for the Stan
#'   model apart from the posterior predictive densities of `pp_check()`. Missing: the
#'   `plot()` methods of BS6.1, or documented examples with bayesplot on `x$fit`.
#' @srrstats {BS6.4} `summary.mc_abc_cv()` averages the cross-validation over folds and
#'   data sets. The other objects have no `summary()` method: the tables `signal`,
#'   `variances` and `species` of `mc_trait_gradient_bayes` and the result of
#'   `predict()` for `mc_abc_fit` are already the posterior summaries.
#' @srrstats {BS7.0} Stan model: `tests/testthat/test-extended-traits.R`
#'   ("trait_gradient_bayes() recovers species means and converges") simulates data
#'   with known species means and requires the 95 percent intervals to cover them at
#'   least `qbinom(0.01, 24, 0.95)` times out of 24, with R-hat below 1.01 and no
#'   divergent transitions. ABC: `tests/testthat/test-extended-sims.R` ("abc_recovery
#'   recovers the parameters of a cheap simulator over seeds") draws the parameters
#'   from the prior and requires R2 above 0.9 and no bias beyond 4 standard errors, and
#'   `tests/testthat/test-abc.R` does a quick version in the regular suite. The
#'   extended tests run with `COMDYNAMITE_EXTENDED_TESTS=true`
#'   (`.github/workflows/extended-tests.yaml`, `tests/README.md`).
#' @srrstats {BS7.1} Stan model: `tests/testthat/test-extended-traits.R`
#'   ("trait_gradient_bayes() follows the prior of lambda when data say little") fits 5
#'   species under the priors beta(2, 8) and beta(8, 2) and requires each posterior
#'   median of lambda to lie inside the central 90 percent interval of its prior and
#'   the two medians to differ by more than 0.3. The data are nearly, not fully,
#'   uninformative: the model is not sampled from the prior alone. ABC: a parameter
#'   that no statistic informs keeps its prior, tested in `tests/testthat/test-abc.R`
#'   ("abc_fit recovers an informative parameter and not an uninformed one",
#'   "abc_recovery separates a recoverable from an unrecoverable parameter").
#' @srrstats {BS7.2} Stan model: the recovery test of BS7.0 also requires the median
#'   width of the intervals of the species means to be below a quarter of the standard
#'   deviation of the trait, far narrower than the prior. ABC:
#'   `tests/testthat/test-extended-sims.R` compares the coverage of the 90 percent
#'   intervals with the nominal level, including a toy model with a known posterior
#'   ("abc_recovery intervals are too narrow with one statistic and small leaves",
#'   "larger leaves bring the abc_recovery coverage towards the nominal level",
#'   "abc_recovery credible intervals are conservative with several statistics").
#' @srrstatsTODO {BS7.3} No test measures how computing time scales with the size of
#'   the input. `tests/testthat/test-extended-sims.R` ("the error of abc_fit falls as
#'   the reference table grows") tests the scaling of accuracy, not of time. Missing:
#'   extended tests that time `abc_fit()` over reference tables of increasing size and
#'   `trait_gradient_bayes()` over increasing numbers of species, and check that the
#'   time grows no faster than expected.
#' @srrstats {BS7.4} The Stan model standardises both responses and returns every
#'   estimate and the replicated data in trait units.
#'   `tests/testthat/test-trait-gradient-bayes.R` checks that the standardised data
#'   transform back to the input values exactly, and
#'   `tests/testthat/test-extended-traits.R` that the intervals of the species means,
#'   in trait units, contain the true means. For ABC, `tests/testthat/test-abc.R`
#'   checks that the estimate is close to the true parameter on its own scale and
#'   inside its interval.
#' @srrstats {BS7.4a} The Stan model assumes nothing about the scale of the trait
#'   because it centres and scales internally; the test "the data of the model are
#'   standardised and indexed by species" in
#'   `tests/testthat/test-trait-gradient-bayes.R` checks that the data sent to Stan
#'   have mean 0 and standard deviation 1 whatever the input, that `center` and `scale`
#'   recover the input, and that a trait that does not vary is refused. ABC makes no
#'   assumption on the scale of parameters or statistics (random forests are invariant
#'   to monotone transformations of the predictors).
#' @noRd
NULL

#' NA_standards
#'
#' Standards of this category that do not apply to this package, with the reasons.
#'
#' @srrstatsNA {BS1.0} The term "hyperparameter" is not used in the documentation,
#'   vignettes or messages of the package.
#' @srrstatsNA {BS1.3b} Each model has one algorithm. `trait_gradient_bayes()` calls
#'   only the `$sample()` method of cmdstanr (NUTS) and ABC uses only the random
#'   forests of abcrf, so there is no choice of sampling algorithm to document.
#' @srrstatsNA {BS2.10} `seed` must be one whole number (`check_count()` and
#'   `check_seed()` refuse a vector; tested for `abc_cv()` in
#'   `tests/testthat/test-validation-inference.R`), so identical seeds cannot be passed
#'   to distinct chains; the chains are told apart by their identifiers (BS2.9).
#' @srrstatsNA {BS2.11} The package defines no argument for starting values; the only
#'   route is `init` of cmdstanr, passed through `...`.
#' @srrstatsNA {BS3.2} There are no predictors, hence no collinear data to process.
#'   The singular correlation matrix of BS3.1 stops with an error before sampling.
#' @srrstatsNA {BS4.1} The package implements no sampler; it uses those of Stan and
#'   abcrf.
#' @srrstatsNA {BS4.4} The sampler of Stan runs chains of a fixed length (`iter_warmup`
#'   and `iter_sampling`) and has no stopping rule on convergence; the checks are made
#'   after sampling. ABC with random forests is not iterative.
#' @srrstatsNA {BS4.6} Sampling is not stopped by a convergence checker (BS4.4), so
#'   there are no two modes to compare.
#' @srrstatsNA {BS4.7} The thresholds of the convergence checks are fixed (R-hat below
#'   1.01, effective sample sizes above 400, no divergent transitions) and are not
#'   arguments.
#' @srrstatsNA {BS6.5} Optional and not provided; the returned `CmdStanMCMC` fit can be
#'   passed to bayesplot for such a graphic.
#' @noRd
NULL
