#' srr_stats: Regression and supervised learning standards
#'
#' rOpenSci statistical standards (`srr`) of this category, with the parts of the
#' package that comply with each one or the reason it does not apply.
#'
#' @srrstats {RE1.0} `trait_flex()` (R/trait_variation.R) takes its model as a one-sided
#'   `formula` over the variables of `env`. The other regressions take matrices
#'   (`node_test()`, `pcc_test()`, `phylo_beta_test()`, `pcps_test()`, `phylo_vector()`,
#'   `trait_partition()`) because their response is not a column of a data frame: it is
#'   computed from the community matrix and the tree (node contrasts, phylogenetic
#'   eigenvectors, a phylogenetic beta diversity matrix), and every predictor of `env` enters
#'   additively. `trait_gradient()`, `population_shifts()` and `drift_test()` fit one fixed
#'   model and take column names. `abc_fit()` and `abc_choose()` regress each parameter (or
#'   the model index) on all statistics of the reference table, chosen with `parameters` and
#'   `statistics`.
#' @srrstats {RE1.1} In `trait_flex()` the right-hand side of `formula` is joined to each
#'   response with `stats::reformulate()` and fitted with `stats::lm()` on the rows of `env`
#'   matched to the sites by name, so the design matrix is that of `stats::model.matrix()`
#'   with the default contrasts. The help of `trait_flex()` states that `formula` is
#'   one-sided, must keep the intercept and uses the variables of `env`.
#' @srrstats {RE1.2} The help pages give the format of the predictors: `env` is a vector,
#'   matrix or data frame with one row per site, matched to `comm` by row names, and the `env`
#'   entries of `node_test()`, `pcc_test()`, `phylo_beta_test()`, `pcps_test()`, `phylo_vector()`
#'   and `trait_partition()` say that columns must be numeric or logical and that factors and
#'   character columns stop with an error, to be coded as dummy variables; `check_env()` and
#'   `check_numeric_columns()` (R/checks.R) enforce it.
#' @srrstats {RE1.3} Names are kept in the results: node numbers and predictor names in the
#'   tables and matrices of `node_test()`, predictor and trait names in `trait_partition()`,
#'   species names in `phylo_vector()` and `phylo_eigenvectors()`, site and species names in the
#'   `scores` and `species` of `pcps_test()`, site, species and trait names in `trait_flex()` and
#'   `trait_gradient()`, trait names in the `W` of `drift_test()`. The `scores` of
#'   `phylo_beta_test()` have the site names of `comm` as row names (`predicted_scores()` in
#'   R/phylo_beta.R), and `predict()` of `abc_fit()` adds a column `name` with the row names of
#'   `newdata` (tests in tests/testthat/test-phylo-beta.R and tests/testthat/test-abc.R).
#' @srrstats {RE1.3a} Where row names cannot be kept as row names, because a result stacks
#'   several rows per observation, they are kept as a column: `name` in `predict()` of
#'   `abc_fit()`, beside the row index `observation`, as its help page says.
#' @srrstats {RE1.4} The help pages state the assumptions and what follows when they fail,
#'   with rates from simulations: exchangeable sites in the permutation tests (`node_test()`,
#'   `pcps_test()`: free permutations reject too often with spatial autocorrelation, and the
#'   restricted null models that replace them); independence of species turnover and
#'   predictors in the F test of `trait_flex()` (conservative otherwise); proportional G and
#'   P matrices and independent components in `drift_test()`; the sources of bias of
#'   `population_shifts()`; approximate intervals in `abc_fit()`. The violations are tested
#'   in tests/testthat/test-extended-perm.R (autocorrelated fields with free, torus,
#'   variogram, patch and rigid null models), test-extended-traits.R (turnover along the
#'   predictor in `trait_flex()`) and test-extended-sims.R (coverage of `abc_recovery()` with
#'   one statistic).
#' @srrstats {RE2.0} The transformations are those of the published methods and are
#'   described in the help pages: `node_test()` and `pcc_test()` standardise contrasts and
#'   eigenvectors with species weights and predictors with site weights; `phylo_beta_test()`
#'   and `pcps_test()` centre the predictors; `population_shifts()` centres each variable
#'   within species; `drift_test()` takes logarithms of variances and eigenvalues;
#'   `trait_gradient()` fills missing trait values with species means; `abc_choose()` adds
#'   linear discriminant axes (avoided with `lda = FALSE`). The standardisations cannot be
#'   switched off because the statistics are defined on the standardised variables and do
#'   not change with the scale of the predictors (tests/testthat/test-phylo-beta.R). No
#'   label is converted to a factor: non-numeric predictors stop with an error
#'   (`check_env()`), and `trait_flex()` passes factors of `env` to `stats::lm()` unchanged.
#' @srrstats {RE2.1} Missing and infinite values are told apart: `check_comm()` and
#'   `check_env()` (R/checks.R) stop on `NA` and on `Inf` with separate messages; `abc_fit()`
#'   drops incomplete statistics and simulations (documented); `trait_flex()`, `drift_test()` and
#'   `population_shifts()` take `na_action = c("omit", "fail")`, where "omit" leaves out
#'   incomplete rows with a message giving their number and "fail" stops (tests in
#'   tests/testthat/test-checks.R and tests/testthat/test-population_shifts.R).
#' @srrstats {RE2.2} Missing values in responses and predictors are handled apart in
#'   `trait_flex()`: individuals without a trait value leave the means with a message, sites
#'   without predictors leave the model with a warning, and `na_action = "fail"` stops on either,
#'   naming the sites with missing predictors. In `population_shifts()` and `drift_test()` a row
#'   needs both response and predictors, so they share one rule.
#' @srrstats {RE2.4} `check_env()` in R/checks.R is the pre-processing routine: it stops on
#'   predictors without variance, found with `no_variance()`, which allows for
#'   floating-point noise, and on predictors that are linear combinations of the others
#'   (rank of the centred matrix). The cases it does not cover are listed in RE2.4a and
#'   RE2.4b.
#' @srrstats {RE2.4a} `check_env()` stops `node_test()`, `pcc_test()`, `phylo_beta_test()`,
#'   `pcps_test()` and `phylo_vector()` on collinear predictors; `trait_partition()` accepts
#'   them, returns `NA` semipartial correlations and says so in its help page.
#'   `check_collinear()` (R/checks.R) is used by `trait_flex()` (on the model matrix of
#'   `formula`) and `population_shifts()` (on the site values of `env`): it stops when a
#'   predictor is a linear combination of the others, naming it, and warns, without removing any
#'   variable, when a variance inflation factor exceeds 10 (tests in tests/testthat/test-checks.R
#'   and tests/testthat/test-population_shifts.R).
#' @srrstats {RE2.4b} `trait_flex()`, `drift_test()`, `node_test()` and `pcps_test()` warn,
#'   through `warn_perfect_fit()` (R/checks.R), when a response is an exact linear function
#'   of the predictors (a residual sum of squares of zero up to rounding error), naming the
#'   models concerned and saying that their F statistics, standard errors and p-values are
#'   not meaningful.
#' @srrstats {RE3.0} Closed-form fits (`stats::lm()`, `stats::lm.fit()`, QR) need no
#'   convergence check, and random forests do not iterate to convergence. Of the iterative fits,
#'   `population_shifts()` warns when the width of the kernel stays at an end of its grid or the
#'   residual sum of squares is flat along it, and `trait_gradient_phylo()` starts the OU fit
#'   from three points, marks estimates at a bound in `at_bound` and warns when `stats::optim()`
#'   returns a convergence code other than 0 (test in tests/testthat/test-trait-gradient.R).
#' @srrstats {RE3.1} The warnings can be silenced with `suppressWarnings()`, and the
#'   information stays in the result: `population_shifts()` keeps the attribute `profile`, with
#'   the residual sum of squares and the estimates at each width, and `trait_gradient_phylo()`
#'   keeps the columns `at_bound` and `converged` of its `models` table.
#' @srrstats {RE3.2} The tolerances are documented: the help of `trait_gradient_phylo()` gives
#'   the bounds of the OU parameter, the use of `stats::optimize()` with its default tolerance
#'   and of `stats::optim()` with L-BFGS-B and its default controls from three starts; that of
#'   `population_shifts()` gives the grid of widths, its refinement by `optimize()` with its
#'   default tolerance, the 5% threshold of a flat profile, and the relative tolerance of 1e-8
#'   and 50 iterations of the variance of populations.
#' @srrstats {RE3.3} `population_shifts()` lets the user fix `width`, which removes the search.
#'   The other tolerances are fixed, with the reason given in the help pages: the surfaces have
#'   one or two parameters on bounded intervals, several starts guard against local optima better
#'   than a tighter tolerance would, and a finer width or variance would not change the estimates
#'   beyond their sampling error; non-convergence is reported (`converged`, warnings).
#' @srrstats {RE4.0} Each analysis returns an object of its own S3 class: `mc_node_test`,
#'   `mc_pcc`, `mc_phylo_beta`, `mc_pcps`, `mc_phylo_vector`, `mc_trait_partition`,
#'   `mc_trait_gradient`, `mc_trait_gradient_test`, `mc_trait_gradient_phylo`, `trait_flex`,
#'   `population_shifts` (a data frame with attributes), `drift_test`, and, for the random
#'   forests, `mc_abc_fit` and `mc_abc_choice`, which hold the `regAbcrf` and `abcrf`
#'   objects of the package abcrf together with the training table.
#' @srrstats {RE4.2} `coef()` methods (R/accessors.R) for `drift_test` (intercept and slope),
#'   `population_shifts` (slopes and displacement), `mc_trait_gradient` (beta, alpha, intercept
#'   and slope per species) and `mc_trait_gradient_phylo` (mu, sigma2 and a per response and
#'   model); `mc_damocles_fit` has its estimates in `confint()`. Tests in tests/testthat/test-
#'   accessors.R, test-population_shifts.R and test-trait-gradient.R.
#' @srrstats {RE4.3} `confint()` methods for `drift_test` (from the \eqn{t} distribution of its
#'   regression, at any `level`), `population_shifts` (the basic bootstrap intervals, with
#'   `nboot > 0`, at the level they were computed) and `mc_damocles_fit`.
#' @srrstats {RE4.4} The only model chosen by the user is the `formula` of `trait_flex()`; the
#'   result keeps it as the element `formula`, returned by `formula()` (R/accessors.R) and shown
#'   by `print()`.
#' @srrstats {RE4.5} `nobs()` methods for `drift_test` (the number of traits, the points of its
#'   regression), `population_shifts` (populations that entered as responses), `trait_flex`
#'   (sites in the model), `mc_trait_gradient` (populations), `mc_trait_gradient_phylo`
#'   (species), `mc_node_test` (sites at the root) and `mc_abc_fit` (simulations in the training
#'   table).
#' @srrstats {RE4.6} `vcov()` methods for `drift_test` (from its `stats::lm()` fit, now kept as
#'   `fit`) and for `population_shifts` with `nboot > 0` (the covariance of the coefficients
#'   among the bootstrap simulations from the fitted model). The permutation tests estimate no
#'   parametric covariance.
#' @srrstats {RE4.7} Where a fit is iterative the result carries what is needed to judge
#'   it: the attribute `profile` of `population_shifts()` (residual sum of squares, R2 and
#'   estimates at each width), `at_bound` and `loglik` in `trait_gradient_phylo()`, the
#'   out-of-bag `error` and `confusion` of `abc_choose()`, and the effective sample size
#'   and Pareto k-hat of the forest weights from `abc_diagnose()`. The convergence code of
#'   `stats::optim()` is still to be added (RE3.0).
#' @srrstats {RE4.8} The responses that the functions compute from the data are returned:
#'   `sites` in `trait_flex()` (the specific, fixed and intraspecific site means), `scores`
#'   in `pcps_test()`, `components` in `drift_test()` (variance among groups per component),
#'   `populations` and `species` in `trait_gradient()`, `data` in `trait_gradient_phylo()`
#'   (values with their standard errors), `training` in `abc_fit()` and `abc_choose()`.
#'   `node_test()`, `pcc_test()` and `phylo_beta_test()` return the statistics of each node,
#'   eigenvector or matrix and not the site-level responses, which `node_contrasts()`,
#'   `phylo_eigenvectors()` and `phylo_beta()` give.
#' @srrstats {RE4.9} `fitted()` methods for `drift_test` (one value per principal component),
#'   `trait_flex` (one per site, trait and component, in the element `model`) and
#'   `population_shifts` (one per half-sample mean of a population, in the attribute `model`);
#'   `predict()` gives modelled values for `mc_abc_fit` and `mc_abc_choice`.
#' @srrstats {RE4.10} `residuals()` methods for the same three classes as RE4.9; their scales
#'   are stated in the help pages: `log(between)` for `drift_test`, trait units for `trait_flex`,
#'   and trait units centred within species for `population_shifts`.
#' @srrstats {RE4.11} Goodness of fit and effect sizes are part of every result: `r2`,
#'   `r2_adjusted` (both weighted as the fit), `r2_original`, `r2_adjusted_original` and the
#'   semipartial correlation of each predictor in `node_test()`; `r2_mean` and `r2_adjusted` in `pcc_test()`; the adjusted R2 and
#'   pseudo-F in `phylo_beta_test()`; `r2` and F per axis in `pcps_test()`; the fractions of
#'   `trait_partition()`; explained and residual sums of squares in `trait_flex()`;
#'   `r_squared` in `population_shifts()`; log-likelihood and AIC in
#'   `trait_gradient_phylo()`; standardised effect sizes in `trait_gradient_test()`; and for
#'   the random forests the R2, bias and coverage of `abc_recovery()`, the coverage and
#'   standardised errors of `abc_cv()` and the out-of-bag error of `abc_choose()`.
#' @srrstats {RE4.13} Predictors are returned where the function changes or subsets them:
#'   `env` in `pcps_test()` (the predictors of the occupied sites), `eigenvalue` in the
#'   `components` of `drift_test()`, `site_mean` in the `populations` of `trait_gradient()`,
#'   and the statistics in the `training` of `abc_fit()` and `abc_choose()`, with the names
#'   of those kept after constant and incomplete ones are dropped in `statistics`. Results
#'   of `node_test()` and `trait_partition()` carry the predictor names as dimnames, and
#'   `permutations` records the null model of sites that was applied to them.
#' @srrstats {RE4.14} `predict()` for `mc_abc_fit` returns, for each new observation, the
#'   posterior median with a credible interval (`low`, `high`) from quantile regression
#'   forests, and for `mc_abc_choice` the posterior probability of the chosen model.
#'   `abc_diagnose()` lists the statistics of an observation that fall outside the reference
#'   table, where the forests extrapolate, and `abc_recovery()` measures the error and
#'   coverage of the predictions on simulations left out.
#' @srrstats {RE4.17} `print()` methods exist for `mc_node_test`, `mc_pcc`, `mc_pcps`,
#'   `mc_phylo_beta`, `mc_phylo_vector`, `mc_trait_partition`, `mc_trait_gradient`,
#'   `mc_trait_gradient_test`, `mc_trait_gradient_phylo`, `drift_test`, `trait_flex` (new, with
#'   the model and the decomposition), `mc_abc_fit` and `mc_abc_choice`, with the numbers of
#'   sites, species or simulations, the statistics and the permutations used; `population_shifts`
#'   prints as the data frame it is.
#' @srrstats {RE4.18} `summary()` methods exist for `mc_node_test` (nodes in order of
#'   `p_adjusted`), `mc_pcc`, `mc_pcps`, `mc_phylo_beta`, `mc_phylo_vector` and
#'   `mc_trait_partition` (help topic `result-methods`, R/methods.R) and for `mc_abc_cv`
#'   (R/abc_cv.R), together with `as.data.frame()` methods. Costly summaries are separate
#'   functions and not `summary()` methods: the bootstrap of `population_shifts()` runs only
#'   with `nboot > 0`, and `abc_recovery()`, `abc_cv()` and `abc_diagnose()` are apart from
#'   `abc_fit()`.
#' @srrstats {RE5.0} `data-raw/benchmark_regression.R` times `node_test()`, `pcc_test()`,
#'   `pcps_test()`, `phylo_beta_test()` and `trait_partition()` over 25 to 200 sites, 25 to
#'   100 species and 99 or 499 permutations; the README (section "Large analyses") gives the
#'   exponents of time on each of these and the time of the largest case.
#' @srrstats {RE6.0} `plot()` methods exist for every regression class: `mc_node_test`,
#'   `mc_pcps`, `mc_phylo_beta`, `mc_trait_partition`, `mc_trait_gradient`,
#'   `mc_trait_gradient_test`, `trait_flex`, `population_shifts` and `drift_test`, and now
#'   `mc_pcc` (adjusted R2 of each eigenvector), `mc_phylo_vector` (correlation of each
#'   eigenvector with the predictor, retained ones marked), `mc_trait_gradient_phylo` (Akaike
#'   weights of the models), `mc_abc_fit` (posterior against prior) and `mc_abc_choice`
#'   (confusion matrix, or votes for an observation). Tests in the test file of each function.
#' @srrstats {RE6.1} Where the plot is drawn by a named function, an S3 method sends the
#'   user to it: `plot.mc_node_test()` calls `plot_node_test()`, `plot.mc_trait_partition()`
#'   calls `plot_trait_partition()` and `plot.mc_phylo_beta()` calls `plot_phylo_beta()`
#'   (R/methods.R, help topic `result-methods`).
#' @srrstats {RE6.2} The plots of the fitted models show the fit: `plot.drift_test()` draws
#'   the points, the fitted line and the line of slope 1 expected under drift;
#'   `plot.mc_trait_gradient()` the populations with the fitted line of each species over
#'   the range of its sites; `plot.population_shifts()` the profile of R2 and of the
#'   estimates along the width, with the bootstrap interval of the displacement;
#'   `plot_phylo_beta()` the community ordered by the scores of the predicted matrix;
#'   `plot.mc_pcps()` the ordination with the predictor. The results of `node_test()`,
#'   `trait_partition()` and `trait_gradient_test()` are tests, so their plots show the
#'   tested statistics (significant nodes on the tree, fractions, null distributions with
#'   the observed value) and not fitted values.
#' @srrstats {RE7.0} Exact relationships among predictors are tested:
#'   tests/testthat/test-phylo-beta.R (a predictor that is twice another),
#'   tests/testthat/test-node.R (collinear predictors in `node_fit()` and `weighted_fit()`,
#'   and a predictor without variance among the sites of a node),
#'   tests/testthat/test-phylo-gradient.R (the same variable in `env` and `space` of
#'   `trait_partition()`) and tests/testthat/test-checks.R (a constant predictor with
#'   floating-point noise).
#' @srrstats {RE7.0a} The same tests confirm the rejection: `phylo_beta_test()` stops with
#'   "collinear", `check_env()` stops with "without variance", `node_fit()` returns `NA`
#'   and `weighted_fit()` `NULL` (so the node gets `NA`), and `trait_partition()` returns
#'   `NA` semipartial correlations.
#' @srrstats {RE7.1} `tests/testthat/test-perfect-fit.R` fits noiseless responses with
#'   `trait_flex()`, `drift_test()`, `node_test()` and `pcps_test()`, checks the warning of
#'   RE2.4b and that the exact relation is recovered, and checks that noisy versions of the
#'   same data raise no warning.
#' @srrstats {RE7.1a} `tests/testthat/test-perfect-fit.R` checks that a fit of
#'   `trait_flex()` on noiseless data takes no more than three times as long as on noisy
#'   data, plus half a second.
#' @srrstats {RE7.2} Tests show that names of the input reach the output:
#'   tests/testthat/test-node.R (predictor and event names as column names of
#'   `semipartial`, defaults for unnamed predictors), tests/testthat/test-phylo-gradient.R
#'   (species names of `pev`, predictor and trait names of the semipartial matrix),
#'   tests/testthat/test-pcc.R (tip labels as row names of the eigenvectors),
#'   tests/testthat/test-trait_variation.R (sites of `env` matched by name in
#'   `trait_flex()`) and tests/testthat/test-drift-test.R (trait names of `W`).
#' @srrstats {RE7.3} The methods are tested: `as.data.frame()` and `summary()` in
#'   tests/testthat/test-node.R and test-pcps.R, `print()` and `plot()` in the test file of each
#'   function, and the accessors `coef()`, `confint()`, `vcov()`, `nobs()`, `fitted()`,
#'   `residuals()` and `formula()` in tests/testthat/test-accessors.R, test-population_shifts.R,
#'   test-trait-gradient.R, test-node.R and test-abc.R.
#' @noRd
NULL

#' NA_standards
#'
#' Standards of this category that do not apply to this package, with the reasons.
#'
#' @srrstatsNA {RE2.3} Centring and offsets are fixed by each method and a parameter would
#'   not change any result: the permutation tests standardise or centre the predictors
#'   internally and are invariant to their location and scale
#'   (tests/testthat/test-phylo-beta.R), `trait_flex()` requires the intercept because its
#'   sums of squares are taken around the mean, and `population_shifts()` removes the
#'   species constants by centring within species.
#' @srrstatsNA {RE4.1} Optional, and not offered: the regressions have a closed form and
#'   cost little before the permutations start. For ABC the costly step, the simulation of
#'   the reference table, is already apart from the fit, in `abc_table()`.
#' @srrstatsNA {RE4.12} No function applies a transformation chosen by the user that would
#'   need an inverse. The fixed transformations (standardisation with weights, logarithms in
#'   `drift_test()`) are internal, and results such as the `components` of `drift_test()`
#'   are returned on the original scale.
#' @srrstatsNA {RE4.15} No function forecasts along a horizon (there is no time-series or
#'   extrapolating model). The only `predict()` methods are those of the ABC random
#'   forests, which predict within the range of the reference table; the help of
#'   `abc_diagnose()` states that forests do not extrapolate beyond it.
#' @srrstatsNA {RE4.16} No model with a `predict()` method has categorical groups among its
#'   predictors: the predictors of the ABC forests are numeric summary statistics. In
#'   `abc_choose()` the categories are the response (the candidate models), and a model
#'   absent from the reference tables cannot be predicted.
#' @srrstatsNA {RE6.3} No `plot()` method draws predictions: `predict()` exists only for
#'   the ABC random forests and returns a table, and no model forecasts (see RE4.15).
#' @srrstatsNA {RE7.4} No function forecasts along a horizon (see RE4.15), so there are no
#'   forecast errors to test.
#' @noRd
NULL
