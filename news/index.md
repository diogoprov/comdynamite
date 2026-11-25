# Changelog

## comdynamite (development version)

- [`population_stats()`](https://provetelab.org/comdynamite/reference/population_stats.md)
  gains `temporal_step` and `temporal_step_slope`: the squared change of
  population means between consecutive surveys, corrected for sampling
  and divided by the interval and the variance within populations, which
  estimates the variance of drift per generation relative to the
  phenotypic variance, and its slope on the inverse of abundance. In
  simulations (`data-raw/abc_temporal_drift.R`), ten consecutive surveys
  raised the R² of Ne/N by 0.07 at most over two surveys 50 steps apart
  (0.39 to 0.56), because drift per generation also depends on G / Ne;
  the statistic tracked G / (Ne (G + V_E)) closely (Spearman 0.83).

- `data-raw/abc_g_ne.R` reparametrises the genetic variance and the
  effective size in the ABC of `data-raw/abc_identifiability.R`: the
  intensity of drift, G / (Ne/N), and the response to selection, G /
  (sigma^2 + G + V_E), were recovered with R² 0.89 to 0.96 in every
  case, against 0.29 to 0.61 for Ne/N, but Ne/N implied by the estimates
  of G and G / (Ne/N) was never better than its own forest. Abundances
  and trait means therefore seem to inform the strength of drift rather
  than the effective size.

- The regional pools of `anurans` can now be cited through a GBIF
  download with a DOI (10.15468/dl.g4gz8u), made with the queries of 11
  November 2026; its counts by species differ slightly from those used,
  as described in `data-raw/anurans_pool/README.md`.

- `data-raw/abc_variable_design.R` tests ABC recovery when the tree and
  the landscape vary among simulations, with and without imperfect
  detection. Emigration, competition, niche breadth and detection were
  still recovered (R² 0.82 to 0.99), while the phylogenetic signal of
  the niche (`lambda`) fell from 0.43 to 0.53 with one tree to 0.32, and
  descriptors of the tree and of the landscape added to the statistics
  did not restore it.
  [`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
  now says so.

- `data-raw/abc_cv.R` repeats the study of
  [`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md)
  with spatial blocks against random folds in grids of 25 and 100 sites,
  with and without a torus. Blocks rejected a wrong model (species
  tracking a hidden, autocorrelated environment) more often than random
  folds only in the grid of 100 sites, by 1.9 percentage points without
  a torus and 1.1 with one, below the 3 points set as the criterion for
  recommending them; part of this gain came from blocks being less
  conservative also when the model was true.
  [`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
  now says so.

- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  accepts `seed_mean` per species, for colonisation from a region where
  some species are commoner than others.
  [`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
  uses it to test limited colonisation as an explanation of the
  competition chosen for the anurans (`data-raw/abc_colonisation.R`):
  without competition, limited colonisation was read as competition in
  half of the simulated data sets, and in up to 89% when regional
  abundances were most uneven, but for the observed ponds competition
  was still chosen, with a lower estimate and a milder misfit of the
  standardised MPD.

- [`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
  reports a test of the hypothesis it raised for the misfit of the
  phylogenetic statistics of the anurans
  (`data-raw/abc_conserved_filter.R`): a filter of tolerance to urban
  cover conserved in the phylogeny did not remove the misfit of the
  standardised MPD, and model choice still favoured competition. In the
  same simulations, the posterior predictive check flagged the MPD in
  23% to 39% of the data sets from a conserved filter without
  competition, against 4% when the model was true.

- [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
  gains `niche_k` and `niche_dist_cor`, which need `env` and `tree` and
  use only what can be observed: Blomberg’s K of the realised niche
  positions of the species present (the mean environment of the sites
  each occupies, weighted by abundance) and the correlation of
  phylogenetic distance with the difference of these positions (the
  Mantel statistic, without the test). In simulations
  (`data-raw/abc_lambda.R`), they raised the recovery of the
  phylogenetic signal of the niche (`lambda`) by ABC from R² 0 to 0.26
  with 30 species and to 0.34 with 100, close to what Blomberg’s K of
  the simulated optima gave (0.40 and 0.37), but not to the criterion of
  R² 0.5; with incidence alone, to 0.05 and 0.17.

- New
  [`population_stats()`](https://provetelab.org/comdynamite/reference/population_stats.md):
  statistics of trait variation among the populations of each species,
  from measured individuals
  ([`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  output or field data), for ABC. Besides the regression of population
  means on the environment, it gives the variance within populations, a
  QST-like divergence beyond the environment, semivariances between near
  and far populations, the relation of residuals to population size and,
  with two surveys, the change of the means between them (the basis of
  the temporal method for the effective size). In simulations
  (`data-raw/abc_identifiability.R`), the variance within populations
  made the environmental variance recoverable (R² 0.72 to 0.96, against
  0.15 to 0.45 from the means alone), and the new statistics raised the
  R² of Ne/N from 0.08 to 0.33 to 0.34 to 0.61, also with unknown
  ancestral optima; the number of sampled populations mattered more than
  the number of individuals per population.
  [`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
  gains a short guide to sampling designs from these results.

- [`damocles_samples()`](https://provetelab.org/comdynamite/reference/damocles_samples.md)
  gains `integrate`: with `integrate = TRUE` it also fits the parameters
  that maximise the likelihood averaged over the configurations
  (posterior samples of regimes or branch rates) and tests them against
  constant rates by a likelihood-ratio test, returned as `$integrated`
  and shown by [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html). In simulations
  with ClaDS posteriors (`data-raw/s3_clads_validation.R`) this was the
  most powerful way to test a branch-rate covariate from posterior
  samples, with the size of the test kept.

- [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md):
  `loci` gives a genetic architecture of a few diploid loci of large,
  additive effect. Each population carries its allele frequencies, its
  genetic mean and variance follow from them, and selection changes them
  by Wright’s (1937) formula, which also removes variance near the
  optimum until mutation balances it. Drift, gene flow, founders, seeds
  and lottery recruitment act on the frequencies, which are returned as
  `allele_frequency`. Tested in the article “Replicating Leibold et
  al. (2019)”: at their mutation rate, the loci give species sorting,
  and no setting gives their local dominance.

- `sample_individuals(source = "evolution")` and
  [`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md)
  now use the genetic variance of each population at the last step when
  the simulation has `mutation_var` (or `loci`). Before, a partial match
  of `genetic_var` read that matrix as one value per species, so the
  variance of the sampled individuals was wrong.

- [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md):
  `variance_selection = "sexual"` keeps the reduction of the genetic
  variance by stabilising selection as linkage disequilibrium, which
  recombination halves each generation (Bulmer 1971), so the variance
  recovers when selection relaxes. `TRUE` (now also `"clonal"`) keeps
  the permanent reduction of clonal inheritance. The disequilibrium of
  each population is returned as `linkage_disequilibrium`. The article
  “Replicating Leibold et al. (2019)” tests it: the transient reduction
  behaves almost like no erosion and does not produce their local
  dominance.

- Internal: the long and complex functions that goodpractice flagged
  (the simulation engines, trait simulation, landscapes, DAMOCLES, the
  trait-variation analyses,
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
  ABC and several plots) are split into smaller internal helpers, each
  at most 50 lines with a cyclomatic complexity of at most 15. Results,
  messages and the sequence of random draws are unchanged: every
  function gives identical output for the same seed. The tests follow
  the lintr style that goodpractice checks.

## comdynamite 0.11.0

- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `r2`, `r2_adjusted` and the semipartial correlations now sum squares
  with the site weights of the fit, so they stay between 0 and 1 (and -1
  and 1). Before, they followed the original code, which sums squares
  without weights and can give an R2 above

  1.  The original values are kept as `r2_original` and
      `r2_adjusted_original`; `r2_weighted` now equals `r2` and will be
      removed in version 1.0.0. `statistic` and p-values do not change.

- [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
  [`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  warn when a response is an exact linear function of the predictors,
  since their F statistics, standard errors and p-values are then not
  meaningful.
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md)
  no longer relays the warning of
  [`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html) in that
  case.

- Error, warning and message texts are now unique; where the same check
  is made in several functions, the text names the function.

- The column `reps` of
  [`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md)
  is an integer.

- The README gives how the time of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  grows with sites, species and permutations (script in
  `data-raw/benchmark_regression.R`).

- Text in plots is larger by default (3.2 to 3.5 mm), and
  [`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md)
  and
  [`plot_assembly()`](https://provetelab.org/comdynamite/reference/plot_assembly.md)
  gain `text_size`.

- Documentation: the intended audience and the priors of the Bayesian
  functions in the README; how to restart
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md)
  from an earlier fit (`init = previous$fit`) and the reference for its
  convergence diagnostics; the tolerances of the optimisers of
  [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
  and
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md);
  that predictors must be numeric or logical; and, in
  [`?comdynamite`](https://provetelab.org/comdynamite/reference/comdynamite-package.md),
  how factors, units and integer tables are converted.

- New accessor methods for the regression results:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) for
  [`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md)
  (which now keeps its [`stats::lm()`](https://rdrr.io/r/stats/lm.html)
  fit as `fit`) and
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  (intervals and covariance from the bootstrap; fitted values and
  residuals of the regression at the chosen width, in the new attribute
  `model`); [`formula()`](https://rdrr.io/r/stats/formula.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and a
  [`print()`](https://rdrr.io/r/base/print.html) method for
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
  whose result now keeps `formula` and the site-level `model`;
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) for
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
  and
  [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md);
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) for
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).

- The [`print()`](https://rdrr.io/r/base/print.html) methods of the
  result classes take `digits`.

- New [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  for
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md),
  [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
  (posterior against prior),
  [`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md)
  and
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md)
  (intervals, or trace plots of the chains with `type = "trace"`).

- The seed is kept with the result of
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md)
  (element `seed`) and of
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
  [`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
  [`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md)
  and
  [`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md)
  (attribute `seed`).

- The `scores` of
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  keep the site names of `comm`, and
  [`predict()`](https://rdrr.io/r/stats/predict.html) of
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
  adds a column `name` with the row names of `newdata`.

- [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
  [`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
  [`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md)
  and
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  gain `na_action`. With the default, `"omit"`, rows with missing values
  are left out as before, but now with a message giving their number
  ([`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md)
  and
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  dropped them silently); `"fail"` stops instead.

- [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md)
  and
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  check the collinearity of their predictors: they stop when one is a
  linear combination of the others (before,
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md)
  dropped it silently and
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  returned `NA` for its slope), and warn when a variance inflation
  factor exceeds 10. No variable is removed.

- Repeated site names now stop `check_sites()`,
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  (row names of `env`, `coords` and `areas`) and
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md)
  (site column of `env`), instead of matching sites by position or to
  the first repeated name.

- [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
  warns when the optimiser of the OU model does not converge and adds
  the column `converged` to `models`.

- The diagnostics of
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md)
  count transitions that hit the maximum tree depth (`max_treedepth`); a
  fit passes only when there are none.

- The neighbour metrics
  ([`mnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`sdnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`cvnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md))
  accept vectors with a class such as `units`, which they drop.

- [`animate_food_web()`](https://provetelab.org/comdynamite/reference/animate_food_web.md)
  adds the suffix `.gif` to `file` when it has none.

- The regional null models of
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  permute the whole region once per randomisation, instead of drawing
  each site from the region independently, as
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  now does. The null distribution of each site is unchanged, but null
  values differ from earlier versions for the same seed, and values that
  combine sites are now valid.

- The eigenvectors returned by
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  and
  [`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md)
  have a fixed sign: the largest absolute value of each is positive. The
  sign of an eigenvector is arbitrary and could differ between
  platforms; the signs of `eigenvector_correlation` follow. Tests and
  gradients (`pev`) are unchanged.

- Functions that implement a published method, such as
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md),
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md),
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  and
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
  show the original papers of the method the first time they run in a
  session, with a request to cite them together with the package.
  `options(comdynamite.cite = FALSE)` turns the message off. The new
  [`method_references()`](https://provetelab.org/comdynamite/reference/method_references.md)
  lists the papers of each function, and `citation("comdynamite")` now
  includes them.

- The null models of T_IC/IR and T_PC/PR in
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  permute the whole region once per randomisation, instead of drawing
  each site from the region independently, so that the sites of a
  randomisation partition the regional pool as the data do. The test of
  the mean over sites in
  [`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)
  now has its nominal size (it rejected in about 8% of null data sets
  for T_PC/PR and in none for T_IC/IR); the site-level null
  distributions are unchanged for T_IC/IR and nearly so for T_PC/PR, but
  null values, effect sizes and p-values differ from earlier versions
  for the same seed. New argument `null_pc`. Its default,
  `"populations"`, permutes population means among the populations of
  the region; `"individuals"` is the null of Taudiere & Violle (2016)
  and of earlier versions, which flags about 80% of sites as “high”
  without any filtering when population sizes differ among sites and
  species differ in their means. Results of `T_PC.PR` therefore change
  with the default.

- [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
  and
  [`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
  gain `min_node_size` (default 5, as before), the minimal size of the
  leaves of the forests. The help of
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
  [`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
  and
  [`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md)
  now says that the intervals of ABC-RF are approximate: conservative
  with several statistics (coverage 0.93 to 0.99 for a nominal 0.9 in
  the extended tests) and too narrow with one very informative statistic
  and small leaves (0.65, 0.89 with leaves of 200 simulations).

- The help of
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md)
  explains that the p-values are F tests of the regressions of the site
  values on the predictors, as in
  [`cati::traitflex.anova()`](https://adrientaudiere.github.io/cati/reference/traitflex.anova.html),
  and that under species turnover along the predictor the intraspecific
  component is underestimated and its test is conservative, because
  species means taken over all sites absorb part of a shared
  intraspecific response. Results are unchanged.

- Extended tests (rOpenSci standards G5.4-G5.12, BS7.0-BS7.2), run with
  `COMDYNAMITE_EXTENDED_TESTS=true` and once a month on GitHub Actions
  (`tests/README.md`): type I error and power of the permutation and
  null-model tests, recovery of known parameters over several seeds,
  convergence as data grow, agreement of the simulators with analytic
  results, and recovery and prior sensitivity of
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md).
  Quick checks that the results do not change with trivial noise in the
  inputs or beyond Monte Carlo error across seeds run with the regular
  tests.

- Input validation (rOpenSci standards G2 and G3): every scalar argument
  of the exported functions is checked for its length, type and range
  before any computation, with a message that names the argument, and
  the help pages say what each one accepts. Infinite values in `comm`,
  `env`, coordinates and traits stop the function, except where `Inf`
  has a meaning (no drift with `ne_ratio = Inf`, no speciation with
  `threshold = Inf`, a permanent change with
  `env_pulse(duration = Inf)`, sampling all individuals with `n = Inf`).
  Checks of zero variance tolerate rounding error. The
  [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) calls on user
  input became messages that say what is expected. Some inputs that were
  silently truncated or gave `NA` results now stop: fractional counts
  (such as `species = 2.5` in
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)),
  `NA` or `Inf` in `env` of
  [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
  character folds in
  [`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
  `seed` beyond about 1e9 in the ABC functions, and parameters that a
  model ignores but that were out of range. The package help page says
  that arguments that choose among options are case-sensitive.

- The Description and the help pages now say that
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md),
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  are translations of the MATLAB code in the supplements of Leibold et
  al. (2010) and Peres-Neto et al. (2012), with credit to their authors,
  and that
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
  was written from the equations and only compared with the MATLAB code
  of McPeek (2017).

- Preparation for rOpenSci review: a contributing guide with the
  development stage of the package, a code of conduct, the ORCID of the
  author, a project status badge, British English as the language of the
  documentation (with `inst/WORDLIST`), a `.lintr` configuration and a
  pkgcheck workflow (run by hand). The README explains how to install
  treeshapesim and cmdstanr.

- [`animate_food_web()`](https://provetelab.org/comdynamite/reference/animate_food_web.md)
  writes to a temporary file by default instead of `food_web.gif` in the
  working directory; the path is returned as before.

- Examples added for
  [`plot.population_shifts()`](https://provetelab.org/comdynamite/reference/plot.population_shifts.md)
  and `pp_check()` on
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md);
  the example of
  [`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md)
  draws only four random partitions and leaves the longer runs in
  `\donttest{}`.

- Internal: closures no longer assign with `<<-`,
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html) calls add to
  earlier ones, and an internal function used only by the tests moved to
  them. Results are unchanged.

## comdynamite 0.10.0

- `kernel_group` (routes of dispersal by group of species) now works
  with immigrant selection and with several evolving traits in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  where the immigrants of each species and the genes they carry follow
  the routes of its group, and in
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
  which takes one group per species at the start: new species stay in
  the group of their mother, the gene flow that keeps the populations of
  a species together follows the routes of its group, and the output has
  the group of each species (`kernel_group`). A run in which every group
  uses the same matrix gives the same results as a vector `kernel_step`
  (with several traits, up to rounding in the last digits).
- Fixed:
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  with point speciation and `immigrant_selection = TRUE` stopped with an
  error (“logical subscript too long”) when all the individuals of a
  species became point mutants in the same generation. Runs that
  finished before are unchanged.
- A niche with correlated axes
  ([`niche_correlation()`](https://provetelab.org/comdynamite/reference/niche_correlation.md))
  now works with `evolution` in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and with speciation in
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
  which refused it. The matrix of Gaussian selection is `Omega = D R D`,
  with `D` the diagonal matrix of the niche breadths and `R` the
  correlation of each species: the growth rate, with the cost of the
  phenotypic variance `sqrt(|Omega| / |Omega + P|)`, the response to
  selection `G (Omega + P)^-1 d`, the fitness of immigrants
  (`immigrant_selection`), the overlap of niches (`niche_overlap`) and
  the load of
  [`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md)
  all use it. A shared correlation applies to new species, and with one
  matrix per species they take that of their mother; `traits` in the
  output keeps the correlation of the extant species.
  `competition_matrix(type = "trait")` also counts the correlation. With
  the identity, the results are those of independent axes to the last
  bit. A niche breadth of 0 cannot be combined with a correlation.
- New article “Persistence in river networks under a flood pulse”: a
  generic study of how network topology, the direction of dispersal and
  seasonal floods change regional richness and persistence, and
  scenarios inspired by the Pantanal (normal pulse, a 12-year drought
  and regulation) with fish guilds that differ in emigration, routes and
  tolerance to drying.
- Dispersal and disturbance by species or group of species in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md):
  `kernel_group` gives each species a group, and a matrix `kernel_step`
  (steps x groups) chooses the matrix of `kernels` that each group
  follows at each step, so that, for example, long-distance migrants and
  sedentary fish of the same river follow different routes;
  `dispersal_site` multiplies emigration site by site and step by step
  (a dispersal pulse only where the floodplain floods); and
  `disturbance_sensitivity` (one value or one per species) is the
  probability that a population is lost when its site is disturbed, so
  that resistant species keep their populations. `dispersal_site` and
  `disturbance_sensitivity` also work in
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
  where new species take the sensitivity of their mother. The defaults
  give the same results as before, to the last bit.
- New
  [`env_seasonal()`](https://provetelab.org/comdynamite/reference/env_seasonal.md)
  (a seasonal cycle with an amplitude per site and anomalies per year)
  and
  [`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md)
  (washout probabilities by flood season, year and river size).
  `disturbance` in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  also takes a matrix of probabilities per site and step.
- New
  [`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md)
  and
  [`env_river()`](https://provetelab.org/comdynamite/reference/env_river.md):
  covariance of the environment along a river network as a sum of
  tail-up (flow-connected sites only, weighted by Shreve magnitude,
  length of channel upstream or any weight per site), tail-down and
  Euclidean exponential components plus a nugget (Ver Hoef et al. 2006;
  Ver Hoef & Peterson 2010; Peterson & Ver Hoef 2010), and Gaussian
  fields drawn from it, fixed or as autoregressive series for
  `env_time`.
  [`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)
  also returns `upstream_length`, a proxy for catchment area.
- New
  [`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)
  and
  [`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md):
  random dendritic river networks grown site by site (or any network
  given by downstream neighbours and reach lengths, such as an optimal
  channel network or a real river), and landscapes on them with
  dispersal along the channel whose scales differ downstream and
  upstream, plus an optional share of straight-line (overland)
  dispersal. The landscape keeps the network, the distances along the
  channel and the distance travelled downstream between every pair of
  sites.
- Simulations with `evolution` record the steps since the environment of
  each site last changed (`evolution$since_change`, at the end and at
  `record_steps`), and
  [`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md)
  gains a column `since_change`, to follow the lag load after a change
  in `env_time`.
- `env_surrogates(lonlat =, strata =)`: great-circle distances, and
  surrogates drawn within groups of sites such as cities.
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  take a list of surrogates as `null`.
- The articles on clades and on the whole tree apply the restricted
  nulls to `daphniids` and `anurans` (`data-raw/real_data_nulls.R`).
- `sample_individuals(visits =, session =)`: several visits to each site
  at the same step, with a column `visit`, and `session = "visit"` for a
  `session_sd` error shared by the individuals of one visit to one site
  rather than by all sites on a date. The help now says that with one
  step and `session = "step"` the error is a constant that vanishes from
  centred statistics. Repeated steps in `step` are an error pointing to
  `visits`.
- `anurans`: counts, one identification (*Scinax tymbamirim*, formerly
  *S. imbegue*), the depth of one pond, coordinates and urban cover as
  corrected by the collector in November
  2026. Urban cover is the share of MapBiomas Collection 9 urban pixels
        within 500 m, not 1 km as documented before.
- New
  [`env_surrogates()`](https://provetelab.org/comdynamite/reference/env_surrogates.md):
  predictors with the same spatial autocorrelation as the observed ones,
  by variogram matching (Viladomat et al. 2014) or a parametric Gaussian
  field matched by rank.
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  take them with `null = "variogram"` or `"parametric"` and `coords`,
  for landscapes with edges, where the torus rejects too often. Both
  still become liberal when the range of the field approaches the size
  of the landscape.
- [`env_surrogates()`](https://provetelab.org/comdynamite/reference/env_surrogates.md)
  also takes a landscape, and returns landscapes with the environment
  (or `env_axes`) replaced and the distances of the landscape kept, as
  null environments for the simulator; with `env_time` and
  `"parametric"`, it moves the whole series of each site.
- The “whole tree” article compares six null models for a continuous,
  autocorrelated predictor in landscapes with edges
  (`data-raw/validation_continuous_nulls.R`). Variogram matching
  (Viladomat et al. 2014) and a parametric surrogate came closest to the
  nominal rate, but every restricted null became liberal when the range
  of the field approached the size of the landscape.
- The validation article has an animation of
  [`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md)
  next to
  [`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md).
- `competition_matrix(type = "trait", traits = )`: competition by the
  overlap of the Gaussian niches of two species, the kernel of
  `niche_overlap` in
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
  with the product over axes for several axes.
- New
  [`niche_correlation()`](https://provetelab.org/comdynamite/reference/niche_correlation.md):
  a correlation between the responses of a species to the axes of its
  niche, so that
  [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md)
  and
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  use a multivariate Gaussian niche (one matrix for all species or one
  per species).
- `landscape_grid(n_axes, axis_cor)`: several correlated environmental
  axes, in `env_axes`, which
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  use with a list of traits when `env` is not given. The first axis is
  the `env` of a single-axis grid with the same seed.
- New
  [`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md):
  the environment of Thompson et al. (2020), a Gaussian field over space
  and time with exponential covariance, held constant during a burn-in
  and converted to ranks, as in the authors’ code. RandomFields is no
  longer on CRAN, so the field is drawn by random Fourier features, with
  the covariance exact in expectation. Unlike
  [`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
  sites close together change together.
- The help of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and the validation article now recommend `immigration = "multinomial"`
  for low dispersal, an environment that changes in time and studies of
  long-term persistence; the Poisson draw stays the default.
- [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  takes `dispersal_time`, `kernels` and `kernel_step`: dispersal that
  changes through time now works with speciation. New species inherit
  the emigration through time of their mother, and the gene flow that
  keeps populations together follows the kernel of each step.
- `immigrant_selection` in
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)
  now works with several traits (fitness under the multivariate Gaussian
  niche) and with speciation, where the gene flow counted for divergence
  includes only the immigrants that survive. `mutation_var` is still
  refused in both cases.
- A matrix of `kernels` whose rows already sum to 1 is no longer
  rescaled, so that the matrix of the landscape given as a kernel
  reproduces a run without `kernels` exactly.

## comdynamite 0.9.0

- New article “Dispersal and disturbance through time”: the shape of the
  dispersal kernel, recurrent pulses with `env_pulse(every =)`,
  `dispersal_time` and `kernels`, and a flood-pulse example that
  compares constant dispersal, dispersal pulses during high water,
  long-distance connectivity during floods and philopatry in the dry
  season.

- `dispersal_time`, `kernels` and `kernel_step` in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  dispersal that changes through time. `dispersal_time` sets the
  emigration probability at each step (multipliers of `dispersal`, or
  probabilities per step and species), and `kernels` with `kernel_step`
  switch among dispersal matrices, for example connected ponds during
  floods and isolated ponds in the dry season.

- `shape = "gaussian"` in
  [`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md)
  and `disp_shape` in
  [`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md)
  and
  [`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md):
  a Gaussian kernel, with a thinner tail than the exponential one (the
  default, unchanged). Kernel weights are now computed on the log scale,
  so that sites far apart relative to `scale` no longer underflow to
  zero.

- `every` in
  [`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md):
  recurrent pulses, one every `every` steps, so that both the duration
  and the frequency of disturbance can be set.

- [`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md):
  pulse (or permanent) changes of the environment of some sites, for the
  `env_time` of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  so that disturbance acts on growth rates through the niche rather than
  by removing individuals. New article “Replicating Hajian-Forooshani &
  Chase (2026)” (`data-raw/replication_hajian2026.R`,
  `data-raw/figures_hajian2026.R`): their pulse-disturbance model of
  scale-dependent diversity change, from the code archived with the
  article, compared setting by setting with their archived simulations
  (`data-raw/compare_hajian2026.R`; correlations of 0.998-0.999), with
  their Figs. 2-4 and a control without disturbance.

- `n_boot` in
  [`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md):
  a parametric bootstrap of the likelihood ratio against constant rates,
  with communities simulated under the constant rates estimated and
  refitted with and without the partition (`p_boot`, `null_boot`); it
  needs no room on the tree for random partitions, so it is available
  for regimes with a single shift, where the test against random
  partitions rarely has 19 distinct ones. `n_rand = 0` runs the
  bootstrap alone. `data-raw/s1_type1_shifts.R` gains `MODE=boot`.

- New article “Nesting map” (`data-raw/figure_nesting_map.R`): seven
  published models mapped onto comdynamite, with the settings that
  approximate each, the mechanisms it lacks and the outcome of the
  replications; comdynamite is presented as a modular framework for
  populations and quantitative genetics in ecological time, not as a
  model that contains the others.

- `recruitment = "lottery"`, `microsites` and `mortality` in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  adults survive with probability `1 - mortality`, offspring are
  produced in proportion to fitness, and the microsites left vacant are
  filled by a multivariate hypergeometric draw from the offspring
  present after dispersal; `mortality < 1` gives overlapping
  generations, with drift and selection scaled to the share of recruits
  (one trait, one step per generation).

- `variance_selection` in
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md):
  with `mutation_var`, the genetic variance also loses the Bulmer (1971)
  reduction by stabilising selection each generation, kept as under
  clonal inheritance, so it approaches a mutation-selection-drift
  balance.

- Sections “Three missing mechanisms” in the articles replicating Urban
  & De Meester (2009) and Leibold et al. (2019)
  (`data-raw/mechanisms_urban2009.R`,
  `data-raw/mechanisms_leibold2019_recruitment.R`): erosion of the
  variance by selection restores the priority effect of standing
  variation, and overlapping generations restore the advantage of the
  better disperser under evolution, in the Urban model; none of the
  three mechanisms produces the local dominance of Leibold et al.

- `pool_from` in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  a sites x species matrix of the first step at which each species can
  occur in each site, for species introduced at different times; until
  then the site is not a destination of its emigrants, which stay home.

- New article “Replicating Urban & De Meester (2009)”
  (`data-raw/replication_urban2009.R`): a species that evolves by
  mutation wins a new site whatever the order of arrival, and evolution
  in both species strengthens the priority effect, as in their model,
  though less; standing variation reproduces their result only with a
  constant genetic variance or with `variance_selection`, since
  `mutation_var` alone carries a variance load; over the dispersal grid,
  evolution makes the poorer disperser win, as in their gene-flow
  maladaptation region, and the better disperser keeps its advantage
  only with overlapping generations.

- `disturbance` in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md)
  and
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md):
  the probability that a site loses all its populations at each step, as
  the patch extinctions of Leibold et al. (2019).

- New article “Replicating Leibold et al. (2019)”
  (`data-raw/replication_leibold2019.R`): with fast enough adaptation,
  both species adapt to both patch types, fitness differences vanish and
  patch type stops predicting composition while it still predicts
  traits, as in their model; local dominance is weaker and regional
  drift removes one species faster. A factorial test
  (`data-raw/mechanism_leibold2019.R`) traced the faster drift mainly to
  drift and founder sampling of the population means, not to the merging
  of immigrants; with a genetic variance that depends on population
  size, no species was lost, but local dominance disappeared at
  intermediate and high dispersal.

- `mutation_var` in
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md):
  each population carries its own genetic variance, with a loss of 1 /
  (2 Ne) and a gain of the mutational variance per generation (Lynch &
  Hill 1986), a founder bottleneck and averaging by gene flow (one
  trait).

- `immigrant_selection` in
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md):
  immigrants from each source arrive in proportion to the fitness of
  their phenotype in the site they reach, before they mix with residents
  (one trait, spatially explicit landscapes).

- Speciation with several traits. With a list of `traits` and an `env`
  matrix (one column per niche axis),
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  evolves one trait per axis in every lineage, with (G) and (E) as
  covariance matrices in
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
  as in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
  Ecological speciation compares populations by their Mahalanobis
  distance under (P = G + E) (`threshold` in units of phenotypic
  standard deviations), `jump_var` of
  [`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)
  takes a covariance matrix for the multivariate normal jump of the
  optimum at speciation, and niche crowding multiplies the overlaps on
  every axis. The outputs carry one table of `traits` per axis,
  `genetic` and the evolutionary state as arrays with a third dimension
  for the traits, and `genetic_<trait>` columns in `series`;
  [`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  and
  [`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md)
  take them, and the last two now add the variance among the lineages of
  a species in a site. `env_time` takes several axes in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  as a list of matrices or an array sites x steps x axes.

- New
  [`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md)
  tests whether divergence among groups is proportional to the variation
  within them, as expected under drift (Lande 1979), by regressing the
  log variance among group means on the log eigenvalues of (W)
  (Ackermann & Cheverud 2004), with (W) the pooled covariance within
  groups or a known (G), after subtracting the sampling variance of the
  means. In a validation in `data-raw` with clades simulated by
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
  the test with a known (G) kept its size under drift (3.6% at 5%) and
  gave slopes below 1 in every run under divergent selection along the
  minor axis of (G), but rejected drift in only 36% of those runs with
  eight traits; with (W) estimated from phenotypes and low heritability
  on the minor axes, it rejected true drift in 57% of runs (Prôa et
  al. 2013).

- New
  [`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md)
  gives, for every population of a simulation run with `evolution`, the
  deviation of its mean phenotype (genetic mean plus plastic response)
  from the environment of its site and the lag load, the share of growth
  lost to it; with several traits, through the Mahalanobis distance
  under (+ P). It takes recorded steps and evolving metacommunities. A
  new article, “Quantitative genetics in the metacommunity”, animates
  migration load along a gradient, evolutionary rescue under warming,
  genetic constraints set by the orientation of (G), and evolving
  priority effects that obscure the relationship between species and the
  environment (Urban & De Meester 2009; Leibold et al. 2019), and
  relates the model to the mechanistic models of Cabral and colleagues.

- Bug fix:
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  with `evolution` and a niche of infinite breadth (`sigma = Inf`) gave
  growth rates of `NaN`.

- Evolution of several traits at once. With a niche of several axes,
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  evolves one trait per axis (its optimum), with `genetic_var` and
  `env_var` of
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)
  as covariance matrices (G) and (E), one for all species or one per
  species: growth falls with the Mahalanobis distance of the phenotype
  to the environment under (+ P), the response to selection is (G (+
  P)^{-1} d) (Lande 1979), drift is multivariate normal with covariance
  (G / N_e), and an axis with `sigma = Inf` carries a trait under no
  selection that changes only by its genetic covariance with the
  selected ones (correlated response) and by drift. The new
  [`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md)
  builds (G) and (E) for each species in three families (proportional to
  a common (P), with no common structure, or with shared principal
  components), and `sample_individuals(source = "evolution")` draws
  individuals from (N(z, G + E)), one column per trait. Behaviour
  change: `evolution` with a niche of several axes used to stop with an
  error; with a single variance it now evolves one trait per axis with
  (G = vI).

- [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  simulates collections on several dates: `step` takes several recorded
  steps, each sampled on its own, with a column `step` in the result,
  `env` takes a matrix with one column per step, and `session_sd` adds
  an error shared by all individuals measured on the same date. With
  `source = "evolution"`, any recorded step can be sampled:
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  keeps the environment at `record_steps`
  (`evolution$environment_history`), and
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  keeps the genetic means, the variance among lineages and the
  environment of the extant species at `record_steps`
  (`evolution$history`, `within_history`, `environment_history`).

- Phase 3b of the block on speciation.
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  logs every change of presence of every lineage in every site
  (`ranges`), and the new
  [`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md)
  maps them on the full or reconstructed tree as a stochastic map (class
  `simmap`) whose states are the combinations of areas occupied (`"A"`,
  `"B"`, `"A+B"`).
  [`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md)
  takes evolving metacommunities, with this history as the true areas of
  the ancestors, and
  [`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
  [`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md)
  and
  [`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md)
  accept ranges of several regions.
  `history_rates(model = "sympatry", overlap = TRUE)` counts the
  lineages whose ranges share an area.

- Phase 3a of the block on speciation.
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  returns a per-generation `series` of every lineage (occupancy,
  abundance, genetic mean, and sites colonised and lost out of those
  exposed) and, with `evolution`, the evolutionary state of the extant
  species at the end (`evolution`). New
  [`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md)
  estimates the rates of colonisation and extirpation of a two-state
  Markov model of presence, and from them the equilibrium occupancy
  (`epsilon`, on the logit scale) and the resilience (`tau`, the log
  turnover rate) of Januario et al. (2025), for evolving metacommunities
  (per lineage, over a `window` of generations), for the records of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  and for arrays of sites x species x times.
  `sample_individuals(source = "evolution")` draws individuals from the
  genetic mean plus the plastic part of each population with the genetic
  and environmental variances, for
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
  and accepts evolving metacommunities. New
  [`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md)
  builds a landscape from a given dispersal matrix, such as a kernel of
  resistance distances.

- [`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)
  gains the modes `"point"` (each individual starts a new species with
  probability `rate` per generation, the neutral limit) and `"fission"`
  (species split at `rate` per generation into the populations nearest a
  random site, `partition = "spatial"`, or a random partition), a
  Gaussian shift of the genetic mean of new species (`jump_var`), and
  protracted speciation (`completion`, fixed or exponential, as in the
  protracted birth-death model). Incipient lineages are separate in the
  dynamics; at the end, the extant lineages are grouped into good
  species, which
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  returns as `abundance`, `tree` and `traits`, with the lineages in
  `lineages` and the times of initiation and completion in `species`
  (phase 2b of the block on speciation). The full tree is now built
  without recursion and the split times are kept only among alive
  lineages, so runs with thousands of lineages are fast. Validations in
  `data-raw/evolving_neutral.R` (point mutation against the Ewens
  sampling formula) and `data-raw/evolving_pbd.R` (protracted fission
  against the gamma statistic and fits of the PBD).

- [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  runs the metacommunity with speciation (phase 2a of the block on
  speciation): species bud when their populations diverge, and die out
  with their last population. The largest group of populations keeps the
  identity and age of the species, and the others become new species,
  branching off in order of divergence (internal branches of length zero
  when several bud at once). `speciation_control("geographic")` counts
  divergence between clusters of populations linked by a migration rate
  of at least `migration_min`, frozen within a cluster and reset only
  after `reset_generations` generations in it, or lowers it in
  proportion to the migrants exchanged (`divergence = "gradual"`).
  `"ecological"` splits groups whose genetic means differ by the
  threshold for `min_generations` generations in a row with a migration
  rate below `migration_max`, so that drift in small populations does
  not split species in a uniform environment. The tree grows out of the
  dynamics (full and reconstructed), competition decays with the
  patristic distance down to `inter_min` and, with
  `niche_overlap = TRUE`, with the difference of phenotypes in each
  site, and
  [`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md)
  turns the output into an `mc_sim`. Without a split the simulation is
  that of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  with the same seed.
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md)
  now runs its preparation, steps and end as internal functions on one
  state, with the same results.

- [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)
  and the argument `evolution` of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md)
  let the niche optimum evolve within each population (phase 1 of the
  block on speciation, with no speciation yet): selection by the local
  environment as the response to Gaussian selection of Lande (1976),
  drift with an effective size proportional to abundance, gene flow by
  the mean of residents and immigrants weighted by their numbers,
  sampling of founders in new populations, a linear reaction norm for
  plasticity whose part is reported apart, and steps per generation.
  Without genetic variance the abundances are the same as without
  evolution for the same seed. The tests check drift against
  Lande (1976) and the divergence of two populations against
  migration-selection balance.

## comdynamite 0.8.0

- [`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
  goes through inference by ABC with random forests with the anurans of
  40 urban ponds in two Brazilian cities (`anurans`, new data set, with
  regional pools from GBIF occurrences and a dated tree pruned from
  Portik et al. 2023): simulator, reference table, recovery, model
  choice, posterior, diagnostics, posterior predictive checks and
  cross-validation over spatial blocks (`data-raw/abc_anurans.R`).

- [`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md)
  builds a landscape from the coordinates of real sites, with
  great-circle distances for longitude and latitude. `daphniids` gains
  `coords`, the position of each lake, and its `tree` now has branch
  lengths from TimeTree 5 (with the topology of TimeTree; four species
  placed at arbitrary ages), while `tree_leibold` keeps the topology of
  Figure 3a of Leibold et al. (2010) for analyses that reproduce the
  article.

- [`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md)
  turns the true abundances of a simulation into survey records: each
  individual is detected with a given probability in each of several
  visits (binomial counts, as in N-mixture models), keeping the largest
  count or presence, for the final community and for each step of the
  series. Inside `simulate` of
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
  it puts observation error in the reference table. In
  `data-raw/abc_detection.R` detection was recovered and the other
  parameters mostly kept their recovery, while ignoring detection biased
  the niche breadth. `data-raw/abc_capacity.R` adds the carrying
  capacity as a parameter: detection times capacity was recovered well,
  each of them only in part.

- [`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md)
  gives, for each observation and parameter of an
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
  the effective sample size and the Pareto k-hat of the forest weights,
  with posterior in Suggests, and the statistics of the observation
  outside the central 99% of the reference table. In
  `data-raw/abc_diagnose.R` k-hat did not predict wrong estimates, while
  a low effective sample size and statistics outside the table flagged
  data from another model.

- Cross-validation of ABC across sites.
  [`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)
  splits the sites into spatial blocks (k-means on the coordinates) or
  random folds,
  [`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md)
  turns a summary function into one that gives the statistics of the
  sites kept and left out in every fold, for use in
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
  and
  [`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md)
  fits
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
  on the sites kept and predicts the statistics of the sites left out,
  with coverage, posterior predictive p-values and standardised errors.
  `data-raw/abc_cv.R` compares spatial blocks with random folds: both
  were calibrated under the right model, and blocks did not flag a wrong
  model more often on a 5 x 5 torus.

- [`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md)
  chooses among models with a classification random forest (Pudlo et
  al. 2016) and reports the out-of-bag error and confusion matrix;
  [`predict()`](https://rdrr.io/r/stats/predict.html) gives the chosen
  model, the votes and its posterior probability.
  [`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md)
  runs posterior predictive checks after
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
  from the reference table or with new simulations, and gives a
  two-sided posterior predictive p-value per statistic.
  `data-raw/abc_model_choice.R` tests both on phylogenetic competition
  against no competition: the models were told apart with an out-of-bag
  error of 0.04, while the check on temporal statistics had little
  power, since the competition model contains the neutral one.

- Inference by approximate Bayesian computation with random forests
  (Raynal et al. 2019), with abcrf in Suggests.
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
  draws parameters from uniform priors and simulates a reference table,
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
  fits one regression forest per parameter,
  [`predict()`](https://rdrr.io/r/stats/predict.html) gives the
  posterior median and a credible interval for observed statistics,
  [`abc_importance()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
  ranks the statistics, and
  [`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
  measures R2, bias and coverage in pseudo-observations left out of the
  table. `data-raw/abc_recovery.R` repeats the design of the ABC pilot:
  the forests recovered the parameters about as well as rejection ABC at
  its best tolerance, without having to choose one, with calibrated or
  conservative intervals.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  has a `clades` argument. It takes either the number of the deepest
  clades of the tree or a vector of clade names, such as families, named
  by species, and draws the species as stars coloured by clade, each
  clade with its convex hull and name, as in Fig. 1 of Debastiani &
  Duarte (2014). With clades, the predictor fills the sites. Sites are
  now drawn below the species.
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  keeps the pruned tree as `tree`.

- The article “Testing the whole tree” cites the origin of the principal
  coordinates of phylogenetic structure and explains why the curve of
  phylogenetic signal of the package PCPS (`pcps.curve()`) is not
  implemented.

- [`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md)
  no longer draws partitions that split the tree as the original one,
  such as the complement of a clade at the root: `max_overlap` now also
  bounds the adjusted Rand index between the random and the original
  partition. In simulations with one shift in the rate of speciation and
  no effect on local rates, such draws made
  [`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md)
  reject in 12% to 40% of replicates.
  [`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md)
  draws the random partitions before fitting them, so its p-value no
  longer depends on `cores`, and counts each distinct random partition
  once in the p-value (repeated draws of a few placements made it too
  small); it returns their number (`n_distinct`) and warns when there
  are fewer than 19, as the p-value then cannot fall below 0.05. Trees
  with a strong shift in diversification often leave that little room.

- The article “Rates and assembly” rejects ES-sim at p ≤ 0.05, the rule
  of the preprint, instead of p \< 0.05. With 999 null draws five cells
  change by one replicate each and no conclusion changes.

## comdynamite 0.7.0

- New article “Does metacommunity assembly leave a trace in branch
  rates?”, a recovery study with planted signals. Branch rates are
  weakly identified. Partner turnover planted in the rates is barely
  recovered among species, but local adaptation and character
  displacement are recovered among populations. The T-statistics
  compensate when both act, and displacement is identifiable only when
  the competition kernel is not much narrower than the spacing of
  species. With fully sampled clades, a PGLS of DR on metapopulation
  persistence, as in Januario et al. (2025), has no power under a causal
  link, whereas ES-sim on turnover has. Simulations are in
  `data-raw/bridge_*.R` and `data-raw/januario_*.R`.

- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md)
  record the colonisation history of every population present at the
  last step: `arrival` (the step of the last founding, reset by
  recolonisation after local extinction; 1 for populations present from
  the start), `propagule` (the number of founders), `immigrants`
  (immigrants and seeds since arrival) and `migrant_load` (the summed
  share of newcomers in the population since arrival), all sites x
  species matrices, plus `n_steps`. No random draws are added, so
  abundances are unchanged for the same seed.
  [`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md)
  returns them as a table.

- [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  has two arrival-dependent population deviations, both off by default
  (which keeps earlier results for the same seed):

  - `founder_sd`: founder effects, a deviation drawn at colonisation
    with standard deviation `founder_sd / sqrt(founders)` that persists
    and is eroded by gene flow by `exp(-gene_flow * migrant_load)`;
  - `monopolization`: evolution-mediated priority effects (Urban & De
    Meester 2009), in which a population closes the gap between the
    optimum of its species and the environment of its site at a rate per
    step since its arrival, toward a balance between selection and
    migration set by `gene_flow`. Early colonists end up better matched
    to their site than late arrivals. With either, the result gains an
    `arrival` column. Neither feeds back on the dynamics.

- [`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md)
  has a `crossed` argument for a second hierarchy crossed with
  `factors`, such as genus \> species crossed with region \> site, with
  a level for the interaction of the finest levels (the population,
  species by site; `interaction = TRUE`). When species occur at several
  sites, a nested `factors = c("region", "site", "species")` makes
  `species` the population and mixes differences among species with
  differences among populations; in a simulation with 82% of the
  variance among species and 9% among populations, the nested model gave
  88% to `species`, the crossed one 78% and 11%. The crossed-factors
  warning now points to `crossed`. Crossed models, and nested ones that
  [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) fails to fit,
  are fitted with
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) (new in
  Suggests); singular fits warn, and the attributes `singular` and
  `engine` record them.

- [`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md)
  has two nulls for conserved niches: `null = "phylo"`, the phylogenetic
  permutations of Lapointe & Garland (2001), in which species means move
  preferably among relatives, and `null = "bm"`, means simulated by
  Brownian motion. In simulations where the trait and the niches both
  evolve along the tree, the phylogenetic permutations bring the error
  from 14–23% (exchange among all species) to 2–11%; `"bm"` is exact
  there but fails for traits without phylogenetic structure.

- [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md)
  fits in Stan (through cmdstanr, in Suggests) a bivariate phylogenetic
  mixed model in which the position of species along the gradient (beta)
  and their mean trait are latent effects estimated from every site. It
  returns the phylogenetic share lambda of beta, alpha and the mean with
  credible intervals, the phylogenetic correlations, species values and
  convergence diagnostics. Lambda has its own prior (`lambda_prior`,
  uniform by default) instead of the U-shaped prior implied by
  independent priors on two standard deviations; `ou_alpha` gives an OU
  structure with a given strength; `pp_check()` works on the result. In
  simulations with about 40 species it recovered lambda near 0 and 1,
  and estimated intermediate lambda low with wide intervals (see
  Limits).

- [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
  runs the trait-gradient analysis of Ackerly & Cornwell (2007): each
  species’ trait splits into its position along the gradient of
  community trait means (beta) and its difference from co-occurring
  species (alpha), with niche breadth and the slope within species.
  Changes from the original:

  - Site means leave out the focal species by default, which removes the
    circularity that inflates the correlation of beta with the species
    mean and with alpha; `leave_one_out = FALSE` reproduces the original
    script.
  - The weights of site means and of species means are chosen
    separately;
    [`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md)
    compares the four combinations.
  - [`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md)
    exchanges the trait means of species among species by default,
    keeping the communities as observed; a curveball null (richness of
    sites and frequency of species fixed) and the original permutation
    are options. In simulations, the variance tests under the curveball
    and the original permutation reject far too often when species turn
    over along the gradient, so they are for the correlation of beta and
    alpha only.
  - [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
    fits a model without phylogeny, Brownian motion and an OU process to
    beta, alpha and the species mean by REML, with the bootstrap
    standard errors of each value added to the covariance, and compares
    them by AIC.
  - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
    gradient, beta against alpha, or the null distributions.

- [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md):

  - The grid of kernel widths no longer goes below the pooled standard
    deviation within populations, and a warning is given when the best
    width is at an end of the grid or the fit barely changes along it.
  - `env` takes several environmental variables, and
    `species_slopes = TRUE` gives each species its own slope.
  - `community` (a site x species matrix of abundances) lets species
    that were present but not measured push the populations measured at
    their sites, with the mean of their populations elsewhere moved
    along the gradient by the fitted slopes.
  - `nboot` runs a parametric bootstrap of the design. It gives a
    p-value for displacement that accounts for the choice of the width
    and is centred on the mean of the simulations without displacement,
    a first-order bias correction of `c`, and intervals for `c` and the
    slopes. The simulations without displacement come from the fit
    without displacement, the species constants from the populations
    that enter the fit, and the variance of populations from a moment
    equation weighted by precision. Taking the slopes of the full fit
    inflated that variance and the bias of `c` under the null, so the
    test rejected too often in large designs (up to 29% at the 5% level
    in the experiments of `data-raw/bridge_experiments.R`). The test is
    now somewhat conservative (about 2–3% at the 5% level in simulations
    without displacement), because the centre of the null simulations
    follows the estimate.
  - The documentation now gives the sign of the bias that deviations of
    whole populations cause in `c`: upward when most competitors are
    closer than `w` in trait, downward when most are farther.
  - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
    profile of R², `c` and `h` along the width.

- All plots now use ggplot2 (in Suggests, version 3.5.0 or later) and
  return a ggplot object, which can be changed further with ggplot2:
  [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)
  (both layouts),
  [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md),
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md),
  [`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md),
  [`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md),
  [`plot_assembly()`](https://provetelab.org/comdynamite/reference/plot_assembly.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
  and the frames of
  [`animate_food_web()`](https://provetelab.org/comdynamite/reference/animate_food_web.md).
  Tables these functions used to return invisibly are kept as attributes
  (`marks`, `values`, `fractions`, `ordered`). Colours are `colours` and
  text sizes `text_size` everywhere; `col`, `cex`, `lwd` and logical
  `labels` (now `tip_labels` or `species_labels`) still work, with a
  warning, until 1.0.0. `legend` is `TRUE` or `FALSE`; a position
  keyword is still accepted. Arguments of base graphics passed in `...`
  are no longer used.

- Renamed for consistency; the old names still work, with a warning,
  until 1.0.0:
  [`traits_ou()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  is now
  [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md)
  (like
  [`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md),
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  and
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md));
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  is now
  [`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md);
  the `nperm` argument of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  and
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  is now `n_perm`, as in the rest of the package;
  `trait_rates(root_value)` is now `root`, as in
  [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md);
  `trait_mk(root)` is now `root_state`, since it is a state;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md)
  takes `colours`, as the other plots.
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md)
  returns the rates as `edge_rates`, as
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  does (`edge_rate` is kept until 1.0.0).

- `tstats_summary(level = )` is now the coverage of the null
  distribution (default 0.95), as in
  [`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md);
  values below 0.5 are read as the size of the test, as before, with a
  warning.

- [`plot_rates()`](https://provetelab.org/comdynamite/reference/plot_rates.md)
  draws any rate per branch on the tree, from
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
  [`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
  [`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
  [`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md)
  or an estimator, with ggplot2.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  draws the sites on two principal coordinates, coloured by a predictor,
  with the species as arrows;
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  now also returns `share`, `species` and `env` for it. The vignette
  “Rates of trait evolution that vary among branches” and the article on
  the whole tree use them.

- Clearer errors: functions that take a tree say so when given something
  else, instead of repeating a failed condition, and
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md),
  [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md),
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  and
  [`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md)
  name the argument at fault.
  [`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
  and the functions that call it stop when no species of `comm` is a tip
  of the tree.

- Help pages link to the other functions of their group (`@family`),
  following the groups of the reference index. The example of
  [`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md)
  uses a tree from treeshapesim.

- Plots for several traits, in base graphics:
  [`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md)
  draws traits from
  [`traits_ou()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  with their tree, as a phylomorphospace (two traits) or a traitgram
  (one trait), with branches coloured by regime and the colour changing
  where a mapped history changes;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  [`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md)
  draws each regime’s rate matrix as an ellipse;
  [`plot_assembly()`](https://provetelab.org/comdynamite/reference/plot_assembly.md)
  draws a community from
  [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)
  in trait space, with the boundary of the filter or the contours of
  limiting similarity in the metric of the assembly. New article,
  “Correlated traits and assembly in trait space”.

- [`traits_ou()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md):
  `theta` can be a matrix with one column per trait, so that a regime
  shifts each trait’s optimum by its own amount; `root` can have one
  value per trait and is checked; a mapped history is matched to the
  tree by its edges, so it may come in another order of edges. A named
  `theta` or `sigma2` of length one is no longer recycled to every
  regime: it must name every regime, as before rate matrices were added.

- [`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md)
  returns an object of class `rate_matrix`, with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
  [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)
  returns an object of class `mc_assembly`, which keeps the rule,
  distance, optimum, traits and the matrix of the metric, with a
  [`print()`](https://rdrr.io/r/base/print.html) method; it stops on
  missing or non-numeric traits instead of keeping such species or
  failing later.

- The limiting-similarity rule of
  [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)
  keeps the nearest neighbour of each species up to date instead of
  scanning the whole distance matrix at each removal: 23 times faster
  for a pool of 2000 species reduced to 100 (from 10.6 to 0.46 s), with
  the same species removed and the same random draws.

- Simulations with a niche of several axes keep their tree.
  [`as_mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  keeps one column of optima per axis,
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
  uses the first axis, and
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  stops with a clear message; both use the `env` that set growth when
  `simulate_metacommunity(env = )` was given one. Axes of `env` without
  names take the names of the list of traits, or `axis1`, `axis2`, …

- Plots of null-model tests:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  and
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  draws, site by site, the null interval (hatched by default) with the
  observed value coloured by the test, or the standardised effect sizes;
  [`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md)
  returns the table behind it.

- [`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md)
  gives trait probability densities of populations and communities
  (Carmona et al. 2016, 2019) with one bandwidth for all populations,
  and
  [`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md)
  draws them with adjustable transparency.
  [`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md)
  warns when the spread within populations rises with their mean (one
  bandwidth then fits badly; the logarithm of the trait usually helps),
  and `plot_population_traits(highlight = )` draws a few species in
  colour and the rest in grey, for sites with many species.

- Classic metrics of trait spacing and spread, for use with
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md):
  [`mnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`sdnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`cvnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`min_nnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`max_nnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`mnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`sdnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`cvnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`mst_length()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`mst_ratio()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`trait_range()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
  [`trait_kurtosis()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md)
  and
  [`nn_distances()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md).
  Zero distances count (cati discards them), and one-trait metrics use
  sorted values instead of distance matrices.

- [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  accepts a named list of indices, computed on the same randomised
  communities, and returns them in one table with an `index` column.

- Sampling effort:
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  gains `cover` (measure only the most abundant species of each site, as
  in Pakeman & Quested 2007) and `min_n`;
  [`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md)
  applies the same schemes to real data, and
  [`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md)
  follows any statistic over the number of individuals per population or
  over `cover`, with a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
  marks the 80% rule of Pakeman & Quested (2007) with a dashed line. Its
  default statistic computes the T-statistics without permutations.

- [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md):
  the null model of `T_IC.IR` now keeps the observed regional variance
  as the denominator. It was computed from the values drawn, which
  repeat individuals across sites; this biased the null upward and made
  the tests anti-conservative (about 6% of sites rejected per tail at a
  nominal 4% in simulations with random traits).

- New functions for analyses of intraspecific trait variation, written
  from the published definitions to replace functions of cati 0.99.6
  that give wrong results (an audit found bugs reported as
  adrientaudiere/cati#13 to
  [\#21](https://github.com/diogoprov/comdynamite/issues/21)):
  [`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)
  summarises the tests of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  (shares of sites below and above the null, test of the mean over
  sites, the sites that deviate), with sites that cannot be tested left
  out instead of given small p-values;
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md)
  splits community-weighted means into species turnover, intraspecific
  variability and covariation (Lepš et al. 2011), with environmental
  variables matched to sites by name, and a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  (ggplot2) with explained and residual sums of squares;
  [`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md)
  gives REML variance components across nested levels (nlme in Suggests)
  and stops when the factors are not ordered from the broadest to the
  finest level, with a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
  draws stacked bars of the shares in Okabe-Ito colours (ggplot2 in
  Suggests);
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  compares any trait index of communities with four null models, with
  all traits of an individual kept together in multivariate indices. The
  standardised effect size is now `NA` when the null values do not vary,
  as for an index that a randomisation cannot change.

- [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  estimates local adaptation (`h`) and character displacement (`c`)
  together from the population means of individuals, by a regression
  within species on the environment of the site and on the push
  predicted from co-occurring competitors, with the individuals of each
  population split in halves against the bias of sampling error and the
  kernel width profiled. `tstats(null_ip = "species")` adds a null model
  for T_IP.IC that exchanges population means among the populations of
  each species, so that its effect size is near zero when the site does
  not matter, positive under convergence and negative under divergence;
  the null model of Taudiere & Violle (2016) is far below zero for any
  set of distinct species. In `data-raw/pilot_tstats_identifiability.R`
  (ABC, `h ~ U(0, 1)`, `c ~ U(0, 10)`, ten individuals per population,
  metacommunities of the pseudo-observations left out), the T-statistics
  alone recover `h` (R2 0.87 to 0.89 over tolerances) but hardly `c`
  (0.56 to 0.58), while the estimates of
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  recover both (0.99, and 0.95 to 0.96; recoverable at tolerance 0.05).
  `group_sum()`, used inside
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
  now gives zero to empty groups, and
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  no longer fails on a trait without values.

- `sample_individuals(displacement = , displacement_width = )` adds
  character displacement: the populations of each site sit at the
  equilibrium between selection toward their mean without competitors
  and a push away from close competitors through a Gaussian competition
  kernel, weighted by their abundance.
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  now has the three null models of Taudiere & Violle (2016) in base R
  (`nperm`, `alternative`), with standardised effect sizes and p-values,
  so cati is needed only for comparison. The null model for T_IP.IC in
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
  0.99.6 puts together permutations of different sites, so its effect
  sizes for T_IP.IC are wrong (adrientaudiere/cati#12); its
  regular-expression match of site names was also reported upstream
  (adrientaudiere/cati#11). In `data-raw/demo_tstats.R`, character
  displacement lowers T_IP.IC and raises T_PC.PR and offsets local
  adaptation: with `local_adaptation = 0.9`, `displacement = 10` takes
  T_IP.IC from 0.61 to 0.011 and T_PC.PR from 0.010 to 0.59.

- New vignette “Rates of trait evolution that vary among branches”
  (`branch-rates`): rates as a model
  ([`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md)),
  from a history along the tree
  ([`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md))
  and from ecological dynamics
  ([`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
  [`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md)),
  and how well two estimators recover them, computed in advance by
  `data-raw/vignette_rates.R`. RRphylo follows the change that happened
  on each branch (mean Spearman rho 0.35 to 0.37) but hardly the true
  rate (0.06 to 0.08, against a ceiling of 0.24 to 0.35 set by one
  realisation per branch); Brownian motion with regimes in mvMORPH
  recovers the ratio between the rates of two habitats (medians 1.02,
  2.12 and 5.46 for true ratios 1, 2 and 5); and the RRphylo rates of
  tip branches follow the turnover of partners weakly (median rho 0.10).

- `fit_damocles(shared = TRUE)` and `damocles_test(shared = TRUE)`
  estimate one pair of DAMOCLES rates for several communities on the
  same tree from the sum of their log-likelihoods. A single community
  rarely holds enough information on colonisation and extinction;
  pooling communities narrows the estimates
  (`data-raw/validation_damocles.R`), provided they can be treated as
  independent realisations of the process.

- [`damocles_samples()`](https://provetelab.org/comdynamite/reference/damocles_samples.md)
  fits DAMOCLES rates over a list of configurations of regimes or edge
  covariates, such as posterior samples of rate shifts or of
  branch-specific rates of speciation, and summarises how often they
  beat constant rates, the parameters that mean the same in every fit,
  and the rate of each edge.
  [`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md)
  tests a partition against random partitions of the same complexity
  placed elsewhere on the tree, drawn by
  [`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md)
  (`max_overlap` excludes placements that nearly coincide with the one
  tested), which asks whether the partition is better placed than
  arbitrary ones rather than only whether rates vary; with the effect in
  another clade, the chi-squared test rejected constant rates in 29 of
  40 replicates and this test in none.
  [`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md)
  also takes node numbers.
  [`dr_rates()`](https://provetelab.org/comdynamite/reference/dr_rates.md)
  computes the tip diversification rate DR of Jetz et al. (2012), per
  tip or per edge, to use as an edge covariate.
  `data-raw/damocles_bamm.R` and `data-raw/damocles_clads.R` show the
  path from BAMM (BAMMtools) and ClaDS (RPANDA) output to these
  functions; with continuous rates, single posterior samples attenuate
  the slope, so the help recommends fitting the posterior summary of
  each rate.

- DAMOCLES rates can change across the tree.
  [`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)
  and
  [`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
  take one rate per edge;
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  and
  [`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md)
  estimate rates that differ among `regimes` of the edges or change with
  an `edge_covariate`, for colonisation, extinction or both (`vary`).
  Regimes can be clades defined a priori
  ([`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md))
  or rate shifts from methods of diversification (BAMM, ClaDS and
  similar), given as plain vectors so that no such package is needed.
  [`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md)
  compares nested fits (constant against varying rates, shared against
  separate, allopatric against estimated `sympatry`) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) gives
  profile-likelihood intervals. Type I error and power of the comparison
  with constant rates are in `data-raw/validation_damocles.R`.

- [`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md),
  [`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  and
  [`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md)
  take `sympatry`, the probability that a locally present ancestor
  leaves both daughters present at speciation (in situ speciation, or
  recolonisation before the daughters diverge), relaxing the strictly
  allopatric rule of Pigot & Etienne (2015), which stays the default
  (`sympatry = 0`). The likelihood remains exact;
  `fit_damocles(sympatry = NULL)` estimates it, which needs several
  communities with shared rates (`data-raw/validation_damocles.R`).

- [`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md)
  turns the rates of trait change of the resource, consumer and predator
  from
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
  into rates per branch for `trait_rates(edge_rates = )`: each tip
  belongs to a trophic level, and over the last `window` units of time
  of its branch the rate is `sigma2 * (1 + k * speed)`, where speed is
  the mean rate of change of the mean trait of that level while present.
  As in
  [`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
  only tip branches are informed; the speeds are directional changes
  under selection, so they give relative rates among levels, not
  Brownian variances.

- Community assembly with several traits.
  `assemble_community(distance = )` measures trait space by Euclidean
  distance (`"euclidean"`, the default and the earlier behaviour),
  Euclidean distance on traits standardised by their standard deviation
  in the pool (`"standardized"`), or Mahalanobis distance with the
  covariance of the pool or a given `trait_cov` (`"mahalanobis"`), in
  both limiting similarity and filtering; filtering by Mahalanobis
  distance is filtering by a multivariate Gaussian niche. The help
  explains when each fits: Euclidean distances count correlated traits
  more than once.
  [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md)
  and
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  take a niche with several environmental axes: a list of traits, one
  per axis, and, in `simulate_metacommunity(env = )`, a matrix with one
  column per axis; the growth rate is the product of the responses to
  each axis, as in
  [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md).

- [`traits_ou()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  simulates several correlated traits: `sigma2` can be an evolutionary
  rate matrix, or a named list with one matrix (or one rate) per
  selective regime, as in the multi-rate models of mvMORPH (Clavel et
  al. 2015). `regimes` can also be a mapped history (class `simmap`,
  such as from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)),
  so that regimes, optima and rates change within branches. With one
  number as `sigma2` the values are the same as before for the same
  seed.
  [`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md)
  builds rate matrices with the structures that `mvMORPH::mvBM()` fits:
  any matrix (`"default"`), independent traits (`"diagonal"`,
  `"equaldiagonal"`), equal rates with correlations (`"equal"`), and,
  across regimes, proportional matrices (`"proportional"`), shared
  eigenvectors (`"shared"`), equal variances (`"variance"`) or equal
  correlations (`"correlation"`). `data-raw/validation_mvtraits.R`
  checks them against mvMORPH 1.2.3: over 20000 simulations with two
  regimes changing within branches, the covariances of the tips stayed
  within 2.9 Monte Carlo standard errors of `sum_k kronecker(R_k, C_k)`
  (mvSIM: 2.8), and under OU within 3.1 of the covariance of Hansen
  (1997); `mvBM()` with the same constraint recovered the matrices (100
  tips, median of 30 fits within 9% of each true element).

- [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  draws individuals from the populations of a simulated metacommunity,
  with trait values made of the optimum of the species, a fraction of
  the gap to the environment of the site (`local_adaptation`), a
  deviation of the population that does not track the environment
  (`population_sd`) and variation among individuals
  (`intraspecific_sd`).
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  computes the T-statistics of Violle et al. (2012) per site, in base R
  for the observed values or through
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
  (Taudiere & Violle 2016, now in Suggests) for its null models; both
  give the same values. Sites are passed to cati under labels of equal
  width, because
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
  matches site names as regular expressions and would pool `site1` with
  `site10`. `data-raw/demo_tstats.R` shows that local adaptation raises
  T_IP.IC (from 0.02 to 1.0 as it goes from none to complete) and lowers
  T_IC.IR and T_PC.PR (from 0.54 and 0.89 to below 0.01), means of 20
  metacommunities.

- [`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md),
  [`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  and
  [`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md)
  implement the dynamic null model of community assembly of Pigot &
  Etienne (2015, DAMOCLES): each lineage colonises the local community
  at rate gamma and goes locally extinct at rate mu along the tree, and
  at each allopatric speciation a present ancestor leaves one present
  daughter. The likelihood follows the pruning algorithm of the article,
  the rates are fitted by maximum likelihood with
  [`optim()`](https://rdrr.io/r/stats/optim.html), and
  [`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md)
  compares MPD and MNTD of each community with communities simulated
  under the fitted rates and with random draws of the same richness. All
  is base R, from the published equations.
  `data-raw/validation_damocles.R` checks it against the DAMOCLES
  package (2.3): log-likelihoods matched `DAMOCLES_loglik()` to 3e-14 in
  180 comparisons (20 trees, three pairs of rates, three priors at the
  root, some trees with a root edge); maximum likelihoods matched
  `DAMOCLES_ML()` to 5e-4 in ten communities, ours equal or higher in
  all; and over 4000 simulated communities per setting the frequency of
  presence of each of 40 tips, the mean richness, MPD and MNTD agreed
  with `DAMOCLES_sim()` (largest difference in a tip frequency 0.032,
  \|z\| \< 3). The rates themselves are poorly identified from one
  community of 60 species (true gamma = mu = 1; estimated mu from zero
  to thousands).

- [`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md)
  counts, for each species of a simulated metacommunity, the entries and
  exits of other species (all of them, or only competitors) in the sites
  it occupies, per site and time step between recorded steps.
  [`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)
  turns this turnover into rates per branch for
  `trait_rates(edge_rates = )`, `sigma2 * (1 + k * turnover)` over the
  last `window` units of time of each tip branch, linking the assembly
  of the metacommunity to rates of trait evolution (evolution speeds up
  as the community around a lineage changes). The simulation runs among
  the tips in ecological time, so only tip branches are informed and the
  match between ecological and branch time is set by `window`.

- [`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md)
  draws the abundances, mean traits and rates of trait change from
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
  in stacked panels, with introductions marked; with `time` it draws the
  state up to that time on the axes of the whole run.
  [`animate_food_web()`](https://provetelab.org/comdynamite/reference/animate_food_web.md)
  turns these frames into a GIF (gifski in `Suggests`); the reference
  page of
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
  shows one, made by `data-raw/food_web_gif.R`.

- [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
  simulates the joint dynamics of abundances and mean traits of a
  resource, a consumer and a predator, the three-trophic-level model of
  chapter 3 of McPeek (2017, *Evolutionary Community Ecology*):
  quadratic intrinsic rates, Holling type II feeding, optional omnivory,
  and attack coefficients that are Gaussian, logistic or
  Michaelis-Menten functions of the traits. Mean traits follow the
  selection gradients, computed numerically from the fitness of a
  phenotype among individuals at the mean; everything is integrated by a
  fixed-step Runge-Kutta method in base R. Levels can be absent or
  introduced at given times, and the output keeps the selection
  gradients and the rates of trait change, for later use as rates per
  lineage.
  [`food_web_params()`](https://provetelab.org/comdynamite/reference/food_web_params.md)
  holds the parameters with the names of the book.
  `data-raw/validation_food_web.R` reproduces figure 3.2 (consumer peak
  of 10.2 at iteration 198, final traits 17.1 and 15.6, attack
  coefficient 0.091, rates of trait change 280 and 47 times lower at
  iteration 3000 than at their maxima), figure 3.3 (the consumer at a
  maximum and the resource at a minimum of their fitness surfaces at the
  end) and the effect of the predator in figure 3.11: with the predator
  the resource evolves a trait closer to its birth optimum of 12 under
  the three attack forms, from 16.8 to 11.8 with Gaussian attack (book:
  15.6 to 12.4), from 13.6 to 12.6 with logistic attack (15.4 to 12.9)
  and from 10.1 to 11.6 with Michaelis-Menten attack (10.1 to 11.5). The
  three-level Michaelis-Menten panel is reproduced closely; the consumer
  and predator traits of the logistic panel cannot be reproduced from
  the parameters of its caption, at which the consumer cannot persist.
  Against the MATLAB code of the resource-consumer models published with
  McPeek (2017, *The American Naturalist*), the fitness functions and
  gradients are the same and the trajectories of figure 3.2 agree within
  1e-4 of the author’s `ode45` integration; that code has no predator,
  and the chapter’s own code was offline and not archived. `trait_floor`
  keeps mean traits at or above a value, as the author’s code does at
  zero (the default with Michaelis-Menten attack).

- [`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md)
  turns the history of a discrete character along the tree, such as the
  region or habitat of each lineage from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  or a stochastic map, into the rate of each branch for
  `trait_rates(edge_rates = )`: a rate per state (`"state"`), a burst
  after each change that decays exponentially (`"burst"`), or a rate
  that depends on the number of lineages in the same state
  (`"sympatry"`, as in the diversity-dependent model with biogeography
  of Drury et al. 2016). The rate of a branch is its mean over time.
  `data-raw/validation_sympatry.R` checks `"sympatry"` against RPANDA:
  the covariance matrices of the tips matched those of the
  diversity-dependent model with biogeography to 3e-15 in 30 stochastic
  maps, and `RPANDA::fit_t_comp()` recovered r from traits simulated
  with these rates (60 tips, three regions; mean estimate -0.30 for a
  true -0.3, 80% of estimates between -0.32 and -0.27, and 0.005 for a
  true 0; mean sigma2 0.98 and 0.91 for a true 1).

- [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md):
  `shift_rate` draws shifts as a Poisson process along the branches, as
  under the prior of BAMM, so they can fall within branches and the rate
  of a branch is its mean over time; `relative = FALSE` gives each new
  regime a rate of its own, as in BAMM and bifrost, instead of
  multiplying the rate before the shift; `shift_size` can be a function
  that draws the multipliers. `data-raw/validation_bamm.R` runs BAMM
  2.5.0 on 20 trees of 100 tips with about seven such shifts each (10
  million generations, lowest effective sample size 62): the mean
  posterior rate per branch correlated with the true rate (Spearman
  0.51, Pearson 0.67 on the log scale), and the marginal probability of
  a shift was 0.18 on branches with a true shift against 0.02 on the
  others.

- [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md)
  simulates a continuous trait whose rate varies among branches and
  keeps the rates, as the known truth for methods that estimate rates
  per branch. Three models: shifts at a few branches, inherited by their
  descendants or not (`"shift"`); a logarithm of the rate that evolves
  by Brownian motion (`"gbm"`); and independent lognormal rates
  (`"iid"`). It gives the trait at tips and nodes, the rate of each
  branch and the rescaled tree.

- [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md):
  `alpha` multiplies the rates of each branch by a gamma factor with
  mean 1, the model of
  [`phytools::fitgammaMk()`](https://rdrr.io/pkg/phytools/man/fitMk.html),
  or `edge_rates` gives the factors; `rate_classes` adds hidden rate
  classes, as in corHMM, and keeps the hidden history in `$hidden`.
  Without them the results are the same as before with the same seed.

- `data-raw/validation_rates.R` compares the truth with RRphylo,
  [`phytools::multirateBM()`](https://rdrr.io/pkg/phytools/man/multirateBM.html),
  traitMedusa (motmot), `fitgammaMk()` and corHMM. With one change per
  branch, rates per branch are weakly identifiable: the Spearman
  correlation between the true change and the true rate of a branch was
  0.29 to 0.40, and RRphylo rates followed the change (0.32 to 0.38) but
  hardly the rate (0.07 to 0.13). multirateBM recovered the log rates at
  nodes with r = 0.62 and slope 0.46. traitMedusa found a true eightfold
  shift in 20 of 30 trees, with a median estimated multiplier of 8.3.
  fitgammaMk gave a median alpha of 0.33 for a true 0.3 (500 tips), and
  corHMM preferred two hidden classes in 10 of 12 trees with classes and
  in none of 12 without.

- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  `immigration = "multinomial"` distributes the emigrants of each site
  among the others by the dispersal kernel and conserves them, as
  mcomsimr does, on spatially explicit landscapes. The default stays
  `"poisson"`, as in the code of Thompson et al. (2020). With
  extirpation 0.02 and dispersal 0.001 on a grid, mean local richness
  was 1.51 with Poisson, 1.61 with multinomial and 1.61 in mcomsimr.

- purrr moved from `Imports` to `Suggests`: it only gives the progress
  bar of long loops in interactive sessions. The package now imports
  only ape and base packages.

## comdynamite 0.6.3

- [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md):
  `optima`, `sigma` and `r_max` take values from elsewhere, such as
  mvMORPH, by species name or in the order of the species. The traits
  not given are drawn as before, and with the same seed they are the
  same as without the external ones. `rescale_optima = FALSE` keeps
  external optima on their own scale.
- Herodotools:
  [`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md)
  gives its inputs from a simulation or an `mc_data` object, with the
  true ancestral areas from a history of
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md);
  [`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md)
  limits each species to the sites of its region;
  [`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md)
  gives the true ages at which lineages arrived in each region. The
  article “Known biogeographic histories for Herodotools” shows the
  steps and the results of `data-raw/validation_herodotools.R`.
- Trait histories along the tree, to validate reconstructions:
  `trait_bm(nodes = TRUE)` and `species_traits(nodes = TRUE)` give the
  values at the internal nodes, drawn after everything else so that the
  other values stay the same with the same seed;
  [`traits_ou()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md)
  keeps them as the attribute `nodes`.
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  simulates a discrete trait under the Mk model and keeps every change,
  in the `simmap` format of phytools.
  [`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md)
  gives the tip-based metrics of Luza et al. (2021), as computed by
  Herodotools, from a known history, plus the number of changes and the
  time since the last one in continuous time.
- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  `pool` limits each species to the sites where it can occur, such as
  its range or a biogeographic reconstruction, and `initial` gives the
  abundances or presences at the start. Without them, results with the
  same seed are the same as before.
- [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md):
  new statistic `trait_k`, Blomberg’s K of the trait over the species
  shared by `trait` and `tree`, as
  [`picante::Kcalc()`](https://rdrr.io/pkg/picante/man/Kcalc.html) gives
  it. It is `NA` without a tree or a trait.

## comdynamite 0.6.2

- [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  and
  [`as_mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md):
  community, tree, traits, environment, coordinates, areas and time
  series in one object, with names checked once. The analyses and
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
  take an `mc_data` object as first argument and fill what the call
  leaves out. Simulations keep the tree their traits came from, as
  `tree`.
- [`run_scenarios()`](https://provetelab.org/comdynamite/reference/run_scenarios.md):
  grids of scenarios with replicates paired by what `setup` draws, a
  seed per cell and optional parallel runs.
- [`equilibrium_check()`](https://provetelab.org/comdynamite/reference/equilibrium_check.md):
  warns when a simulation is still changing at the end of its record.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for
  `mc_sim`: long table from `history`.
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  with `site_permutations`: each set of predictors now moves alone in
  its tests, with the response and the other set in place: `env` for
  `ab`, `a` and its semipartials, `space` for `bc`, `c` and its
  semipartials; `abc` still moves both. Up to 0.6.1 the moves took both
  sets together in every test, which for `c` compared space with an
  environment that no longer explained anything: on a grid it rejected
  22% of data sets with a strong environment without structure and never
  with a structured one. In `data-raw/validation_fractions_null.R`
  (choice on one generator, check on another), moving `env` alone
  rejected `a` in 3.7% to 7.1% of data sets where the environment had no
  effect, where free permutations reached 20% with an environment
  structured outside `space`, and moving `space` alone rejected `c` in
  3.0% to 8.0% of null data sets, where free permutations reached 24%
  with a strong, structured environment. When the environment acts on a
  trait other than the one analysed, the test of sites of `a` still
  rejects, in 45% to 75%, as with free permutations: read `p_max`, which
  stayed at 2% to 4%. The power of that test of `c` was complete for a
  field in the span of `space` and 25% to 28% for a field that `space`
  holds only in part, since the maps moved around the grid explain such
  a field about as well: it asks whether the maps in `space` matter
  beyond the same maps elsewhere, not whether the remainder has any
  spatial structure. The result records `permutations$moved` and the
  print says which set moves. Free permutations are unchanged.
- Permutations of
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  run in parallel under a plan of the future package, which is optional.
- Install from r-universe, as the README now shows.
- [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md):
  the same statistics for simulated and observed data, covering
  diversity, occupancy, abundance distribution, SES MPD and MNTD,
  distance decay, trait-environment relation and, with more than one
  step, temporal variation.
- [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md):
  principal coordinates of phylogenetic structure of Duarte et
  al. (2016), the method of the package PCPS, with tests of all the
  dissimilarity and of each axis, and `site_permutations` for a
  restricted null model, which that package does not have. Its
  statistics match those of PCPS 1.0.9 (`data-raw/validation_pcps.R`).
  With autocorrelated fields unrelated to each other, free permutations
  of sites rejected 29-80% at the 5% level; the torus rejected 4-7% when
  fields wrap around the grid and 6-25% when they stop at its edges.
- [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md):
  rectangular grids, with `side` of length two, and new `symmetries` and
  `n_perm`.

## comdynamite 0.6.1

- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  the tests of fractions `a` and `c`, and of the semipartials, now use
  the share of the component that the other set leaves unexplained as
  statistic, which is monotone in the partial F, instead of the fraction
  itself, which the original code tests. The permutations are the same.
  The old test lost its power when the other set explained much: in
  `data-raw/validation_fractions.R` it rejected `a` in none of 500 data
  sets where the environment had no effect and space a strong one (new:
  5.6% to 5.8%), and in 15% where the environment had a weak effect
  beside space (new: 100%). `fractions` has a new column `statistic`,
  and p-values are `NA` when a fraction cannot be estimated. See the
  validation article for the limits: a free permutation of sites rejects
  `a` too often when `env` has structure at scales outside `space`, and
  `site_permutations` are nearly degenerate for `c`.
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  test of the shared fraction `b` with `shared_null`, by Moran spectral
  randomisation of `env` from `coords` (adespatial) or by moves of `env`
  alone with `site_permutations`, as in
  [`adespatial::envspace.test()`](http://adeverse.github.io/adespatial/reference/envspace.test.md).
  `b` is tested only when `env` is related to `space` (`level`), and the
  result has `shared`.
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md)
  colours `b` by its test when there is one.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  the results of the analyses.
- [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md):
  restricted null model for the events of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  that moves the whole event by a rotation or a reflection and a
  translation, so it keeps the shape of the sides and the distance
  between them. It is better than
  [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
  for sides with room to move, and gives few or no placements for events
  that span the landscape. Its `independent` argument covers those at
  the price of an assumption on the range of autocorrelation, and
  remains liberal in some layouts. See its help for the rates.
- [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
  off a regular lattice now enumerates the moves, on a grid of rotations
  and reflections with a translation to each site, instead of drawing
  them at random, so which placements exist no longer changes between
  runs. In `data-raw/validation_offgrid.R` it rejected 4.7% and 5.9% of
  the tests at the 5% level for adjacent areas off a lattice, and 6.8%
  to 8.0% for two blocks placed at random, where it is liberal.
- [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
  and
  [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md):
  the default `gap` follows the fourth nearest site, the step of a full
  grid and the spacing of sparse sites, instead of the nearest one.
  [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
  thus no longer stops on a grid with holes, and the holes between sites
  count as inside the landscape while the rim beyond an edge does not.
  The default of
  [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
  is unchanged on a full grid where most sites are off the edge, and
  larger with sparse or irregular sites, where it changes the patches
  and the results: with 25% to 40% of the sites of a square and areas in
  opposite corners it rejected 5.6% to 5.9%, against 15% to 17% before
  (`data-raw/validation_gap.R`).
- [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
  draws the placements at random when an event has more than `n_perm` of
  them. Before, with fewer than twice `n_perm`, it took the first ones
  in the order of enumeration.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `permutations$areas` now gives the number of permutations of each
  event, which may differ, and
  [`print()`](https://rdrr.io/r/base/print.html) shows them.
- [`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md):
  community matrix with sites ordered by
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  and species by the tree, as Peres-Neto & Kembel (2015) suggest.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `r2_weighted`, which stays between 0 and 1, unlike `r2`, inherited
  from the original code; `n_draws`, the valid draws of sites. Errors in
  parallel workers now stop the function. The default of `n_perm` is
  999, as in the other functions.
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  `site_permutations`.
- [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  accept empty sites with `site_permutations`, since these sites are
  part of the landscape.
- The significance level is now `level` in
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md),
  [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md),
  [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)
  and
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md),
  because `alpha` also names the competition matrix. `alpha` still
  works, with a warning.
- Names of results, which break code that reads them: `p_taxa` is now
  `p_species`, `p_partial_taxa` is `p_partial_species`, and `p_partial`,
  the larger of two p-values, is `p_partial_max`, as in `p_max`. `r2` of
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  is `r2_mean` and that of
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  is `r2_adjusted`, which say what they hold. In
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md),
  `r2` is unchanged. `statistic` stays, as the test statistic of each
  function.
- [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md):
  `predictor` is now `env`, as in the other functions, and
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md):
  `test = "taxa"` is now `"species"`. The old forms work, with a
  warning.
- [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md):
  evaluated with autocorrelation up to the scale of the landscape, and
  advice on `gap` for sparse sites. The rates given in 0.6.0 averaged
  over random places of areas. With fixed areas the test is liberal for
  sides far apart (6% to 10% at the 5% level) and conservative for
  neighbouring sides, and the help now says so.
- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  matches `alpha` to `traits` by name. Before, a matrix with species in
  another order gave wrong results without warning.
- [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
  [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)
  and
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  accept traits as a matrix or a tibble, and species without names are
  called sp1, sp2, …
  [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)
  matches sets of traits by name.
  [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)
  warns when it invents names.
- Predictors must vary among occupied sites with `site_permutations`.
- Help: weights of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  and what the tests of fractions of
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  ask.
- [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  gives an `r2` of 0, not `NaN`, when no adjusted R2 is positive.
- [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
  accepts clades with branches of length zero, and
  [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)
  accepts a predictor and an event with the same name.
- Clearer errors in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md)
  with one site,
  [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md)
  and
  [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)
  without tested nodes.

## comdynamite 0.6.0

- [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md):
  restricted null model for the events of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  in landscapes with edges, which moves the patches of each area and
  keeps their number and size. It rejected 3.0% to 6.3% of tests at the
  5% level in `data-raw/validation_edges.R`, against 44% to 55% with
  free permutations.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `area_permutations` takes one matrix of permutations per event.
- New title and description, since the package now covers analyses as
  well as simulations.
- [`?comdynamite`](https://provetelab.org/comdynamite/reference/comdynamite-package.md)
  and the README map questions to functions, and the reference and
  articles on the website are grouped by question.
- New articles: “Preparing your data” and “Testing the whole tree”, on
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  and
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md).

## comdynamite 0.5.2

- [`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md):
  events of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  from the area of each site, as pairs of areas stated in advance or as
  each area against the others.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `site_permutations` gives a restricted null model. For events, it
  moves the areas across the sites.
- [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md):
  biogeographic events are triangles, as in Leibold & Chase (2018).
- [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md):
  two overlapping bars along the total variation, as in Figure 1 of
  Desdevises et al. (2003), so that the shared fraction is seen as
  shared. New `labels`.
- Progress bars with an estimate of the time left in
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md),
  in interactive sessions.
- Examples, README and vignettes use trees of treeshapesim, in Suggests.
  The results of
  [`vignette("phylogenetic-structure")`](https://provetelab.org/comdynamite/articles/phylogenetic-structure.md)
  and
  [`vignette("evaluating-methods")`](https://provetelab.org/comdynamite/articles/evaluating-methods.md)
  were run again.
- The help of functions adapted from the MATLAB code of the supplements
  says so.

Corrections that change results:

- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  draws are now ranked on F statistics. The statistic of each predictor,
  taken from the original code, depended on the scale of the response,
  and neither statistic took into account that draws of species occupy
  different numbers of sites. `p_partial`, `p_partial_joint`, `p_taxa`,
  `p_max` and `p_joint` differ from those of 0.5.1, and so do the same
  p-values of events. `statistic`, `p_site` and `p_partial_site` of
  predictors do not.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  a node with a clade without individuals is no longer tested, as it was
  a repeat of another node.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  with `cores` above 1, each call now uses its own random streams.
- [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  and
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md):
  species without individuals leave the analysis. In
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  they could turn numerical noise into an axis.
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  species without individuals leave the analysis, which could stop it.

Checks of data:

- All analysis functions match sites by name when the data have row
  names, and stop if the names differ. Before, misordered predictors
  gave a wrong result without warning.
- Missing values, predictors without variance and columns that are not
  numeric now stop with a message that names them.
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  gave p = 1 in these cases.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  checks `biogeography` before the permutations and explains a tree in
  `nodes`.
- Results keep the numbers of permutations and the null model, and
  [`print()`](https://rdrr.io/r/base/print.html) shows them with the
  smallest possible p-value. Nodes that were not tested are reported,
  and
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  says when no vector was retained.

## comdynamite 0.5.1

- [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md):
  `branch_lengths = TRUE` gives the coding with branch lengths of the
  supplement of Leibold et al. (2010).
- [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md):
  with `type = "fan"`, the names of the variables follow the radii and
  no longer overlap the names of the species.
- [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md):
  semipartial correlations of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  as a table of nodes by variables.
- [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md)
  now uses the values of the code of Ives & Helmus (2011): the default
  of `intercept_sd` is 0.5, and slopes have a standard deviation of
  `slope_sd` among species.
- [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md):
  in scenario 3, the error is added to `y` with unit variance and only
  in `y_test`, and in scenarios 5 and 6 `x` is in random order along the
  transect.
- New article on analysing clades and traits. The article on validation
  now compares Table 5 of Peres-Neto et al. (2012).

## comdynamite 0.5.0

- [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md):
  static metacommunities from the Gaussian response of Dray & Legendre
  (2008), also used by Peres-Neto et al. (2012), as a baseline without
  population dynamics. Data can be Poisson counts, binomial presences or
  the expected values. With several gradients, for example environment
  and space, the expected abundance is the product of the responses.
- [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md)
  and
  [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md):
  presences and absences from logistic responses to two gradients, with
  the seven scenarios of Appendix C of Peres-Neto et al. (2012).
- [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  phylogenetic gradient analysis of Peres-Neto et al. (2012), tested
  against the original code in Octave. The test between a trait and the
  environment stays in
  [`ade4::fourthcorner()`](https://adeverse.github.io/ade4/reference/fourthcorner.html).
- Functions that take a community and a phylogeny, traits or contrasts
  now match species by name and, as
  [`picante::match.phylo.comm()`](https://rdrr.io/pkg/picante/man/data.checking.html),
  drop with a warning the species found in only one of them. Repeated
  names are an error.
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  semipartial correlations of each predictor, with their tests,
  following Appendix A of Peres-Neto et al. (2012).
- [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md):
  results of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  on the tree, as in Figure 3 of Leibold et al. (2010), as a phylogram
  or a fan.
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md):
  fractions of
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md).
- `daphniids`: data of Leibold et al. (2010), rebuilt from the public
  survey.
- [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
  now keeps nodes that join two tips, as the article does.
  `cherries = FALSE` gives the original code.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `semipartial` replaces `partial_r2`, following the article, and there
  are `p_adjusted` (false discovery rate), biogeographic events, the
  joint permutation of the article (`p_joint`), 9999 permutations by
  default and `cores`.
- [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  and
  [`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md):
  phylogenetic community composition of Peres-Neto & Kembel (2015), with
  site and species permutations.
- [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  and
  [`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md):
  phylogenetic beta-diversity framework of the same article, a
  distance-based redundancy analysis of the mean phylogenetic distance
  between communities.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md):
  regression of clades on the environment, node by node, from Leibold et
  al. (2010), with site and species permutations.
- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  stochastic extirpation of populations (`extirpation`), used by
  Thompson et al. (2020) in their competition-colonisation scenario.
- New article on the validation of the model, with the differences
  between the article of Thompson et al. (2020) and the code archived
  with it.
- New article with animations of metacommunity dynamics.

## comdynamite 0.1.0

- First version: the model of Thompson et al. (2020) with phylogenetic
  traits and competition, spatially implicit and grid landscapes, torus
  permutations and time-varying environments. Code from the pcps-sim
  project.

- Rule-based assembly
  ([`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)):
  limiting similarity (GAUSE algorithm of Colwell & Winkler 1984),
  habitat filtering and neutral assembly, as in Kraft et al. (2007).

- [`traits_ou()`](https://provetelab.org/comdynamite/reference/comdynamite-deprecated.md):
  traits under an Ornstein-Uhlenbeck process with one or several optima
  (multi-peak OU; Brownian motion when `alpha = 0`), with
  [`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md)
  to assign regimes to clades.

- New vignettes:
  [`vignette("phylogenetic-structure")`](https://provetelab.org/comdynamite/articles/phylogenetic-structure.md),
  on how dispersal, niche conservatism and competition shape
  phylogenetic community structure, and
  [`vignette("evaluating-methods")`](https://provetelab.org/comdynamite/articles/evaluating-methods.md),
  on type I error and power with simulations.

- [`vignette("evaluating-methods")`](https://provetelab.org/comdynamite/articles/evaluating-methods.md)
  shows tree shape as a factor of the design, with trees from
  treeshapesim (in Suggests).
