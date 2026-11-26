# Estimate parameters by ABC with random forests

Fits one regression random forest per parameter, with the statistics as
predictors, following Raynal et al. (2019) and the package abcrf. Random
forests need no tolerance and no choice of statistics, cope with many
correlated statistics and rank their importance.
[`predict()`](https://rdrr.io/r/stats/predict.html) gives the posterior
median and a credible interval from quantile regression forests
(Meinshausen 2006). Statistics with more than 5% of missing or infinite
values and constant statistics are dropped, and then simulations with
any missing statistic.

## Usage

``` r
abc_fit(
  table,
  parameters = NULL,
  statistics = NULL,
  ntree = 500,
  ncores = 1,
  min_node_size = 5
)

# S3 method for class 'mc_abc_fit'
predict(object, newdata, level = 0.9, ncores = 1, calibration = NULL, ...)

abc_importance(x)
```

## Arguments

- table:

  An `mc_abc_table` from
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
  or a data frame with the parameters and the statistics.

- parameters:

  Names of the parameters to estimate. By default those of `table`.

- statistics:

  Names of the statistics to use. By default all the others.

- ntree:

  Number of trees per forest. One whole number, 1 or more.

- ncores:

  Number of threads for the forests. One whole number, 1 or more.

- min_node_size:

  Smallest number of simulations in a node for it to be split, the
  `min.node.size` of
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).
  One whole number, 1 or more. The default, 5, is that of Raynal et al.
  (2019); larger values give larger leaves, which widen intervals that
  [`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
  shows to be too narrow, at the cost of more smoothing.

- object:

  An `mc_abc_fit`.

- newdata:

  Named numeric vector of observed statistics, as from
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md),
  or a data frame with one observation per row.

- level:

  Probability of the credible interval. One number greater than 0 and
  less than 1.

- calibration:

  `NULL` (the default), or an `mc_abc_calibration` from
  [`abc_calibrate()`](https://provetelab.org/comdynamite/reference/abc_calibrate.md)
  for the same fit, whose quantiles of the probability integral
  transform replace the nominal probabilities of `low` and `high`.

- ...:

  Not used.

- x:

  An `mc_abc_fit`.

## Value

An `mc_abc_fit`, a list with `forests` (one per parameter), `training`,
`parameters`, `statistics` and `prior`.

For [`predict()`](https://rdrr.io/r/stats/predict.html), a data frame
with `observation` (the row of `newdata`; with row names, also `name`),
`parameter`, `estimate` (posterior median), `low`, `high` and `mean`.
The interval from `low` to `high` is approximate; see Details.

## Details

Check with
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
that a parameter can be recovered from simulated data before trusting
its estimate from real data. A good estimate does not validate the model
either.

The credible intervals are approximate and are not calibrated by
construction. The forest gives the parameters of the simulations that
share leaves with the observation, a sample of the posterior given
statistics near the observed ones and not given the observed statistics
themselves. With several statistics the leaves are wide along each of
them, the posterior is smoothed and the intervals are wider than the
exact ones: Raynal et al. (2019) report fatter tails and slightly
conservative intervals, and in our extended tests (10 statistics, two
well identified parameters) the 90% intervals held the truth in 92% to
99% of 200 pseudo-observations. With one or two statistics the leaves
are narrow and the sample is small instead: with one statistic and the
default `min_node_size` about 5 simulations carry the weight, and the
90% interval of a parameter with a known posterior held the truth in
about 65% of the cases, whatever the number of simulations. The
consistency of quantile regression forests asks for leaves that grow
with the number of simulations (Meinshausen 2006); there, with 5000
simulations, a `min_node_size` of 50, 200 and 500 brought the coverage
to 0.85, 0.89 and 0.91.
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
measures the coverage that the intervals have for a given table,
[`abc_calibrate()`](https://provetelab.org/comdynamite/reference/abc_calibrate.md)
corrects it on average over the prior, and
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md)
gives the effective number of simulations behind a posterior.

## Functions

- `predict(mc_abc_fit)`: Posterior median and credible interval for each
  observation and parameter.

- `abc_importance()`: Importance of each statistic for each parameter,
  as the decrease in impurity scaled to sum to one.

## References

Meinshausen, N. (2006). Quantile regression forests. *Journal of Machine
Learning Research*, 7, 983–999.

Raynal, L., Marin, J.-M., Pudlo, P., Ribatet, M., Robert, C. P., &
Estoup, A. (2019). ABC random forests for Bayesian parameter inference.
*Bioinformatics*, 35, 1720–1728.
[doi:10.1093/bioinformatics/bty867](https://doi.org/10.1093/bioinformatics/bty867)

## See also

Other inference:
[`abc_calibrate()`](https://provetelab.org/comdynamite/reference/abc_calibrate.md),
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_sequential()`](https://provetelab.org/comdynamite/reference/abc_sequential.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`nobs.mc_node_test()`](https://provetelab.org/comdynamite/reference/nobs.mc_node_test.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`plot.mc_abc_choice()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_choice.md),
[`plot.mc_abc_fit()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_fit.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
table <- data.frame(a = runif(300), b = runif(300))
table$s1 <- table$a + rnorm(300, 0, 0.05)
table$s2 <- table$b * table$a + rnorm(300, 0, 0.05)
fit <- abc_fit(table, parameters = c("a", "b"), ntree = 100)
predict(fit, c(s1 = 0.5, s2 = 0.2))
#>   observation parameter  estimate       low      high      mean
#> 1           1         a 0.5083821 0.3529379 0.5872987 0.5063800
#> 2           1         b 0.4057038 0.2934010 0.7191584 0.4410285
```
