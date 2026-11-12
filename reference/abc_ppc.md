# Posterior predictive check after ABC

Draws parameter sets from the approximate posterior of an
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
and compares the observed statistics with those of the data the model
predicts. Statistics that did not enter the fit make the check stronger,
since the fit was tuned to reproduce the others (Gelman et al. 2013). A
model can give good estimates and still fail this check (Pontarp et al.
2019).

## Usage

``` r
abc_ppc(
  fit,
  observed,
  table,
  statistics = NULL,
  n = 200,
  simulate = NULL,
  setup = NULL,
  summarise = function(sim, shared) summary_stats(sim, shared$tree),
  seed = 1
)
```

## Arguments

- fit:

  An `mc_abc_fit`.

- observed:

  Named numeric vector of observed statistics, with those of `fit` and
  those to check.

- table:

  The reference table `fit` was fitted on, needed without `simulate` and
  to draw the parameters.

- statistics:

  Names of the statistics to check. By default those of `table` that did
  not enter `fit`, or all if every statistic did.

- n:

  Number of draws from the posterior.

- simulate, setup, summarise:

  Optional functions to simulate new data, as in
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md).

- seed:

  Integer that sets all the seeds.

## Value

A data frame with one row per statistic: `observed`, the median and the
5% and 95% quantiles of the predicted statistic, and `p`, the two-sided
posterior predictive p-value (twice the smaller tail). Small p-values
point to statistics the model does not reproduce. The attribute `draws`
holds the drawn parameters and predicted statistics.

## Details

Random forests give one set of weights over the reference table per
parameter, not a joint posterior. Whole rows of the table are drawn with
the mean of these weights, which keeps the correlation among parameters
within a row but is an approximation of the joint posterior. Without
`simulate`, the statistics of the drawn rows are used, which costs
nothing; with `simulate`, new data are simulated from the drawn
parameters.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian data analysis* (3rd ed.). Chapman &
Hall/CRC.

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
table <- data.frame(a = runif(400))
table$s1 <- table$a + rnorm(400, 0, 0.05)
table$s2 <- table$a^2 + rnorm(400, 0, 0.05)
fit <- abc_fit(table, "a", "s1", ntree = 100)
abc_ppc(fit, c(s1 = 0.5, s2 = 0.25), table, n = 100)
#>   statistic observed    median        low      high    p   n
#> 1        s2     0.25 0.2600238 0.08515926 0.3168431 0.66 100
```
