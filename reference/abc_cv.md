# Cross-validation of an ABC model across sites

For each fold, fits
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
on the statistics of the sites kept, draws from the approximate
posterior given the observed statistics of those sites, and predicts the
statistics of the sites left out from the same rows of the reference
table, as in
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md).
A model that fits the sites it saw but not the others gives observed
test statistics outside the predictive intervals. Use spatial blocks
from
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)
so that the test sites are not neighbours of the training sites.

## Usage

``` r
abc_cv(
  table,
  observed,
  parameters = NULL,
  statistics = NULL,
  level = 0.9,
  n = 500,
  ntree = 500,
  ncores = 1,
  seed = 1
)
```

## Arguments

- table:

  An `mc_abc_table` from
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
  with `summarise = fold_summarise(...)`.

- observed:

  Named numeric vector with the same statistics for the observed data,
  from the same
  [`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md)
  function, or a data frame with one data set per row. The forests are
  fitted once per fold for all rows.

- parameters:

  Names of the parameters. By default those of `table`.

- statistics:

  Names of the statistics without the prefix, such as `"alpha"`. By
  default all.

- level:

  Probability of the predictive interval.

- n:

  Number of draws from the posterior of each fold.

- ntree:

  Number of trees per forest.

- ncores:

  Number of threads for the forests.

- seed:

  Integer that sets the seed.

## Value

A data frame with one row per observed data set, fold and statistic:
`observed`, the median and the interval of the predicted statistic,
`inside` (whether the observed value is in the interval), `p` (two-sided
posterior predictive p-value) and `z`, the observed value minus the
predicted median over the predictive standard deviation.
[`summary()`](https://rdrr.io/r/base/summary.html) averages over folds
and data sets (`n` rows per statistic): the coverage, the mean of `z`,
the root mean square of `z` and the fraction of p-values below 0.05.

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
table <- data.frame(a = runif(400))
for (f in 1:2) {
  table[[paste0("f", f, "_train_s")]] <- table$a + rnorm(400, 0, 0.05)
  table[[paste0("f", f, "_test_s")]] <- table$a + rnorm(400, 0, 0.1)
}
observed <- c(f1_train_s = 0.5, f1_test_s = 0.55, f2_train_s = 0.5, f2_test_s = 0.9)
cv <- abc_cv(table, observed, parameters = "a", ntree = 100)
cv
#>   observation fold statistic observed    median       low      high inside   p
#> 1           1    1         s     0.55 0.5007813 0.2110433 0.6163250   TRUE 0.7
#> 2           1    2         s     0.90 0.5203422 0.4107620 0.5585261  FALSE 0.0
#>           z
#> 1 0.4384356
#> 2 6.7852324
summary(cv)
#>   statistic n coverage   z_mean   z_rms rejected
#> 1         s 2      0.5 3.611834 4.80789      0.5
```
