# How many simulations does an ABC posterior rest on?

The posterior of
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
for an observation is a weighted sample of the reference table, with
weights from each random forest. For each observation and parameter this
gives the effective sample size of the weights (Kish) and the Pareto
k-hat of their upper tail (Vehtari et al. 2024), computed by
[`posterior::pareto_khat()`](https://mc-stan.org/posterior/reference/pareto_khat.html)
on the positive weights, and for each observation the statistics that
fall outside the central `coverage` of the reference table.

## Usage

``` r
abc_diagnose(fit, newdata, coverage = 0.99)
```

## Arguments

- fit:

  An `mc_abc_fit`.

- newdata:

  Named numeric vector of observed statistics or a data frame with one
  observation per row.

- coverage:

  Central probability of the reference table within which a statistic
  counts as inside.

## Value

A data frame with `observation`, `parameter`, `n_positive` (simulations
with positive weight), `ess`, `khat`, `khat_threshold`, `n_outside`
(statistics of the observation outside the central `coverage` of the
table) and `outside` (their names).

## Details

The weights of random forests are not importance ratios, so the
thresholds of k-hat, derived for importance sampling, do not carry over.
In our checks (`data-raw/abc_diagnose.R`), k-hat and the effective
sample size did not predict which estimates were wrong within the model,
while data from another model had about half the effective sample size
and more statistics outside the table. Read them as warnings that the
observation is unlike the simulations, not as measures of the precision
of an estimate. An observation outside the range of the simulated
statistics is the clearer warning, since the forests then extrapolate
from the edge of the table without saying so.

## References

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72), 1–58.

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
table <- data.frame(a = runif(500))
table$s1 <- table$a + rnorm(500, 0, 0.05)
fit <- abc_fit(table, "a", "s1", ntree = 100)
abc_diagnose(fit, data.frame(s1 = c(0.5, 3)))
#> Warning: Number of tail draws cannot be less than 5. Changing to 5.
#> Warning: Number of tail draws cannot be less than 5. Changing to 5.
#>   observation parameter n_positive      ess khat khat_threshold n_outside
#> 1           1         a         13 5.532480   NA             NA         0
#> 2           2         a         11 4.602561   NA             NA         1
#>   outside
#> 1        
#> 2      s1
```
