# Recalibrate the credible intervals of ABC with random forests

Measures, on simulations that were not used to fit, where the true value
of each parameter falls in its estimated posterior, and uses this to
correct the credible intervals of
[predict()](https://provetelab.org/comdynamite/reference/abc_fit.md).
For each simulation of `pseudo`, the probability integral transform
(PIT) is the posterior distribution function at the true value, read
from 999 quantiles of the quantile regression forest. If the intervals
were calibrated, the PIT would be uniform between 0 and 1 (Prangle et
al. 2014); intervals that are too wide give PIT values crowded in the
middle, and intervals that are too narrow give values crowded at the
ends. With `calibration`,
[predict()](https://provetelab.org/comdynamite/reference/abc_fit.md)
then asks the forest for the quantiles of the PIT instead of the nominal
probabilities: for a `level` of 0.9, the 5% and 95% quantiles of the PIT
values, so that 90% of the calibration simulations would fall inside.

## Usage

``` r
abc_calibrate(fit, pseudo, ncores = 1)
```

## Arguments

- fit:

  An `mc_abc_fit` from
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).

- pseudo:

  Data frame of simulations from the same prior and model that were not
  used to fit, with the parameters and statistics of `fit` (for example
  a second
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
  with another seed). Rows with missing or infinite values are dropped;
  at least 20 must remain. A few hundred give quantiles of the PIT
  precise enough for a 90% interval.

- ncores:

  Number of threads for the forests. One whole number, 1 or more.

## Value

An `mc_abc_calibration`: a list with `pit` (one vector per parameter),
`coverage` (a data frame with the fraction of simulations of `pseudo`
inside the nominal 50%, 80%, 90% and 95% intervals of each parameter,
before recalibration) and `n`, the number of simulations used.

## Details

This is a global version of the recalibration of Rodrigues et al.
(2018), who weight the calibration simulations by their closeness to the
observation. It corrects the coverage on average over the prior, as
measured by
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
and not the interval of one observation, which can still be too wide or
too narrow where the posterior differs from the average. It cannot make
an uninformative posterior informative either: for a parameter that the
statistics do not inform, the PIT is already about uniform and nothing
changes. Nor can it widen an interval beyond the parameters of the
simulations that share leaves with the observation: with one statistic
and the default `min_node_size`, the truth fell outside all of them in
about 13% of the simulations, so recalibrated 90% and 95% intervals both
held it in about 87% (`data-raw/abc_recalibration.R`); a larger
`min_node_size` of
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
is the remedy there.

## References

Prangle, D., Blum, M. G. B., Popovic, G., & Sisson, S. A. (2014).
Diagnostic tools for approximate Bayesian computation using the coverage
property. *Australian & New Zealand Journal of Statistics*, 56, 309–329.
[doi:10.1111/anzs.12087](https://doi.org/10.1111/anzs.12087)

Rodrigues, G. S., Prangle, D., & Sisson, S. A. (2018). Recalibration: A
post-processing method for approximate Bayesian computation.
*Computational Statistics & Data Analysis*, 126, 53–66.
[doi:10.1016/j.csda.2018.04.004](https://doi.org/10.1016/j.csda.2018.04.004)

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
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
make <- function(n) {
  d <- data.frame(a = runif(n))
  d$s <- d$a + rnorm(n, 0, 0.05)
  d
}
fit <- abc_fit(make(1000), parameters = "a", ntree = 100)
#> abc_fit() implements the method of:
#>   Raynal, L., Marin, J.-M., Pudlo, P., Ribatet, M., Robert, C. P., & Estoup, A. (2019). ABC random forests for Bayesian parameter inference. Bioinformatics, 35, 1720-1728. https://doi.org/10.1093/bioinformatics/bty867
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
calibration <- abc_calibrate(fit, make(200))
calibration$coverage
#>   parameter level coverage
#> 1         a  0.50    0.375
#> 2         a  0.80    0.600
#> 3         a  0.90    0.650
#> 4         a  0.95    0.730
predict(fit, c(s = 0.5), calibration = calibration)
#>   observation parameter  estimate       low      high      mean
#> 1           1         a 0.4934904 0.4349082 0.5394022 0.4883662
```
