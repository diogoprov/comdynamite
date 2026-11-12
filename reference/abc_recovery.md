# Can the parameters be recovered by ABC?

Leaves `n_test` simulations out of the reference table as
pseudo-observations with known parameters, fits
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
on the rest and compares the estimates with the truth. A parameter that
is not recovered from simulated data should not be inferred from real
data (Pontarp et al. 2019).

## Usage

``` r
abc_recovery(
  table,
  parameters = NULL,
  statistics = NULL,
  n_test = 200,
  level = 0.9,
  ntree = 500,
  ncores = 1,
  seed = 1
)
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

- n_test:

  Number of pseudo-observations.

- level:

  Probability of the credible interval.

- ntree:

  Number of trees per forest.

- ncores:

  Number of threads for the forests.

- seed:

  Integer for the choice of the pseudo-observations.

## Value

A data frame with one row per parameter: `n`, `r2` (1 minus the squared
error over the variance of the truth), `bias` and `width` (mean width of
the interval) as fractions of the prior range, and `coverage`, the
fraction of intervals that contain the truth, with its standard error.
The attribute `estimates` holds every estimate. A coverage near `level`
says the intervals are calibrated; R2 says how much the data narrow the
prior.

## References

Pontarp, M., Brännström, Å., & Petchey, O. L. (2019). Inferring
community assembly processes from macroscopic patterns using dynamic
eco-evolutionary models and Approximate Bayesian Computation (ABC).
*Methods in Ecology and Evolution*, 10, 450–460.
[doi:10.1111/2041-210X.13129](https://doi.org/10.1111/2041-210X.13129)

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
table <- data.frame(a = runif(400), b = runif(400))
table$s1 <- table$a + rnorm(400, 0, 0.05)
table$s2 <- rnorm(400)
abc_recovery(table, parameters = c("a", "b"), n_test = 50, ntree = 100)
#>   parameter  n         r2        bias coverage coverage_se    width
#> 1         a 50  0.9523572  0.01580127     0.86  0.04907138 0.223698
#> 2         b 50 -0.1949136 -0.01829474     0.70  0.06480741 0.798200
```
