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
  seed = 1,
  min_node_size = 5
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

  Number of pseudo-observations. One whole number, 2 or more and at
  least 11 fewer than the rows of `table`.

- level:

  Probability of the credible interval. One number greater than 0 and
  less than 1.

- ntree:

  Number of trees per forest. One whole number, 1 or more.

- ncores:

  Number of threads for the forests. One whole number, 1 or more.

- seed:

  Integer for the choice of the pseudo-observations. One finite whole
  number.

- min_node_size:

  Smallest number of simulations in a node for it to be split, the
  `min.node.size` of
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).
  One whole number, 1 or more. The default, 5, is that of Raynal et al.
  (2019); larger values give larger leaves, which widen intervals that
  `abc_recovery()` shows to be too narrow, at the cost of more
  smoothing.

## Value

A data frame with one row per parameter: `n`, `r2` (1 minus the squared
error over the variance of the truth), `bias` and `width` (mean width of
the interval) as fractions of the prior range, and `coverage`, the
fraction of intervals that contain the truth, with its standard error.
The attribute `estimates` holds every estimate. A coverage near `level`
says the intervals are calibrated, one above it that they are
conservative and one below it that they are too narrow (see Details); R2
says how much the data narrow the prior.

## Details

The credible intervals of
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
are approximate, so their coverage is something to measure and not to
assume, and `coverage` is that measure for this table: the fraction of
pseudo-observations, drawn from the prior, whose interval holds the
truth. Expect it at or above `level` with several statistics, since the
forests smooth the posterior (Raynal et al. 2019 report slightly
conservative intervals; 0.92 to 0.99 at a `level` of 0.9 in our extended
tests), and equal to `level` for a parameter the statistics say nothing
about, whose posterior is the prior. A coverage below `level` by more
than about two `coverage_se` says the intervals are too narrow, which
happens with one or two statistics and small leaves (about 0.65 with one
statistic and the default `min_node_size`); a larger `min_node_size`
widens them. The coverage is an average over the prior and does not say
that the interval of one observation is right.

## References

Pontarp, M., Brännström, Å., & Petchey, O. L. (2019). Inferring
community assembly processes from macroscopic patterns using dynamic
eco-evolutionary models and Approximate Bayesian Computation (ABC).
*Methods in Ecology and Evolution*, 10, 450–460.
[doi:10.1111/2041-210X.13129](https://doi.org/10.1111/2041-210X.13129)

Raynal, L., Marin, J.-M., Pudlo, P., Ribatet, M., Robert, C. P., &
Estoup, A. (2019). ABC random forests for Bayesian parameter inference.
*Bioinformatics*, 35, 1720–1728.
[doi:10.1093/bioinformatics/bty867](https://doi.org/10.1093/bioinformatics/bty867)

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
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
