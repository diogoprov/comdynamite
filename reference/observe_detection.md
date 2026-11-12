# Imperfect detection of a community

Turns the true abundances of a simulation into what a survey would
record. Each individual is detected in each of `visits` visits with
probability `detection`, independently, so the count of a species at a
site in a visit is binomial and a species with more individuals is more
likely to be found, as in the models of Royle & Nichols (2003) and Royle
(2004). The record keeps the largest count over visits, or only whether
the species was detected. The same applies to every step of `series`,
with new visits at each step. Use it inside `simulate` in
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
with `detection` drawn from a prior, so that the reference table carries
the same observation error as the data.

## Usage

``` r
observe_detection(x, detection, visits = 3, record = c("counts", "presence"))
```

## Arguments

- x:

  An `mc_sim` object or an
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  object with true abundances.

- detection:

  Probability that an individual is detected in a visit: one number, or
  one per species, named by species.

- visits:

  Number of visits to each site.

- record:

  `"counts"` for the largest count over visits or `"presence"` for
  detected at least once (1) or not (0).

## Value

An
[`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
object with the recorded `comm` and `series`, and the attribute
`detection`.

## References

Royle, J. A., & Nichols, J. D. (2003). Estimating abundance from
repeated presence–absence data or point counts. *Ecology*, 84, 777–790.
[doi:10.1890/0012-9658(2003)084\[0777:EAFRPA\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282003%29084%5B0777%3AEAFRPA%5D2.0.CO%3B2)

Royle, J. A. (2004). N-mixture models for estimating population size
from spatially replicated counts. *Biometrics*, 60, 108–115.
[doi:10.1111/j.0006-341X.2004.00142.x](https://doi.org/10.1111/j.0006-341X.2004.00142.x)

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
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
tree <- treeshapesim::sim_tree_shape(8)
sim <- simulate_metacommunity(landscape_grid(side = 4), species_traits(tree, 0.5),
                              competition_matrix(tree, "none"), dispersal = 0.05,
                              n_steps = 40)
seen <- observe_detection(sim, detection = 0.1, visits = 3)
c(true = sum(sim$abundance > 0), recorded = sum(seen$comm > 0))
#>     true recorded 
#>       91       66 
```
