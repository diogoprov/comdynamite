# Folds of sites for cross-validation

Splits the sites into `k` folds, either as spatial blocks, clusters of
neighbouring sites found by k-means on the coordinates, or at random.
With spatial autocorrelation, sites left out at random have close
neighbours among the sites kept, so random folds can make a model look
better at predicting new sites than it is; spatial blocks hold out whole
regions (Roberts et al. 2017).

## Usage

``` r
spatial_folds(coords, k = 5, method = c("blocks", "random"), seed = 1)
```

## Arguments

- coords:

  Matrix or data frame with one row per site and the coordinates in
  columns, such as `landscape$coords` or the `coords` of an
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md).

- k:

  Number of folds.

- method:

  `"blocks"` for spatial blocks or `"random"`.

- seed:

  Integer that sets the seed.

## Value

Integer vector with the fold of each site.

## References

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J.,
Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder,
B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann,
C. F. (2017). Cross-validation strategies for data with temporal,
spatial, hierarchical, or phylogenetic structure. *Ecography*, 40,
913–929. [doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881)

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
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md)

## Examples

``` r
grid <- landscape_grid(side = 6)
table(spatial_folds(grid$coords, k = 4))
#> 
#> 1 2 3 4 
#> 9 9 9 9 
```
