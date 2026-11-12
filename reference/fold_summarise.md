# Statistics by fold of sites

Turns a function that summarises a data set into one that summarises,
for every fold, the sites kept (`train`) and the sites left out
(`test`). Use it as `summarise` in
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
so that one simulation gives the statistics of every fold, and apply it
to the observed data before
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md).
The statistics are named `f1_train_alpha`, `f1_test_alpha` and so on.

## Usage

``` r
fold_summarise(folds, summarise = function(x, shared) summary_stats(x))
```

## Arguments

- folds:

  Integer vector with the fold of each site, as from
  [`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md).

- summarise:

  Function of an
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  object and `shared` that gives a named numeric vector. The default is
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md).

## Value

A function of `sim`, a simulation or an
[`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
object, and `shared`.

## See also

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
tree <- treeshapesim::sim_tree_shape(8)
grid <- landscape_grid(side = 4)
sim <- simulate_metacommunity(grid, species_traits(tree, 0.5), competition_matrix(tree, "none"),
                              dispersal = 0.05, n_steps = 40)
by_fold <- fold_summarise(spatial_folds(grid$coords, k = 2),
                          function(x, shared) summary_stats(x, n_null = 9)[1:3])
by_fold(sim, list())
#> f1_train_alpha f1_train_gamma  f1_train_beta  f1_test_alpha  f1_test_gamma 
#>       5.000000       7.000000       1.400000       3.375000       6.000000 
#>   f1_test_beta f2_train_alpha f2_train_gamma  f2_train_beta  f2_test_alpha 
#>       1.777778       3.375000       6.000000       1.777778       5.000000 
#>  f2_test_gamma   f2_test_beta 
#>       7.000000       1.400000 
```
