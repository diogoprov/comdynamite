# Numbers of observations of the node test and of ABC

Numbers of observations of the node test and of ABC

## Usage

``` r
# S3 method for class 'mc_node_test'
nobs(object, ...)

# S3 method for class 'mc_abc_fit'
nobs(object, ...)
```

## Arguments

- object:

  Result of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  or
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).

- ...:

  Not used.

## Value

For
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
the number of sites at the root, which holds all the sites analysed
(each node uses the sites where its clade occurs, given in the column
`n_sites` of `nodes`). For
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
the number of complete simulations in the training table.

## See also

[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).

Other inference:
[`abc_calibrate()`](https://provetelab.org/comdynamite/reference/abc_calibrate.md),
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_sequential()`](https://provetelab.org/comdynamite/reference/abc_sequential.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`plot.mc_abc_choice()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_choice.md),
[`plot.mc_abc_fit()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_fit.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
out <- node_test(daphniids$comm, node_contrasts(daphniids$tree), daphniids$env[, 1],
                 n_perm = 9)
#> node_test() implements the method of:
#>   Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity phylogenetics: separating the roles of environmental filters and historical biogeography. Ecology Letters, 13, 1290-1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
#>   Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. Ecology, 93, S14-S30. https://doi.org/10.1890/11-0494.1
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
nobs(out)
#> [1] 239
```
