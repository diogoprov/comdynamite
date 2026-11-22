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
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`plot.mc_abc_choice()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_choice.md),
[`plot.mc_abc_fit()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_fit.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)
