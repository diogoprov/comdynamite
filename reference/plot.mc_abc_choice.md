# Plot an ABC model choice

Without `newdata`, draws the out-of-bag confusion matrix of
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md)
as the share of the simulations of each true model (rows) assigned to
each model (columns), which shows before looking at the data which
models can be told apart. With `newdata`, draws the votes of the trees
for each model for one observation, with the chosen model in `colour`
and its posterior probability in the subtitle, which is more reliable
than the share of votes.

## Usage

``` r
# S3 method for class 'mc_abc_choice'
plot(x, newdata = NULL, ntree = 500, colour = "#0072B2", ...)
```

## Arguments

- x:

  An `mc_abc_choice` from
  [`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md).

- newdata:

  `NULL`, or a named numeric vector of observed statistics or a data
  frame with one row.

- ntree:

  Number of trees for the posterior probability, passed to
  [`predict()`](https://rdrr.io/r/stats/predict.html). One whole number,
  1 or more.

- colour:

  Colour of the chosen model, or of the largest shares.

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md)

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
[`nobs.mc_node_test()`](https://provetelab.org/comdynamite/reference/nobs.mc_node_test.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`plot.mc_abc_fit()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_fit.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
a <- data.frame(s1 = rnorm(200, 0), s2 = rnorm(200))
b <- data.frame(s1 = rnorm(200, 1.5), s2 = rnorm(200))
choice <- abc_choose(list(neutral = a, niche = b), ntree = 100)
plot(choice)

plot(choice, c(s1 = 1.2, s2 = 0), ntree = 100)
```
