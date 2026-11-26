# Plot the posterior of an ABC fit against the prior

For one observation, draws for each parameter the histogram of the
simulations of the reference table weighted by its random forest, an
approximate sample from the posterior, over the histogram of the
unweighted simulations, which are draws from the prior. When the fit has
the `prior` of
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
each panel spans the prior range. Each parameter has its own weights, so
the panels are marginal posteriors.

## Usage

``` r
# S3 method for class 'mc_abc_fit'
plot(x, newdata, bins = 30, colour = "#0072B2", ...)
```

## Arguments

- x:

  An `mc_abc_fit` from
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).

- newdata:

  Named numeric vector of observed statistics, as from
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md),
  or a data frame with one row.

- bins:

  Number of bins of the histograms. One whole number, 1 or more.

- colour:

  Colour of the posterior.

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md)
for the effective number of simulations behind the posterior.

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
[`plot.mc_abc_choice()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_choice.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
table <- data.frame(a = runif(300), b = runif(300))
table$s1 <- table$a + rnorm(300, 0, 0.05)
table$s2 <- table$b * table$a + rnorm(300, 0, 0.05)
fit <- abc_fit(table, parameters = c("a", "b"), ntree = 100)
plot(fit, c(s1 = 0.5, s2 = 0.2))
```
