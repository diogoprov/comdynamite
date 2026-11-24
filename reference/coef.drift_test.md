# Coefficients, intervals and fitted values of the test of drift

Accessors for the result of
[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
which regresses the logarithm of the variance among group means on the
logarithm of the eigenvalues of `W`, one point per principal component.

## Usage

``` r
# S3 method for class 'drift_test'
coef(object, ...)

# S3 method for class 'drift_test'
confint(object, parm = NULL, level = 0.95, ...)

# S3 method for class 'drift_test'
vcov(object, ...)

# S3 method for class 'drift_test'
nobs(object, ...)

# S3 method for class 'drift_test'
fitted(object, ...)

# S3 method for class 'drift_test'
residuals(object, ...)
```

## Arguments

- object:

  Result of
  [`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md).

- ...:

  Not used.

- parm:

  Names of the coefficients (`"intercept"`, `"slope"`); all by default.

- level:

  Confidence level, one number between 0 and 1, both excluded.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named vector with
`intercept` and `slope`.
[`confint()`](https://rdrr.io/r/stats/confint.html): a matrix with one
row per coefficient and the lower and upper limits, from the \\t\\
distribution with `n_traits - 2` degrees of freedom.
[`vcov()`](https://rdrr.io/r/stats/vcov.html): their covariance matrix.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): the number of points of
the regression, the number of traits.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html): one value per
principal component, on the scale of `log(between)`.

## See also

[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md).

Other trait variation among individuals:
[`coef.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/coef.mc_trait_gradient.md),
[`coef.population_shifts()`](https://provetelab.org/comdynamite/reference/coef.population_shifts.md),
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.mc_trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient_bayes.md),
[`plot.mc_trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient_phylo.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
[`population_stats()`](https://provetelab.org/comdynamite/reference/population_stats.md),
[`print.trait_flex()`](https://provetelab.org/comdynamite/reference/print.trait_flex.md),
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(300 * 4), 300, 4) %*% diag(c(3, 2, 1.5, 1))
group <- rep(1:30, each = 10)
x <- x + (matrix(rnorm(30 * 4), 30, 4) %*% diag(c(3, 2, 1.5, 1)))[group, ]
res <- drift_test(x, group)
#> drift_test() implements the method of:
#>   Ackermann, R. R., & Cheverud, J. M. (2004). Detecting genetic drift versus selection in human evolution. Proceedings of the National Academy of Sciences, 101, 17946-17951. https://doi.org/10.1073/pnas.0405919102
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
coef(res)
#>  intercept      slope 
#> 0.06371541 1.01559120 
confint(res)
#>                2.5 %    97.5 %
#> intercept -0.6257980 0.7532288
#> slope      0.5247076 1.5064748
nobs(res)
#> [1] 4
```
