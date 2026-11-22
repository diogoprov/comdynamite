# Coefficients, intervals and fitted values of population shifts

Accessors for the result of
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md).

## Usage

``` r
# S3 method for class 'population_shifts'
coef(object, ...)

# S3 method for class 'population_shifts'
confint(object, parm = NULL, level = attr(object, "level"), ...)

# S3 method for class 'population_shifts'
vcov(object, ...)

# S3 method for class 'population_shifts'
nobs(object, ...)

# S3 method for class 'population_shifts'
fitted(object, ...)

# S3 method for class 'population_shifts'
residuals(object, ...)
```

## Arguments

- object:

  Result of
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md).

- ...:

  Not used.

- parm:

  Names of the coefficients (the slopes `local_adaptation...` and
  `displacement`); all by default.

- level:

  Confidence level. The intervals are those computed by
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  with its `level`; another level stops with an error, since the
  bootstrap must be run again.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named vector with the
slopes on each environmental variable and the strength of displacement
`c`. [`confint()`](https://rdrr.io/r/stats/confint.html): a matrix with
one row per coefficient and the lower and upper limits of the basic
bootstrap intervals (needs `nboot > 0`).
[`vcov()`](https://rdrr.io/r/stats/vcov.html): the covariance of the
coefficients among the bootstrap simulations from the fitted model
(needs `nboot > 0`). [`nobs()`](https://rdrr.io/r/stats/nobs.html): the
number of populations that entered as responses.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html): one value per
half-sample mean of a population, centred within species, in the order
of the attribute `model` of `object`, at the chosen width.

## See also

[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md).

Other trait variation among individuals:
[`coef.drift_test()`](https://provetelab.org/comdynamite/reference/coef.drift_test.md),
[`coef.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/coef.mc_trait_gradient.md),
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
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = 2)
fit <- population_shifts(ind, width = 0.5, nboot = 19)
#> population_shifts() implements the method of:
#>   Goldberg, E. E., & Lande, R. (2006). Ecological and reproductive character displacement on an environmental gradient. Evolution, 60, 1344-1357. https://doi.org/10.1554/05-696.1
#>   Adams, D. C., & Collyer, M. L. (2007). Analysis of character divergence along environmental gradients and other covariates. Evolution, 61, 510-515. https://doi.org/10.1111/j.1558-5646.2007.00063.x
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
coef(fit)
#> local_adaptation     displacement 
#>        0.3165995      -16.5543448 
confint(fit)
#>                        2.5 %     97.5 %
#> local_adaptation   0.2949466  0.3257822
#> displacement     -41.8578167 17.2543841
nobs(fit)
#> [1] 91
```
