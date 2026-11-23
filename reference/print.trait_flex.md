# Model, fitted values and residuals of the decomposition of community trait means

Accessors and print method for the result of
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md).

## Usage

``` r
# S3 method for class 'trait_flex'
print(x, digits = 3, ...)

# S3 method for class 'trait_flex'
formula(x, ...)

# S3 method for class 'trait_flex'
fitted(object, ...)

# S3 method for class 'trait_flex'
residuals(object, ...)

# S3 method for class 'trait_flex'
nobs(object, ...)
```

## Arguments

- x, object:

  Result of
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md).

- digits:

  Number of significant digits printed. One whole number, 1 or more.

- ...:

  Not used.

## Value

[`formula()`](https://rdrr.io/r/stats/formula.html): the one-sided
formula of the environmental model.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html): the column of
the same name of `object$model`, one value per site, trait and component
(`turnover`, `intraspecific`, `total`), in the units of the trait;
without a model the fitted value is the mean over sites.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): the number of sites with
values of all variables of the model.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md).

Other trait variation among individuals:
[`coef.drift_test()`](https://provetelab.org/comdynamite/reference/coef.drift_test.md),
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
ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
                  species = sample(c("a", "b", "c"), 120, replace = TRUE))
env <- data.frame(site = paste0("s", 1:6), moisture = 1:6)
ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] +
  env$moisture[match(ind$site, env$site)] + rnorm(120)
res <- trait_flex(ind, formula = ~ moisture, env = env)
res
#> Decomposition of community trait means (Lepš et al. 2011): 6 sites, 1 trait
#> Model: ~moisture 
#>  trait     component explained residual total        p
#>  trait      turnover      5.63   2.4400  8.07 3.84e-02
#>  trait intraspecific     17.20   0.1480 17.30 2.73e-05
#>  trait   covariation     19.70   0.0833 19.80       NA
#>  trait         total     42.50   2.6700 45.20 1.34e-03
formula(res)
#> ~moisture
#> <environment: 0x5a9976c965c0>
head(residuals(res))
#> [1]  0.9149881 -1.2164700  0.1542713 -0.1912190  0.2105633  0.1278662
```
