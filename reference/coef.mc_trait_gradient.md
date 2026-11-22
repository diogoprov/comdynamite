# Coefficients and numbers of observations of the trait-gradient analysis

Accessors for the results of
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
and
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md).

## Usage

``` r
# S3 method for class 'mc_trait_gradient'
coef(object, ...)

# S3 method for class 'mc_trait_gradient'
nobs(object, ...)

# S3 method for class 'mc_trait_gradient_phylo'
coef(object, ...)

# S3 method for class 'mc_trait_gradient_phylo'
nobs(object, ...)
```

## Arguments

- object:

  Result of
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
  or
  [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md).

- ...:

  Not used.

## Value

For
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`coef()`](https://rdrr.io/r/stats/coef.html) gives a matrix with one
row per species and the columns `beta`, `alpha`, `intercept` and `slope`
(the last two `NA` for species with too few sites), and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of populations
(site by species). For
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`coef()`](https://rdrr.io/r/stats/coef.html) gives a matrix with one
row per response and model (named `response:model`) and the columns
`mu`, `sigma2` and `a`, and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of species in
the fits.

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md).

Other trait variation among individuals:
[`coef.drift_test()`](https://provetelab.org/comdynamite/reference/coef.drift_test.md),
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
tg <- trait_gradient(sample_individuals(sim, n = 3, intraspecific_sd = 1))
#> trait_gradient() implements the method of:
#>   Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to community assembly: partitioning of species trait values into within- and among-community components. Ecology Letters, 10, 135-145. https://doi.org/10.1111/j.1461-0248.2006.01006.x
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
head(coef(tg))
#>         beta       alpha intercept        slope
#> t1  63.55145   3.4999605  70.97109 -0.061677261
#> t10 61.81770  -9.8716077  51.76083  0.002996902
#> t11 62.54357  -0.5798133  60.08829  0.029986505
#> t12 63.50042   0.6812822        NA           NA
#> t13 64.62860 -38.7360921  18.94203  0.107544884
#> t15 59.64996 -60.0514282        NA           NA
nobs(tg)
#> [1] 113
```
