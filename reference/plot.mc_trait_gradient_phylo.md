# Plot the phylogenetic models of a trait-gradient analysis

Draws, for each response of
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
the Akaike weight of each model of trait evolution, with the best model
in `colour` and the AIC difference from it above each bar. Check
`at_bound` and `converged` in `x$models` before reading an OU fit.

## Usage

``` r
# S3 method for class 'mc_trait_gradient_phylo'
plot(x, colour = "#0072B2", ...)
```

## Arguments

- x:

  Result of
  [`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md).

- colour:

  Colour of the best model of each response.

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)

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
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance")
plot(trait_gradient_phylo(tg, tree, n_boot = 99))
#> trait_gradient_phylo() implements the method of:
#>   Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to community assembly: partitioning of species trait values into within- and among-community components. Ecology Letters, 10, 135-145. https://doi.org/10.1111/j.1461-0248.2006.01006.x
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
```
