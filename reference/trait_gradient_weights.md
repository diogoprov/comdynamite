# Trait-gradient analysis with the two choices of weights

Runs
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
with each combination of weights of the site means and of the species
means (abundance or occurrence) and puts the values of \\\beta\\ and
\\\alpha\\ of each species side by side, with the correlations among the
combinations. Large differences mean that the result depends on the
dominant species.

## Usage

``` r
trait_gradient_weights(data, abundance, ...)
```

## Arguments

- data:

  Data frame with one row per individual or per population (species in a
  site), such as from
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).
  Rows of one population are averaged.

- abundance:

  Name of the column with abundances (needed here).

- ...:

  Further arguments to
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md).

## Value

A list with `beta` and `alpha`, data frames with one row per species and
one column per combination (`abundance/abundance`,
`abundance/occurrence`, `occurrence/abundance`, `occurrence/occurrence`,
sites first), and `correlation`, the correlations of \\\beta\\ and of
\\\alpha\\ between each combination and the first.

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
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
ind <- sample_individuals(sim, n = 5, intraspecific_sd = 1)
trait_gradient_weights(ind, abundance = "abundance")$correlation
#>                 weights      beta     alpha
#> 1   abundance/abundance 1.0000000 1.0000000
#> 2  abundance/occurrence 0.7747809 0.9483062
#> 3  occurrence/abundance 0.7130104 0.9773422
#> 4 occurrence/occurrence 0.3973584 0.9702739
```
