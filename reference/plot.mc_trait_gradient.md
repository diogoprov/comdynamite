# Plot a trait-gradient analysis

`type = "gradient"` draws the trait value of each population against the
mean of its site, with the line of each species over the range of its
sites (its slope \\b_i\\) and the species means at \\(\beta_i, \bar
t_i)\\; the dashed line is \\y = x\\, where \\\alpha = 0\\.
`type = "beta_alpha"` draws \\\alpha\\ against \\\beta\\. For
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
the null distributions of the three statistics with the observed value.

## Usage

``` r
# S3 method for class 'mc_trait_gradient'
plot(
  x,
  type = c("gradient", "beta_alpha"),
  species_labels = FALSE,
  colour = "#0072B2",
  point_size = 1.5,
  ...
)

# S3 method for class 'mc_trait_gradient_test'
plot(x, ...)
```

## Arguments

- x:

  Result of
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
  or
  [`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md).

- type:

  `"gradient"` or `"beta_alpha"`.

- species_labels:

  Write the names of species in `type = "beta_alpha"`.

- colour:

  Colour of the species lines and means.

- point_size:

  Size of the points.

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md)

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
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
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance")
plot(tg)

plot(tg, type = "beta_alpha", species_labels = TRUE)

plot(trait_gradient_test(tg, n_null = 99))
```
