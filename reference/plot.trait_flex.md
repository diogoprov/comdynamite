# Plot the decomposition of community trait means

Draws, for each trait, the sums of squares among sites of the community
mean due to species turnover, intraspecific variability, their
covariation and the total, from
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md).
With an environmental model each bar is split into the part explained by
the model and the residual part. A negative covariation (turnover and
intraspecific variability moving the community mean in opposite
directions) is drawn below zero. Needs the package ggplot2.

## Usage

``` r
# S3 method for class 'trait_flex'
plot(x, relative = TRUE, colours = c("#0072B2", "grey80"), ...)
```

## Arguments

- x:

  Result of
  [`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md).

- relative:

  If `TRUE` (the default), sums of squares are divided by the total sum
  of squares of each trait, so that traits on different scales can be
  compared.

- colours:

  Two colours, for the explained and the residual parts (only the first
  is used without a model). The default is the blue of the Okabe-Ito
  palette and a light grey, as in
  [`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md).

- ...:

  Not used.

## Value

A ggplot object, which can be changed further with ggplot2 functions.

## See also

[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md).

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
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
ind <- data.frame(site = rep(paste0("s", 1:8), each = 20),
                  species = sample(c("a", "b", "c"), 160, replace = TRUE))
env <- data.frame(site = paste0("s", 1:8), moisture = 1:8)
shift <- env$moisture[match(ind$site, env$site)]
ind$height <- c(a = 10, b = 20, c = 30)[ind$species] + shift + rnorm(160)
ind$leaf <- c(a = 3, b = 2, c = 1)[ind$species] - 0.3 * shift + rnorm(160)
plot(trait_flex(ind, traits = c("height", "leaf"), formula = ~ moisture, env = env))
```
