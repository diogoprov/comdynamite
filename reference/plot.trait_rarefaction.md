# Plot a rarefaction of a statistic

Draws the mean of each statistic from
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md)
against sampling effort, with the interval over subsamples as a band,
one panel per statistic. Dashed vertical lines mark reference levels of
effort: for a rarefaction over `cover`, by default the published rule of
measuring at least the species that make up 80% of the vegetation
(Pakeman & Quested 2007). No such rule is drawn by default for the
number of individuals per population; give it in `reference` when the
field or study has one.

## Usage

``` r
# S3 method for class 'trait_rarefaction'
plot(
  x,
  reference = NULL,
  colours = "#0072B2",
  alpha = 0.25,
  scales = "free_y",
  ...,
  colour = NULL
)
```

## Arguments

- x:

  Result of
  [`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md).

- reference:

  Levels of effort marked by dashed lines; named values are written next
  to the lines. `NULL` (the default) marks the 80% rule for `cover` and
  nothing for `n`; `FALSE` draws no line.

- colours:

  Colour of the line and band; the default is the blue of the Okabe-Ito
  palette.

- alpha:

  Transparency of the band.

- scales:

  Passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

- ...:

  Not used.

- colour:

  Deprecated; use `colours`, as in the other plots.

## Value

A ggplot object.

## References

Pakeman, R. J., & Quested, H. M. (2007). Sampling plant functional
traits: what proportion of the species need to be measured? *Applied
Vegetation Science*, 10, 91–96.
[doi:10.1111/j.1654-109x.2007.tb00507.x](https://doi.org/10.1111/j.1654-109x.2007.tb00507.x)

## See also

[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md).

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
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
ind <- data.frame(site = rep(paste0("s", 1:5), each = 40),
                  species = sample(c("a", "b", "c", "d"), 200, TRUE, prob = c(5, 3, 1, 1)))
ind$trait <- c(a = 10, b = 12, c = 15, d = 18)[ind$species] + rnorm(200)
plot(trait_rarefaction(ind, n = 2:10, reps = 10))

plot(trait_rarefaction(ind, cover = seq(0.3, 1, by = 0.1)))

# a rule of the study, for the number of individuals
plot(trait_rarefaction(ind, n = 2:10, reps = 10), reference = c("protocol" = 5))
```
