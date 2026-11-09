# Variance components of traits across nested levels

Splits the variance of each trait among nested levels, such as regions,
sites, species and populations, with a linear mixed model of random
intercepts nested from the broadest level to the finest, fitted by
restricted maximum likelihood (Messier et al. 2010).

## Usage

``` r
nested_variance(data, traits, factors)
```

## Arguments

- data:

  Data frame with the traits and the factors.

- traits:

  Names of the columns with trait values.

- factors:

  Names of the columns with the nested factors, from the broadest level
  to the finest.

## Value

A data frame with one row per factor and one for the variance within the
finest level (`within`), and one column per trait with the share of
variance of each level. The attribute `variance` holds the variances. It
has class `nested_variance`, with a
[plot()](https://provetelab.org/comdynamite/reference/plot.nested_variance.md)
method.

## Details

The levels must be nested: each level of a factor must lie within a
single level of the factor before it. Two checks guard against the
common mistakes. A factor that does not split any group of the factors
before it (for instance species after populations, since each population
holds a single species) adds no level of variation, and the function
stops, because the order of the factors is then most likely reversed. A
factor whose levels occur under several groups of the factors before it
(for instance genera, found at several sites, after sites) is crossed
with them rather than nested; the model then treats each combination as
a separate group, which is rarely what is meant, and the function warns.

## References

Messier, J., McGill, B. J., & Lechowicz, M. J. (2010). How do traits
vary across ecological scales? A case for trait-based ecology. *Ecology
Letters*, 13, 838–848.
[doi:10.1111/j.1461-0248.2010.01476.x](https://doi.org/10.1111/j.1461-0248.2010.01476.x)

## See also

[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md).

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
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
[`trait_gradient_brms()`](https://provetelab.org/comdynamite/reference/trait_gradient_brms.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
d <- data.frame(region = rep(c("north", "south"), each = 60),
                site = rep(paste0("s", 1:6), each = 20),
                species = sample(c("a", "b", "c"), 120, replace = TRUE))
d$population <- paste(d$site, d$species)
d$trait <- c(a = 10, b = 20, c = 30)[d$species] + rnorm(120)
res <- nested_variance(d, "trait", c("region", "site", "population"))
res
#>        level        trait
#> 1     region 1.741575e-08
#> 2       site 1.125857e-17
#> 3 population 9.862820e-01
#> 4     within 1.371800e-02
if (requireNamespace("ggplot2", quietly = TRUE)) plot(res)
```
