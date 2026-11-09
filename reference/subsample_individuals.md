# Subsamples of individuals

Draws a subsample of the individuals of each population of real data,
with the sampling schemes of
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md):
a fixed number of individuals per population (`n`), a share of the
individuals (`fraction`), only the most abundant species of each site
(`cover`), and a minimum number of individuals per population (`min_n`).
It shows how much a result depends on sampling effort, alone or through
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md).

## Usage

``` r
subsample_individuals(
  individuals,
  n = NULL,
  fraction = NULL,
  cover = NULL,
  min_n = 1,
  site = "site",
  species = "species",
  abundance = NULL
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- n:

  Number of individuals drawn from each population, or `NULL` to keep
  all.

- fraction:

  Probability that each individual is kept; used when `n` is `NULL`.

- cover:

  Share (0 to 1) of the abundance of each site covered by the species
  kept, from the most abundant down (Pakeman & Quested 2007), or `NULL`
  to keep every species. Abundances are those of the `abundance` column,
  or the numbers of individuals in the data. Species tied in abundance
  are taken in the order of the data, and populations with an abundance
  of zero are never kept.

- min_n:

  Minimum number of individuals of a population after sampling.

- site, species:

  Names of the columns with the site and the species.

- abundance:

  Name of a column with the abundance of each population (the same value
  for all its individuals), or `NULL` to use the numbers of individuals
  in the data.

## Value

The rows of `individuals` kept, in their original order.

## Details

The schemes are applied in this order: species outside `cover` are
removed, individuals are drawn without replacement (`n` from each
population, all of them when the population is smaller, or each one with
probability `fraction`), and populations left with fewer than `min_n`
individuals are removed.

## References

Pakeman, R. J., & Quested, H. M. (2007). Sampling plant functional
traits: what proportion of the species need to be measured? *Applied
Vegetation Science*, 10, 91–96.
[doi:10.1111/j.1654-109x.2007.tb00507.x](https://doi.org/10.1111/j.1654-109x.2007.tb00507.x)

## See also

[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).

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
ind <- data.frame(site = rep(c("s1", "s2"), each = 50),
                  species = sample(c("a", "b", "c"), 100, TRUE, prob = c(6, 3, 1)),
                  trait = rnorm(100))
table(subsample_individuals(ind, n = 5)[c("site", "species")])
#>     species
#> site a b c
#>   s1 5 5 4
#>   s2 5 5 2
# only the species making up 80% of the individuals of each site
table(subsample_individuals(ind, cover = 0.8)[c("site", "species")])
#>     species
#> site  a  b
#>   s1 27 19
#>   s2 30 18
```
