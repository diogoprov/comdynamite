# Trait indices of communities against null models

Computes any index of the trait values of the individuals of each site
(community), such as the variance, the range or a nearest-neighbour
distance, and compares it with randomised communities.

## Usage

``` r
community_index(
  individuals,
  index,
  traits = "trait",
  site = "site",
  species = "species",
  null = c("local", "regional.ind", "regional.pop", "regional.pop.prab"),
  n_perm = 99,
  alternative = c("two.sided", "less", "greater"),
  multivariate = FALSE,
  nperm = NULL
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- index:

  Function of the trait values of a site (its first argument: a vector,
  or a matrix with `multivariate = TRUE`) and optionally of their
  species (an argument named `species`), returning a single number, such
  as [`stats::var`](https://rdrr.io/r/stats/cor.html),
  [`mnnd()`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md)
  or `function(x, species) stats::var(tapply(x, species, mean))`; or a
  named list of such functions, all computed on the same observed and
  randomised communities, such as
  `list(variance = stats::var, MNND = mnnd, range = trait_range)`.

- traits:

  Names of the columns with trait values.

- site, species:

  Names of the columns with the site and the species.

- null:

  Null model (see Details).

- n_perm:

  Number of randomisations; 0 for the observed values only.

- alternative:

  Direction of the test: `"two.sided"` (the default), `"less"` or
  `"greater"`.

- multivariate:

  If `TRUE`, the index receives all traits at once.

- nperm:

  Deprecated; use `n_perm`.

## Value

A data frame of class `community_index` with one row per site, trait (or
`"multivariate"` with `multivariate = TRUE`) and index, and the columns
`site`, `trait`, `index` (the name in the list, or the name of the
function), `observed`, and with `n_perm > 0` also `ses` and `p`; the
attribute `null` holds the null values (randomisation x row), and
[plot()](https://provetelab.org/comdynamite/reference/plot.tstats.md)
draws the tests. Individuals with a missing value are left out for that
trait (for any trait with `multivariate = TRUE`).

## Details

The null models follow Taudiere & Violle (2016), written from their
description:

- `"local"`: trait values are permuted among the individuals of each
  site; the values of the site are kept and only their link with species
  is broken, so it is useful for indices that use `species`;

- `"regional.ind"`: each site receives as many trait values as it has
  individuals, drawn without replacement from all individuals of the
  data;

- `"regional.pop"`: each individual takes the mean of its population,
  and each site receives as many of these values as it has individuals,
  drawn without replacement from all individuals. The observed index
  uses the values of the individuals, so the null removes the variation
  within populations as well as the species composition of the site;

- `"regional.pop.prab"`: the index is computed on the means of the
  populations of each site (one value per population, observed and null
  alike), and each site receives as many population means as it has
  populations, drawn without replacement from all populations.

In all cases the individuals (or populations) keep their species, so an
index with an argument named `species`, such as `function(x, species)`,
receives the species of the values of the site. With
`multivariate = TRUE` the index receives the matrix of all traits of the
site, and randomisations move whole rows, so that the traits of an
individual stay together. Sites and species are matched by name, so the
order of the rows does not matter. An index that the null model cannot
change (for example [`stats::var`](https://rdrr.io/r/stats/cor.html)
under `"local"`) gets `ses = NA`.

## References

Taudiere, A., & Violle, C. (2016). cati: an R package using functional
traits to detect and quantify multi-level community assembly processes.
*Ecography*, 39, 699–708.
[doi:10.1111/ecog.01433](https://doi.org/10.1111/ecog.01433)

## See also

[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md).

Other trait variation among individuals:
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
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:5), each = 20),
                  species = sample(letters[1:4], 100, replace = TRUE))
ind$trait <- c(a = 1, b = 2, c = 3, d = 4)[ind$species] + rnorm(100, sd = 0.3)
# range of trait values within sites, against draws from the whole region
community_index(ind, function(x) diff(range(x)), null = "regional.ind", n_perm = 99)
#>   site trait index observed         ses    p
#> 1   s1 trait index 3.155201 -1.79804612 0.12
#> 2   s2 trait index 3.695020  0.53805219 0.68
#> 3   s3 trait index 3.843434  1.13858880 0.40
#> 4   s4 trait index 3.297534 -0.97698903 0.42
#> 5   s5 trait index 3.579527 -0.02346576 0.94
# variance among the population means of a site
pop_var <- function(x, species) stats::var(tapply(x, species, mean))
community_index(ind, pop_var, null = "local", n_perm = 99)
#>   site trait   index observed      ses    p
#> 1   s1 trait pop_var 1.476106 4.219795 0.02
#> 2   s2 trait pop_var 1.577807 7.631109 0.02
#> 3   s3 trait pop_var 2.090731 7.206326 0.02
#> 4   s4 trait pop_var 1.661834 8.237044 0.02
#> 5   s5 trait pop_var 1.713591 4.415371 0.02
# several indices at once, on the same randomised communities
community_index(ind, list(variance = stats::var, MNND = mnnd, kurtosis = trait_kurtosis),
                null = "regional.ind", n_perm = 99)
#>    site trait    index    observed          ses    p
#> 1    s1 trait variance  0.96184073 -1.204549411 0.24
#> 2    s2 trait variance  1.01906139 -1.315648430 0.24
#> 3    s3 trait variance  1.55696893  0.955877885 0.34
#> 4    s4 trait variance  0.88096775 -1.394588992 0.16
#> 5    s5 trait variance  1.72488384  2.015358652 0.06
#> 6    s1 trait     MNND  0.08888870 -0.721103852 0.50
#> 7    s2 trait     MNND  0.11849284  0.847592936 0.46
#> 8    s3 trait     MNND  0.12324692  1.028942015 0.38
#> 9    s4 trait     MNND  0.04634544 -2.976977563 0.02
#> 10   s5 trait     MNND  0.09741112 -0.370322230 0.72
#> 11   s1 trait kurtosis -1.05229048 -0.007918248 0.90
#> 12   s2 trait kurtosis -0.10003062  3.187360206 0.06
#> 13   s3 trait kurtosis -1.23233713 -0.420626901 0.76
#> 14   s4 trait kurtosis -0.68087734  0.906284555 0.36
#> 15   s5 trait kurtosis -1.55901734 -1.373313414 0.02
```
