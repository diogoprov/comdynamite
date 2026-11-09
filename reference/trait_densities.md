# Trait probability densities of populations and communities

Estimates the trait probability density (TPD; Carmona et al. 2016, 2019)
of each population with a Gaussian kernel, and that of each community as
the sum of the densities of its populations weighted by their relative
abundances, so that both integrate to one.

## Usage

``` r
trait_densities(
  individuals,
  trait = "trait",
  site = "site",
  species = "species",
  abundance = NULL,
  bw = NULL,
  n = 512
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- trait:

  Name of the column with the trait.

- site, species:

  Names of the columns with the site and the species.

- abundance:

  Name of a column with the abundance of each population (the same value
  for all its individuals), or `NULL` to weight populations by the
  number of individuals measured.

- bw:

  Bandwidth of the kernel, or `NULL` for the rule above.

- n:

  Number of points of the grid of trait values.

## Value

A data frame with the columns `site`, `species` (`NA` for communities),
`level` (`"population"` or `"community"`), `weight` (relative abundance
of the population in the site, 1 for communities), `trait` (grid value)
and `density`. The bandwidth is in the attribute `bw`.

## Details

All populations share one bandwidth, so that their densities can be
compared and summed. By default it follows Silverman's rule of thumb (as
[`stats::bw.nrd0()`](https://rdrr.io/r/stats/bandwidth.html)) applied to
the deviations of individuals from the mean of their population, pooled
over all populations with two or more individuals (the standard
deviation with \\n_i - 1\\ degrees of freedom per population, and the
interquartile range of the deviations rescaled by \\\sqrt{n_i / (n_i -
1)}\\, so that small populations do not shrink the bandwidth), with the
mean number of individuals per population as the sample size. A
population of one individual then gets a kernel of that width around its
value, instead of no density.

One bandwidth suits populations whose spread is about the same. For
traits whose spread grows with their mean, such as body size across
species of very different sizes, the bandwidth follows the many small
populations and the densities of large ones break into one peak per
individual; the logarithm of the trait usually evens out the spread.
With the default bandwidth, a warning is given when the trait is
positive and the standard deviation of populations of three or more
individuals rises with their mean (Spearman correlation above 0.7 over
at least ten populations).

## References

Carmona, C. P., de Bello, F., Mason, N. W. H., & Lepš, J. (2016). Traits
without borders: integrating functional diversity across scales. *Trends
in Ecology & Evolution*, 31, 382–394.
[doi:10.1016/j.tree.2016.02.003](https://doi.org/10.1016/j.tree.2016.02.003)

Carmona, C. P., de Bello, F., Mason, N. W. H., & Lepš, J. (2019). Trait
probability density (TPD): measuring functional diversity across scales
based on TPD with R. *Ecology*, 100, e02876.
[doi:10.1002/ecy.2876](https://doi.org/10.1002/ecy.2876)

## See also

[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md).

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(c("s1", "s2"), each = 30),
                  species = sample(c("a", "b"), 60, replace = TRUE))
ind$size <- c(a = 10, b = 14)[ind$species] + rnorm(60)
d <- trait_densities(ind, "size")
# each density integrates to one
tapply(d$density, paste(d$site, d$species), sum) * diff(d$trait[1:2])
#>     s1 NA      s1 a      s1 b     s2 NA      s2 a      s2 b 
#> 0.9999378 0.9999067 1.0000000 0.9999521 0.9999966 0.9999263 
```
