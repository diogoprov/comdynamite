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
  nperm = 99,
  alternative = c("two.sided", "less", "greater"),
  multivariate = FALSE
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- index:

  Function of the trait values of a site (its first argument: a vector,
  or a matrix with `multivariate = TRUE`) and optionally of their
  species (an argument named `species`), returning a single number, such
  as [`stats::var`](https://rdrr.io/r/stats/cor.html) or
  `function(x, species) stats::var(tapply(x, species, mean))`.

- traits:

  Names of the columns with trait values.

- site, species:

  Names of the columns with the site and the species.

- null:

  Null model (see Details).

- nperm:

  Number of randomisations; 0 for the observed values only.

- alternative:

  Direction of the test: `"two.sided"` (the default), `"less"` or
  `"greater"`.

- multivariate:

  If `TRUE`, the index receives all traits at once.

## Value

A data frame with one row per site and trait (or per site with
`multivariate = TRUE`, trait `"multivariate"`) and the columns `site`,
`trait`, `observed`, and with `nperm > 0` also `ses` and `p`; the
attribute `null` holds the null values (randomisation x row).
Individuals with a missing value are left out for that trait (for any
trait with `multivariate = TRUE`).

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

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:5), each = 20),
                  species = sample(letters[1:4], 100, replace = TRUE))
ind$trait <- c(a = 1, b = 2, c = 3, d = 4)[ind$species] + rnorm(100, sd = 0.3)
# range of trait values within sites, against draws from the whole region
community_index(ind, function(x) diff(range(x)), null = "regional.ind", nperm = 99)
#>   site trait observed         ses    p
#> 1   s1 trait 3.155201 -1.79804612 0.12
#> 2   s2 trait 3.695020  0.53805219 0.68
#> 3   s3 trait 3.843434  1.13858880 0.40
#> 4   s4 trait 3.297534 -0.97698903 0.42
#> 5   s5 trait 3.579527 -0.02346576 0.94
# variance among the population means of a site
pop_var <- function(x, species) stats::var(tapply(x, species, mean))
community_index(ind, pop_var, null = "local", nperm = 99)
#>   site trait observed      ses    p
#> 1   s1 trait 1.476106 4.219795 0.02
#> 2   s2 trait 1.577807 7.631109 0.02
#> 3   s3 trait 2.090731 7.206326 0.02
#> 4   s4 trait 1.661834 8.237044 0.02
#> 5   s5 trait 1.713591 4.415371 0.02
```
