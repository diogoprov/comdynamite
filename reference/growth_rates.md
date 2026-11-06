# Growth rates of each species at each site

Density-independent growth rate of species j at site i,
`r_max[j] * exp(-(E[i] - optimum[j])^2 / (2 * sigma[j]^2))`: a Gaussian
response whose standard deviation is the niche breadth `sigma`. This
`sigma` equals `sqrt(2)` times the niche breadth parameter of Thompson
et al. (2020).

## Usage

``` r
growth_rates(env, traits)
```

## Arguments

- env:

  Environment of each site.

- traits:

  Data frame from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  or any data frame with columns `optimum`, `sigma` and `r_max`.

## Value

Numeric matrix with one row per site and one column per species: the
density-independent growth rate of each species at each site.

## Examples

``` r
set.seed(1)
g <- growth_rates(seq(0, 100, 10), species_traits(3))
round(g, 2)
#>        sp1  sp2  sp3
#>  [1,] 3.13 0.00 7.61
#>  [2,] 3.94 0.00 0.00
#>  [3,] 4.28 0.00 0.00
#>  [4,] 4.00 0.00 0.00
#>  [5,] 3.23 0.00 0.00
#>  [6,] 2.25 0.01 0.00
#>  [7,] 1.35 0.07 0.00
#>  [8,] 0.70 0.52 0.00
#>  [9,] 0.31 2.15 0.00
#> [10,] 0.12 5.02 0.00
#> [11,] 0.04 6.66 0.00
```
