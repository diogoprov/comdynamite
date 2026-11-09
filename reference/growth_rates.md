# Growth rates of each species at each site

Density-independent growth rate of species \\j\\ at site \\i\\, \\r_j
\exp(-(E_i - \theta_j)^2 / (2 \sigma_j^2))\\, with \\r_j\\ = `r_max`,
\\\theta_j\\ = `optimum` and \\\sigma_j\\ = `sigma`: a Gaussian response
whose standard deviation is the niche breadth `sigma`. This `sigma`
equals \\\sqrt{2}\\ times the niche breadth parameter of Thompson et al.
(2020).

## Usage

``` r
growth_rates(env, traits)
```

## Arguments

- env:

  Environment of each site; for several axes, a matrix, data frame or
  list with one axis per column or element, in the order of `traits`.

- traits:

  Data frame from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  or any data frame with columns `optimum`, `sigma` and `r_max`; for
  several axes, a list of them with the same species.

## Value

Numeric matrix with one row per site and one column per species: the
density-independent growth rate of each species at each site.

## Details

With several environmental axes, each species has a response to each
axis and the growth rate is the product of the responses, as in
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md):
a Gaussian niche in several dimensions whose axes are independent (no
correlation between the responses to different axes), with the optimum
and breadth of each axis in its own set of traits. The maxima multiply
too, so give the traits of all but one axis `r_max` of 1, for example
with `species_traits(rmax_range = c(1, 1))`.

## See also

Other simulate metacommunities:
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)

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
