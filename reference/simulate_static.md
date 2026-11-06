# Simulate a static metacommunity

Draws species abundances directly from unimodal responses to an
environmental gradient, without population dynamics, competition or
dispersal. The expected abundance of species j at site i is
`r_max[j] * exp(-(env[i] - optimum[j])^2 / (2 * sigma[j]^2))`, the
Gaussian response of Dray & Legendre (2008), also used by Peres-Neto et
al. (2012). By default, observed abundances are Poisson counts with that
mean. This is our addition: in those simulations the abundances are the
expected values, with `r_max` between 0.5 and 1, which `poisson = FALSE`
returns. Because nothing but the gradient structures the data, this is a
baseline for
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
and a fast generator for evaluating methods.

## Usage

``` r
simulate_static(traits, env, incidence = FALSE, poisson = TRUE)
```

## Arguments

- traits:

  Data frame from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  or any data frame with columns `optimum`, `sigma` and `r_max`, or a
  list of them, one per gradient, with the same species. Here `r_max` is
  the maximum expected abundance, not a growth rate, so set `rmax_range`
  accordingly. With several gradients the maxima multiply, so a range
  such as `c(1, 1)` in all but one of them keeps abundances on the scale
  of the first.

- env:

  Value of the gradient at each site, or a list or data frame with one
  gradient per element, in the order of `traits`.

- incidence:

  Return presences and absences instead of abundances? A species is
  present where its Poisson count is positive, so this needs
  `poisson = TRUE`. Dray & Legendre (2008) instead draw presences from a
  binomial distribution with the expected value as the probability.

- poisson:

  Draw Poisson counts? If `FALSE`, return the expected abundances.

## Value

List of class `mc_sim` with `abundance`, a matrix with one row per site
and one column per species, plus `env` and `traits`. Sites and species
with no individuals are kept, so remove them if a method requires it.

## Details

With several gradients, for example an environmental and a spatial one,
each species has its own response to each gradient and the expected
abundance is the product of the responses, as in Peres-Neto et al.
(2012).

## References

Dray, S., & Legendre, P. (2008). Testing the species traits-environment
relationships: the fourth-corner problem revisited. *Ecology*, 89,
3400–3412. [doi:10.1890/08-0349.1](https://doi.org/10.1890/08-0349.1)

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

## Examples

``` r
set.seed(1)
tr <- ape::rcoal(20)
traits <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
sim <- simulate_static(traits, env = runif(100, 0, 100))
sim
#> <mc_sim> 100 sites x 20 species
#> species persisting: 19 of 20 
#> mean local richness: 9 

# an environmental and a spatial gradient
space <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
                        rmax_range = c(1, 1))
gradients <- data.frame(env = runif(100, 0, 100), space = 1:100)
simulate_static(list(traits, space), gradients)
#> <mc_sim> 100 sites x 20 species
#> species persisting: 19 of 20 
#> mean local richness: 3.58 
```
