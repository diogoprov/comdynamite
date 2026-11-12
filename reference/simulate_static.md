# Simulate a static metacommunity

Draws species abundances directly from unimodal responses to an
environmental gradient, without population dynamics, competition or
dispersal. The expected abundance of species \\j\\ at site \\i\\ is
\\r_j \exp(-(E_i - \theta_j)^2 / (2 \sigma_j^2))\\, with \\r_j\\ =
`r_max`, \\\theta_j\\ = `optimum`, \\\sigma_j\\ = `sigma` and \\E_i\\
the environment, the Gaussian response of Dray & Legendre (2008), also
used by Peres-Neto et al. (2012). By default, observed abundances are
Poisson counts with that mean, as in Peres-Neto & Kembel (2015) and
Duarte et al. (2016), where `r_max` lies between 0 and 30. In the two
earlier studies `r_max` lies between 0.5 and 1, abundances are the
expected values (`sampling = "expected"`), and Dray & Legendre (2008)
draw presences from a binomial distribution with the expected value as
the probability (`sampling = "binomial"`). Because nothing but the
gradient structures the data, this is a baseline for
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
and a fast generator for evaluating methods.

## Usage

``` r
simulate_static(
  traits,
  env,
  incidence = FALSE,
  sampling = c("poisson", "binomial", "expected")
)
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

  Return presences and absences instead of Poisson counts? A species is
  present where its count is positive. Only for `sampling = "poisson"`.

- sampling:

  How to obtain the data from the expected values: `"poisson"` draws
  counts, `"binomial"` draws presences and absences, which needs
  expected values of at most 1, and `"expected"` returns the expected
  values.

## Value

List of class `mc_sim` with `abundance`, a matrix with one row per site
and one column per species, plus `env`, `traits` and, when the traits
came from
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)
with a tree, `tree`. Sites and species with no individuals are kept, so
remove them if a method requires it.

## Details

With several gradients, for example an environmental and a spatial one,
each species has its own response to each gradient and the expected
abundance is the product of the responses, as in Peres-Neto et al.
(2012).

## References

Dray, S., & Legendre, P. (2008). Testing the species traits-environment
relationships: the fourth-corner problem revisited. *Ecology*, 89,
3400–3412. [doi:10.1890/08-0349.1](https://doi.org/10.1890/08-0349.1)

Duarte, L. D. S., et al. (2016). Dissecting phylogenetic fuzzy
weighting: theory and application in metacommunity phylogenetics.
*Methods in Ecology and Evolution*, 7, 937–946.
[doi:10.1111/2041-210X.12547](https://doi.org/10.1111/2041-210X.12547)

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

## See also

Other simulate metacommunities:
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
sim <- simulate_static(traits, env = runif(100, 0, 100))
sim
#> <mc_sim> 100 sites x 20 species
#> species persisting: 19 of 20 
#> mean local richness: 8.35 

# an environmental and a spatial gradient
space <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
                        rmax_range = c(1, 1))
gradients <- data.frame(env = runif(100, 0, 100), space = 1:100)
simulate_static(list(traits, space), gradients)
#> <mc_sim> 100 sites x 20 species
#> species persisting: 18 of 20 
#> mean local richness: 4.01 
```
