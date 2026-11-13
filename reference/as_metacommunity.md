# Convert to a simulated metacommunity

Turns the output of
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
into an `mc_sim` object, as from
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
with the extant species, their reconstructed tree and their traits (the
mean genetic optimum over their populations), for
[`as_mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md),
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
and the other functions of the package.

## Usage

``` r
as_metacommunity(x, ...)
```

## Arguments

- x:

  An `mc_evolving` object.

- ...:

  Unused.

## Value

An `mc_sim` object.

## See also

Other simulate metacommunities:
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md),
[`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md),
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

## Examples

``` r
set.seed(1)
ls <- landscape_grid(5, disp_scale = 0.5)
ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
sim <- simulate_evolving_metacommunity(ls, ancestor, dispersal = 0.002, n_steps = 300,
                                       speciation = speciation_control("geographic", 80))
as_metacommunity(sim)
#> <mc_sim> 25 sites x 86 species
#> species persisting: 86 of 86 
#> mean local richness: 6.4 
```
