# Colonisation history of a simulated metacommunity

The arrival step, number of founders, immigrants since arrival and
migrant load of every population present at the last step of a
simulation, recorded by
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
(see its section Colonisation history).

## Usage

``` r
colonisation_history(x)
```

## Arguments

- x:

  An `mc_sim` object from
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

## Value

Data frame with one row per population present at the last step and the
columns `site`, `species`, `abundance`, `arrival`, `age`
(`n_steps - arrival`), `propagule`, `immigrants` and `migrant_load`,
ordered by site and species.

## See also

[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
whose `founder_sd` and `monopolization` use this history.

Other simulate metacommunities:
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(10)
sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree),
                              competition_matrix(tree, "phylo"), n_steps = 100,
                              extirpation = 0.01)
head(colonisation_history(sim))
#>    site species abundance arrival age propagule immigrants migrant_load
#> 1 site1      t5         1     100   0         1          0    0.0000000
#> 2 site1      t9         1     100   0         1          0    0.0000000
#> 3 site1     t10        17      48  52         1          7    0.6214777
#> 4 site2      t5        99      12  88         1          7    0.2206160
#> 5 site2      t6        48      92   8         1          2    0.5000000
#> 6 site3      t3         1     100   0         3          0    0.0000000
```
