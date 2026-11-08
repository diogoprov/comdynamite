# Turnover of interaction partners in a simulated metacommunity

Counts, for each species, how often the composition around it changes:
the number of other species that enter or leave the sites where the
species is present, between successive recorded steps of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
The count is divided by the number of sites and time steps over which
the species was present at both ends of an interval, so it is a rate of
change of partners per site and per time step. It measures the entries
and exits of partners of scenario 5 of the notes on branch rates, and
feeds
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md).

## Usage

``` r
partner_turnover(x, partners = c("all", "competitors"))
```

## Arguments

- x:

  An `mc_sim` object with `history`, that is, run with `record_steps`.

- partners:

  `"all"` counts every other species as a partner; `"competitors"`
  counts only the species with a non-zero competition coefficient with
  the focal species in `x$alpha` (in either direction).

## Value

Named numeric vector, one value per species. Species never present at
two successive recorded steps get zero.

## See also

[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md)

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(10)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tr),
                              competition_matrix(tr, "none"), n_steps = 300,
                              record_steps = seq(100, 300, 20), extirpation = 0.02)
partner_turnover(sim)
#>         t1         t2         t3         t4         t5         t6         t7 
#> 0.05564516 0.08333333 0.05588235 0.09000000 0.04500000 0.07142857 0.00000000 
#>         t8         t9        t10 
#> 0.00000000 0.06093750 0.06106870 
```
