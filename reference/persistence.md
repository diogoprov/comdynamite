# Metapopulation persistence

Estimates, for each species, the rates of colonisation \\c\\ and local
extirpation \\e\\ of a two-state Markov model of presence in sites over
time, and from them the two axes of metapopulation persistence of
Januario et al. (2025): the equilibrium occupancy \\\epsilon =
\mathrm{logit}(c / (c + e))\\ and the resilience \\\tau = \log(c + e)\\,
the log of the rate of turnover. The estimates are the maximum
likelihood ones: \\c\\ is the share of the site-intervals in which an
absent species arrived, and \\e\\ the share of those in which a present
one was lost.

## Usage

``` r
persistence(x, window = NULL, ...)
```

## Arguments

- x:

  An `mc_evolving` or `mc_sim` object, or an array (sites x species x
  times) of abundances or presences.

- window:

  Number of the last intervals (generations, or intervals between
  records) to use; `NULL` for all.

- ...:

  Unused.

## Value

A data frame with one row per species (lineage, for `mc_evolving`):
`species`, the number of `intervals` used, `colonisation`,
`extirpation`, the equilibrium occupancy \\c / (c + e)\\ as `occupancy`,
`epsilon` and `tau`. For `mc_evolving`, also `extant` and the good
species each extant lineage belongs to (`good_species`). Rates are `NA`
without any interval in the corresponding state, and `epsilon` is
infinite when one of the rates is 0.

## Details

For
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
transitions come from its `series`, between consecutive generations, for
every lineage (extinct ones included, with the sites lost in the
generation they died out); a site that a mother passed to a new species
counts neither as an extirpation nor as exposure. For
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
they come from the abundances at `record_steps` (`history`), between
consecutive records, which should be evenly spaced (seeding from the
pool, `seed_steps`, counts as colonisation). An array of sites x species
x times can be given directly.

## References

Januario, M., Pinsky, M. L., & Rabosky, D. L. (2025). The metapopulation
bridge to macroevolutionary speciation rates: a conceptual framework and
empirical test. *Ecology Letters*, 28, e70021.
[doi:10.1111/ele.70021](https://doi.org/10.1111/ele.70021)

## See also

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
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
sim <- simulate_evolving_metacommunity(ls, ancestor, dispersal = 0.01, extirpation = 0.02,
                                       n_steps = 200,
                                       speciation = speciation_control("fission", rate = 0.02))
head(persistence(sim, window = 50))
#>   species intervals colonisation extirpation occupancy     epsilon       tau
#> 1     sp1        50   0.07993967  0.09466437 0.4578340 -0.16906551 -1.745234
#> 2     sp2        50   0.09105691  0.08346457 0.5217519  0.08706251 -1.745707
#> 3     sp3        50   0.09843750  0.08789386 0.5282927  0.11329183 -1.680229
#> 4     sp4        50   0.08345979  0.06643952 0.5567723  0.22807282 -1.897791
#> 5     sp5         0           NA          NA        NA          NA        NA
#> 6     sp6        50   0.04871795  0.07659574 0.3887680 -0.45249400 -2.076935
#>   extant good_species
#> 1   TRUE          sp1
#> 2   TRUE          sp2
#> 3   TRUE          sp3
#> 4   TRUE          sp4
#> 5  FALSE         <NA>
#> 6   TRUE          sp6
```
