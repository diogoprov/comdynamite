# Seasonal floods along a river network

Probabilities that each site of a river network loses all its
populations at each step because of floods, for `disturbance` in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
Each period (a year) has a flood season of `flood_steps` steps, starting
at step `start` of the period, and a flood happens in a given year with
probability `year_prob`, so that years differ. During a flood, the
probability of washout at a site is proportional to its size relative to
the largest site raised to `size_effect`, scaled so that the most
exposed site has probability `probability`: positive values make floods
harsher in large rivers, negative values in headwaters (flashier floods
in small, steep streams), and 0 the same everywhere.

## Usage

``` r
river_floods(
  landscape,
  n_steps,
  period,
  flood_steps = 1,
  start = 1,
  probability = 0.5,
  size_effect = 1,
  size = "upstream_length",
  year_prob = 1
)
```

## Arguments

- landscape:

  A landscape from
  [`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md).

- n_steps:

  Number of steps: one whole number, 1 or more.

- period:

  Steps per year: one whole number, 1 or more.

- flood_steps:

  Steps of the flood season in each year: one whole number, from 1 up to
  `period`.

- start:

  Step of the year at which the flood season begins: one whole number,
  from 1 up to `period`.

- probability:

  Probability of washout during a flood at the most exposed site: one
  number between 0 and 1.

- size_effect:

  Exponent of the relative size of each site: one finite number.

- size:

  Size of each site: `"upstream_length"` (the default, a proxy for
  catchment area and discharge), `"magnitude"` (Shreve magnitude) or a
  positive numeric vector with one value per site.

- year_prob:

  Probability that a flood happens in a year: one number between 0 and
  1.

## Value

A matrix of probabilities with one row per site and one column per step,
with an attribute `flood_years`, the years in which a flood happened.

## Details

Floods here remove populations; to change the environment during floods
instead, use
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md)
on the sites and steps of the flood.

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
[`env_river()`](https://provetelab.org/comdynamite/reference/env_river.md),
[`env_seasonal()`](https://provetelab.org/comdynamite/reference/env_seasonal.md),
[`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md),
[`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md),
[`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md),
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)

## Examples

``` r
set.seed(1)
land <- landscape_river(river_network(30))
floods <- river_floods(land, n_steps = 200, period = 20, flood_steps = 2,
                       probability = 0.3, year_prob = 0.7)
attr(floods, "flood_years")
#> [1] 1 3 4 5 7 8 9
# the mainstem floods more than the headwaters
tapply(apply(floods, 1, max), land$network$headwater, mean)
#>      FALSE       TRUE 
#> 0.08136869 0.01041697 
traits <- data.frame(optimum = 50, sigma = 30, r_max = 3, row.names = "sp1")
sim <- simulate_metacommunity(land, traits, matrix(0.05, dimnames = list("sp1", "sp1")),
                              n_steps = 200, disturbance = floods)
```
