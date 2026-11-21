# Seasonal cycle of the environment

Adds a seasonal cycle to the environment of each site, \\a_s \sin(2 \pi
(t - \phi) / P)\\, with an amplitude per site (for instance larger in
shallow ponds or in headwaters) and, optionally, an anomaly per year
shared by all sites (wet and dry years). The environment can already
change in time (an `env_time` matrix, e.g. from
[`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md)
with `ranks = FALSE` or
[`env_river()`](https://provetelab.org/comdynamite/reference/env_river.md)
with `n_steps`), which then gives variation among years that differs
among sites; floods can be added with
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md)
or, as losses of populations, with
[`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md).

## Usage

``` r
env_seasonal(env, n_steps = NULL, period, amplitude, phase = 0, year_sd = 0)
```

## Arguments

- env:

  Environment of each site (a vector), or a matrix (sites x steps).

- n_steps:

  Number of steps. Taken from `env` when it is a matrix.

- period:

  Steps per cycle (a year).

- amplitude:

  Amplitude of the cycle: one value, or one per site.

- phase:

  Step at which the cycle crosses its mean on the way up.

- year_sd:

  Standard deviation of an anomaly drawn for each year (each block of
  `period` steps) and added to all sites; 0 (the default) for none.

## Value

Matrix with one row per site and one column per step, for `env_time` in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
[`env_river()`](https://provetelab.org/comdynamite/reference/env_river.md),
[`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md),
[`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md),
[`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md),
[`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md),
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)

## Examples

``` r
set.seed(1)
e <- env_seasonal(c(10, 20, 30), n_steps = 100, period = 25, amplitude = c(1, 3, 6),
                  year_sd = 1)
matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
```
