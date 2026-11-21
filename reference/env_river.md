# Environment along a river network

Draws environments for the sites of a river network from a Gaussian
field with the covariance of
[`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md),
for
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md):
a single value per site or, with `n_steps`, a series that changes in
time as a first-order autoregressive process with the same spatial
covariance at every step, as `env_time`.

## Usage

``` r
env_river(landscape, n_steps = NULL, mean = 0, temporal_scale = Inf, ...)
```

## Arguments

- landscape:

  A landscape from
  [`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md).

- n_steps:

  `NULL` (the default) for one value per site, or the number of steps of
  a series.

- mean:

  Mean of the environment (one value, or one per site, such as a
  gradient from the outlet to the headwaters).

- temporal_scale:

  Steps over which the correlation in time falls by a factor of \\e\\
  (the correlation between consecutive steps is \\\exp(-1 /\\
  `temporal_scale`)); `Inf` keeps the first draw at every step and 0
  draws every step anew.

- ...:

  Components of the covariance, passed to
  [`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md)
  (`tail_up`, `tail_down`, `euclidean`, `nugget` and `weights`).

## Value

A numeric vector with one value per site, or with `n_steps` a matrix of
sites x steps.

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
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
land <- landscape_river(river_network(40))
# water temperature: carried by the flow, warmer towards the outlet
temp <- env_river(land, tail_up = c(4, 10), mean = 25 - land$network$distance / 2)
land$env <- temp
# a series that changes over about 50 steps, for env_time
series <- env_river(land, n_steps = 200, tail_up = c(4, 10), nugget = 0.5,
                    mean = 25 - land$network$distance / 2, temporal_scale = 50)
dim(series)
#> [1]  40 200
```
