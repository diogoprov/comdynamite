# Spatially implicit landscape

Sites without coordinates, with environmental values drawn from a
uniform distribution. Dispersal among sites is uniform: emigrants can
reach any site, so sites are exchangeable and permuting them is a valid
null model.

## Usage

``` r
landscape_implicit(n_sites, env_range = c(0, 100))
```

## Arguments

- n_sites:

  Number of sites: one whole number, 2 or more.

- env_range:

  Range of the uniform distribution of the environment: two finite
  numbers, the first smaller.

## Value

An `mc_landscape` object: a list with `env` (one value per site), and
`coords`, `dist` and `dispersal` set to `NULL`.

## See also

[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md)
for a spatially explicit landscape.

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
[`env_river()`](https://provetelab.org/comdynamite/reference/env_river.md),
[`env_seasonal()`](https://provetelab.org/comdynamite/reference/env_seasonal.md),
[`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md),
[`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md),
[`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md),
[`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md),
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)

## Examples

``` r
set.seed(1)
landscape_implicit(20)
#> <mc_landscape> 20 sites; spatially implicit, uniform dispersal 
#> environment:   6.18 to 99.19 
```
