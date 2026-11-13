# Landscape from site coordinates

Builds a spatially explicit landscape from the coordinates and
environment of real sites, such as surveyed lakes, so that simulations
share their spatial layout. Emigrants reach the other sites with
probabilities from
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md).

## Usage

``` r
landscape_points(coords, env, disp_scale, lonlat = FALSE)
```

## Arguments

- coords:

  Matrix or data frame with one row per site and two columns of
  coordinates. With `lonlat = TRUE`, longitude and latitude in decimal
  degrees.

- env:

  Numeric vector with the environment of each site.

- disp_scale:

  Scale of the exponential dispersal kernel, in the units of the
  distances (kilometres with `lonlat = TRUE`).

- lonlat:

  Are the coordinates longitude and latitude? Distances are then
  great-circle distances in kilometres, with a spherical Earth of radius
  6371 km; otherwise Euclidean.

## Value

An `mc_landscape` object, as from
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
without a torus.

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md)

## Examples

``` r
xy <- cbind(lon = c(-73, -72.5, -71), lat = c(43, 43.2, 44))
ls <- landscape_points(xy, env = c(10, 20, 80), disp_scale = 50, lonlat = TRUE)
round(ls$dist)
#>      [,1] [,2] [,3]
#> [1,]    0   46  196
#> [2,]   46    0  150
#> [3,]  196  150    0
```
