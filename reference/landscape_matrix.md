# Landscape from a dispersal matrix

Builds a spatially explicit landscape from a given matrix of movement
among sites, for connectivity that distances alone do not describe: a
kernel of resistance (least-cost or circuit) distances, estimated
movement rates, or a network of ponds linked by corridors. Rows are
rescaled to sum to 1 and the diagonal is ignored, so the matrix may hold
relative weights.

## Usage

``` r
landscape_matrix(dispersal, env, coords = NULL, dist = NULL)
```

## Arguments

- dispersal:

  Square, non-negative matrix: row s gives the relative probabilities
  that an emigrant from site s reaches each other site. Every site needs
  at least one destination. For resistance distances `r`,
  `dispersal_kernel(r, scale)` gives an exponential kernel.

- env:

  Numeric vector with the environment of each site.

- coords:

  Optional matrix or data frame with two columns of site coordinates,
  kept for plotting and for functions that need them.

- dist:

  Optional matrix of distances among sites (for instance, the resistance
  distances behind `dispersal`), kept with the landscape.

## Value

An `mc_landscape` object, as from
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md).

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md)

## Examples

``` r
# four ponds in a row, where a road cuts the movement between ponds 2 and 3
resistance <- as.matrix(dist(1:4))
resistance[1:2, 3:4] <- resistance[1:2, 3:4] + 10
resistance[3:4, 1:2] <- resistance[3:4, 1:2] + 10
ls <- landscape_matrix(dispersal_kernel(resistance, scale = 1), env = c(10, 20, 60, 70),
                       dist = resistance)
round(ls$dispersal, 3)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    0    0
#> [2,]    1    0    0    0
#> [3,]    0    0    0    1
#> [4,]    0    0    1    0
```
