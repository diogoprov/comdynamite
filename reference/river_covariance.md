# Covariance of the environment along a river network

Covariance among the sites of a river network as a sum of exponential
components, the variance-component approach of Peterson and Ver Hoef
(2010) for stream networks:

## Usage

``` r
river_covariance(
  landscape,
  tail_up = NULL,
  tail_down = NULL,
  euclidean = NULL,
  nugget = 0,
  weights = "magnitude"
)
```

## Arguments

- landscape:

  A landscape from
  [`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md).

- tail_up, tail_down, euclidean:

  `NULL` (the default) to leave the component out, or a numeric vector
  with the partial sill (variance) and the range \\\alpha\\, named
  `sill` and `range` or in that order. Ranges are in the units of the
  reach lengths (tail-up and tail-down) or of the coordinates
  (Euclidean).

- nugget:

  Variance independent among sites: one number, 0 or more.

- weights:

  Weights of the branches for the tail-up component: `"magnitude"`,
  `"upstream_length"` or a positive numeric vector with one value per
  site.

## Value

A symmetric sites x sites covariance matrix.

## Details

- `tail_up`: correlation only between flow-connected sites (one
  downstream of the other), \\\sigma^2 \exp(-h / \alpha) \sqrt{A_u /
  A_d}\\, where \\h\\ is the distance along the channel and \\A_u\\ and
  \\A_d\\ are the additive function at the upstream and downstream site
  (Ver Hoef et al. 2006). It suits variables carried by the flow, such
  as water chemistry or temperature, which mix at confluences in
  proportion to the flow of each branch.

- `tail_down`: correlation between all sites, \\\sigma^2 \exp(-h /
  \alpha)\\, which for the exponential model is the same for
  flow-connected and flow-unconnected pairs at the same distance along
  the channel (Ver Hoef and Peterson 2010). It suits processes that
  spread up and down the channel.

- `euclidean`: \\\sigma^2 \exp(-d / \alpha)\\ with the straight-line
  distance \\d\\, for drivers of the terrestrial landscape and the
  atmosphere, such as climate or land use, that ignore the channel.

The additive function comes from weights of the branches: at each
confluence the moving average of the tail-up model splits among the
branches upstream in proportion to their weights, and the additive
function of a site is the product of these shares down to the outlet,
which keeps the variance equal at all sites. With
`weights = "magnitude"` (the default) the shares follow the Shreve
magnitude, with `"upstream_length"` the length of channel upstream (a
proxy for catchment area or discharge), and a numeric vector gives any
other positive measure per site, such as discharge.

## References

Peterson, E. E., & Ver Hoef, J. M. (2010). A mixed-model moving-average
approach to geostatistical modeling in stream networks. *Ecology*, 91,
644–651. [doi:10.1890/08-1668.1](https://doi.org/10.1890/08-1668.1)

Ver Hoef, J. M., Peterson, E., & Theobald, D. (2006). Spatial
statistical models that use flow and stream distance. *Environmental and
Ecological Statistics*, 13, 449–464.
[doi:10.1007/s10651-006-0022-8](https://doi.org/10.1007/s10651-006-0022-8)

Ver Hoef, J. M., & Peterson, E. E. (2010). A moving average approach for
spatial statistical models of stream networks. *Journal of the American
Statistical Association*, 105, 6–18.
[doi:10.1198/jasa.2009.ap08248](https://doi.org/10.1198/jasa.2009.ap08248)

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
[`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md),
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)

## Examples

``` r
set.seed(1)
land <- landscape_river(river_network(30))
v <- river_covariance(land, tail_up = c(sill = 1, range = 5), nugget = 0.1)
# flow-unconnected sites are uncorrelated under the tail-up model
unconnected <- land$downstream > 0 & t(land$downstream) > 0
range(v[unconnected])
#> [1] 0 0
```
