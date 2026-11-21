# Landscape on a river network

Builds a landscape whose sites lie on a dendritic river network, with
dispersal along the channel that can differ upstream and downstream, as
in the riverine metacommunities of Carrara et al. (2012) and Altermatt
(2013). The probability that an emigrant from site s reaches site t is
proportional to \$\$\exp\left(-\frac{d\_{\downarrow}}{s\_{\downarrow}} -
\frac{d\_{\uparrow}}{s\_{\uparrow}} \right),\$\$ where \\d\_\downarrow\\
is the distance it travels downstream (from s to the confluence where
the paths of s and t meet) and \\d\_\uparrow\\ the distance it then
travels upstream (from the confluence to t). A `scale_up` smaller than
`scale_down` gives the drift of passive dispersers, and equal scales a
symmetric kernel along the channel. A share `overland` of emigrants
instead moves in a straight line, with an exponential kernel of the
straight-line distance (scale `overland_scale`), for organisms with an
aerial or terrestrial stage that can cross between branches.

## Usage

``` r
landscape_river(
  network = river_network(50),
  env = NULL,
  scale_down = 2,
  scale_up = scale_down,
  overland = 0,
  overland_scale = 1,
  env_range = c(0, 100)
)
```

## Arguments

- network:

  A data frame with one row per site and the columns `downstream` (row
  of the downstream neighbour, `NA` for the single outlet) and `length`
  (of the reach to the downstream neighbour), and optionally `x` and
  `y`; for instance from
  [`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md).
  Without `x` and `y`, coordinates are laid out as in
  [`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md).

- env:

  Numeric vector with the environment of each site, or `NULL` (the
  default) for the distance from the outlet along the channel rescaled
  to `env_range`, a gradient from the lowlands to the headwaters.

- scale_down, scale_up:

  Scales of the kernel for the distance travelled downstream and
  upstream, in the units of `length`.

- overland:

  Share of emigrants that move in a straight line rather than along the
  channel (between 0 and 1).

- overland_scale:

  Scale of the exponential kernel of straight-line dispersal.

- env_range:

  Range of the environment when `env` is `NULL`.

## Value

An `mc_landscape` object, as from
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
whose `dist` holds the distances along the channel and which also keeps
`network` (the network completed as in
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md))
and `downstream` (sites x sites matrix of the distance travelled
downstream from the site in the row to the site in the column; its
transpose is the distance travelled upstream).

## References

Altermatt, F. (2013). Diversity in riverine metacommunities: a network
perspective. *Aquatic Ecology*, 47, 365–377.
[doi:10.1007/s10452-013-9450-3](https://doi.org/10.1007/s10452-013-9450-3)

Carrara, F., Altermatt, F., Rodriguez-Iturbe, I., & Rinaldo, A. (2012).
Dendritic connectivity controls biodiversity patterns in experimental
metacommunities. *Proceedings of the National Academy of Sciences*, 109,
5761–5766.
[doi:10.1073/pnas.1119651109](https://doi.org/10.1073/pnas.1119651109)

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
[`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md),
[`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md),
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)

## Examples

``` r
set.seed(1)
land <- landscape_river(river_network(40), scale_down = 3, scale_up = 1)
land
#> <mc_landscape> 40 sites; river network (9 headwaters) 
#> environment:    0 to 100 
# emigrants from a headwater go mostly downstream
head(order(-land$dispersal[which.max(land$network$distance), ]))
#> [1] 18 11  9 16  6 19
traits <- data.frame(optimum = c(20, 50, 80), sigma = 20, r_max = 4,
                     row.names = paste0("sp", 1:3))
alpha <- matrix(0.05, 3, 3, dimnames = list(rownames(traits), rownames(traits)))
sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.1, n_steps = 100)
```
