# Dispersal kernel

The probability that an emigrant from site s reaches site t is
proportional to exp(-d / `scale`) for the exponential kernel or exp(-d^2
/ (2 `scale`^2)) for the Gaussian one, where d is the distance between
the sites. Emigrants do not stay at their own site: the share that stays
is set by the emigration probability (`dispersal` in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)).
For the same `scale`, the Gaussian kernel sends emigrants to near sites
more often and to far sites much less often (a thin tail), whereas the
exponential kernel keeps more long-distance dispersal. Any other kernel
can be passed as a matrix to
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md).

## Usage

``` r
dispersal_kernel(dist, scale, shape = c("exponential", "gaussian"))
```

## Arguments

- dist:

  Square matrix of distances among sites.

- scale:

  Scale of the kernel, in the units of `dist`: the mean distance of the
  exponential kernel and the standard deviation of the Gaussian one.

- shape:

  `"exponential"` (the default) or `"gaussian"`.

## Value

Row-stochastic matrix with a zero diagonal: row s gives the
probabilities that an emigrant from site s reaches each site.

## See also

Other landscapes:
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md)

## Examples

``` r
d <- as.matrix(dist(cbind(1:5, 0)))
round(dispersal_kernel(d, scale = 1), 2)
#>      1    2    3    4    5
#> 1 0.00 0.64 0.24 0.09 0.03
#> 2 0.40 0.00 0.40 0.15 0.05
#> 3 0.13 0.37 0.00 0.37 0.13
#> 4 0.05 0.15 0.40 0.00 0.40
#> 5 0.03 0.09 0.24 0.64 0.00
round(dispersal_kernel(d, scale = 1, shape = "gaussian"), 2)
#>      1    2    3    4    5
#> 1 0.00 0.81 0.18 0.01 0.00
#> 2 0.45 0.00 0.45 0.10 0.01
#> 3 0.09 0.41 0.00 0.41 0.09
#> 4 0.01 0.10 0.45 0.00 0.45
#> 5 0.00 0.01 0.18 0.81 0.00
```
