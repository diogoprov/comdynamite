# Exponential dispersal kernel

The probability that an emigrant from site s reaches site t is
proportional to exp(-d / `scale`), where d is the distance between the
sites. Emigrants do not stay at their own site.

## Usage

``` r
dispersal_kernel(dist, scale)
```

## Arguments

- dist:

  Square matrix of distances among sites.

- scale:

  Scale of the kernel, in the units of `dist`.

## Value

Row-stochastic matrix with a zero diagonal: row s gives the
probabilities that an emigrant from site s reaches each site.

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
```
