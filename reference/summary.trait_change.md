# Summarise an eco-evolutionary partition of trait change

Averages the partition of
[`trait_change()`](https://provetelab.org/comdynamite/reference/trait_change.md)
over sites (community level) and over populations (population level) for
each pair of steps.

## Usage

``` r
# S3 method for class 'trait_change'
summary(object, ...)
```

## Arguments

- object:

  Output of
  [`trait_change()`](https://provetelab.org/comdynamite/reference/trait_change.md).

- ...:

  Not used.

## Value

A list with `community`, a data frame with one row per pair of steps,
the number of sites with both steps occupied (`sites`) and the mean of
each component over those sites, and `growth`, the same over populations
(with `populations` and the mean absolute value of each component,
`abs_genetic`, `abs_plastic` and `abs_environment`) or `NULL`.

## Examples

``` r
d <- data.frame(site = "a", species = rep(c("x", "y"), 2), step = rep(1:2, each = 2),
                abundance = c(10, 10, 5, 15), trait = c(1, 3, 1.5, 3))
summary(trait_change(d))
#> $community
#>   from to sites total sorting turnover within covariance genetic plastic
#> 1    1  2     1 0.625     0.5        0   0.25     -0.125      NA      NA
#>   cov_genetic cov_plastic
#> 1          NA          NA
#> 
#> $growth
#> NULL
#> 
```
