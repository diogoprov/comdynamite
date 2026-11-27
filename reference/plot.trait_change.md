# Plot an eco-evolutionary partition of trait change

Draws, for each pair of steps, the components of the change of the
community-weighted mean trait (averaged over sites) as stacked bars,
with the total change as a point; or the same for the growth rate of
populations. Components that lower the mean are drawn below zero. Needs
the package ggplot2.

## Usage

``` r
# S3 method for class 'trait_change'
plot(x, level = c("community", "growth"), colours = NULL, ...)
```

## Arguments

- x:

  Output of
  [`trait_change()`](https://provetelab.org/comdynamite/reference/trait_change.md).

- level:

  `"community"` (the default) or `"growth"`.

- colours:

  Colours of the components, in the order in which they are drawn:
  sorting, turnover, genetic, plastic and covariance (or within and
  covariance without genetic means) for `"community"`; genetic, plastic
  and environment for `"growth"`. `NULL` uses a palette safe for colour
  blindness.

- ...:

  Not used.

## Value

A ggplot object, which can be changed further with ggplot2 functions.

## Examples

``` r
d <- data.frame(site = "a", species = rep(c("x", "y"), 3), step = rep(1:3, each = 2),
                abundance = c(10, 10, 5, 15, 2, 18), trait = c(1, 3, 1.5, 3, 2, 3.5))
plot(trait_change(d))
#> trait_change() implements the method of:
#>   Ellner, S. P., Geber, M. A., & Hairston, N. G. (2011). Does rapid evolution matter? Measuring the rate of contemporary evolution and its impacts on ecological dynamics. Ecology Letters, 14, 603-614. https://doi.org/10.1111/j.1461-0248.2011.01616.x
#>   Hairston, N. G., Ellner, S. P., Geber, M. A., Yoshida, T., & Fox, J. A. (2005). Rapid evolution and the convergence of ecological and evolutionary time. Ecology Letters, 8, 1114-1127. https://doi.org/10.1111/j.1461-0248.2005.00812.x
#>   Price, G. R. (1970). Selection and covariance. Nature, 227, 520-521. https://doi.org/10.1038/227520a0
#>   Lepš, J., de Bello, F., Šmilauer, P., & Doležal, J. (2011). Community trait response to environment: disentangling species turnover vs intraspecific trait variability effects. Ecography, 34, 856-863. https://doi.org/10.1111/j.1600-0587.2010.06904.x
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
```
