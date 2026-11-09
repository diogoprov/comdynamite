# Plot the profile of a fit of population shifts

Draws, along the width `w` of the competition kernel (on a log scale),
the \\R^2\\ of the fit, the estimate of displacement `c` and the slopes
`h` of
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
with the chosen width as a dashed line and the lower end of the default
grid (the pooled standard deviation within populations, when it sets it)
as a dotted line. A flat \\R^2\\ means that `w` is not identified, and
estimates that change along it are then not to be trusted. With a
bootstrap, the interval of `c` is drawn at the chosen width.

## Usage

``` r
# S3 method for class 'population_shifts'
plot(x, ...)
```

## Arguments

- x:

  Result of
  [`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
  with the width estimated (more than one width in its profile).

- ...:

  Not used.

## Value

A ggplot object.
