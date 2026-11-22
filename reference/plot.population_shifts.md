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

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = 2)
fit <- population_shifts(ind)
#> population_shifts() implements the method of:
#>   Goldberg, E. E., & Lande, R. (2006). Ecological and reproductive character displacement on an environmental gradient. Evolution, 60, 1344-1357. https://doi.org/10.1554/05-696.1
#>   Adams, D. C., & Collyer, M. L. (2007). Analysis of character divergence along environmental gradients and other covariates. Evolution, 61, 510-515. https://doi.org/10.1111/j.1558-5646.2007.00063.x
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
plot(fit)
```
