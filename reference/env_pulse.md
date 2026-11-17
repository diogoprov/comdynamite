# Pulse disturbance of the environment

Sets the environment of some sites to a new value from step `start` for
`duration` steps, after which it returns to its previous course, as the
pulse disturbances of Hajian-Forooshani & Chase (2026). Disturbance then
acts on growth rates through the niche of each species, not by removing
individuals as `disturbance` in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
does. Their disturbed environment (0.1) lies below the optima of all
species (0.5 to 1), so that no species grows well in disturbed sites; a
value inside the range of optima instead shifts which species are
favoured, and `duration = Inf` makes the change permanent.

## Usage

``` r
env_pulse(env, n_steps = NULL, sites, start, duration = Inf, value)
```

## Arguments

- env:

  Environment of each site (a vector, e.g. the `env` of an
  `mc_landscape`), or a matrix (sites x steps) with the environment at
  each step, e.g. from
  [`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md).

- n_steps:

  Number of time steps. Taken from `env` when it is a matrix.

- sites:

  Disturbed sites: indices or a logical vector with one value per site.
  With no site (e.g. `integer(0)`), `env` is returned unchanged, as an
  undisturbed control.

- start:

  First disturbed step.

- duration:

  Number of disturbed steps; `Inf` for a permanent change.

- value:

  Environment of the disturbed sites: one value, or one per disturbed
  site.

## Value

Matrix with one row per site and one column per time step, for the
`env_time` argument of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

## References

Hajian-Forooshani, Z., & Chase, J. M. (2026). Disturbance interacts with
dispersal and niche breadth to shape scale-dependent diversity change in
metacommunities. *The American Naturalist*, 207, 838–850.
[doi:10.1086/739955](https://doi.org/10.1086/739955)

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md)

## Examples

``` r
set.seed(1)
env <- runif(6, 0.5, 1)
e <- env_pulse(env, n_steps = 30, sites = 1:3, start = 10, duration = 10, value = 0.1)
matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
```
