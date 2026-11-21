# Environment autocorrelated in space and time

Simulates the changing environment of Thompson et al. (2020): a Gaussian
random field over space and time with exponential covariance, held
constant during a burn-in, rescaled to \[0, 1\] and converted to ranks.
Unlike
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
which gives each site its own AR(1) series, sites close together change
together.

## Usage

``` r
env_space_time(
  landscape,
  n_steps,
  burn_in = 200,
  spatial_scale = 50,
  temporal_scale = 500,
  ranks = TRUE,
  min_range = 0.6,
  n_features = 2000
)
```

## Arguments

- landscape:

  An `mc_landscape` with coordinates, from
  [`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md)
  or
  [`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md),
  or a two-column matrix or data frame of site coordinates.

- n_steps:

  Number of steps, burn-in included: one whole number, larger than
  `burn_in`.

- burn_in:

  Number of initial steps with a constant environment: one whole number,
  0 or more.

- spatial_scale, temporal_scale:

  Scales \\a_s\\ (in units of the coordinates) and \\a_t\\ (in steps) of
  the exponential covariance: each one number greater than 0.

- ranks:

  `TRUE` or `FALSE`. If `TRUE`, convert the environment to ranks over
  all sites and steps, as in their code; if `FALSE`, keep the rescaled
  field.

- min_range:

  Minimum range of the rescaled environment at the first step after the
  burn-in; the field is redrawn until it is exceeded: one number in \[0,
  1). Use 0 to accept the first draw.

- n_features:

  Number of cosine waves of the spectral method: one whole number, 1 or
  more.

## Value

Matrix with one row per site and one column per step, with values in
\[0, 1\], for the `env_time` argument of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
Multiply it to put it on the scale of the traits.

## Details

The covariance between two site-steps is \\\exp(-\sqrt{(d/a_s)^2 +
(\Delta t/a_t)^2})\\, with \\d\\ the distance between the sites,
\\\Delta t\\ the number of steps between them, and \\a_s\\ and \\a_t\\
the spatial and temporal scales. The authors' code (Meta_com_framework
v1.0, Zenodo 10.5281/zenodo.3708350) draws this field with
`RandomFields::RFsimulate(RMexp(scale = 500))` on coordinates multiplied
by 10 and steps in units of one, which gives the defaults \\a_s = 50\\
and \\a_t = 500\\ for patches with coordinates between 1 and 100.
RandomFields is no longer on CRAN, so the field is drawn here by the
spectral method (random Fourier features): the sum of `n_features`
cosine waves whose frequencies are drawn from the spectral density of
the exponential covariance, a multivariate Student t with one degree of
freedom. The covariance is exact in expectation, and the field is close
to Gaussian for large `n_features`.

As in their code, the field is drawn for the steps after the burn-in,
the environment of the first step after it is repeated over the burn-in,
the whole field is rescaled to \[0, 1\] and redrawn until the
environment of that first step spans more than `min_range`, and the
values of all sites and steps, burn-in included, are converted to ranks,
so that the environment is uniform on \[0, 1\]. The field is not on a
torus, even when the landscape is.

## References

Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik,
R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A
process-based metacommunity framework linking local and regional scale
community ecology. *Ecology Letters*, 23, 1314–1329.
[doi:10.1111/ele.13568](https://doi.org/10.1111/ele.13568)

## See also

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
[`env_river()`](https://provetelab.org/comdynamite/reference/env_river.md),
[`env_seasonal()`](https://provetelab.org/comdynamite/reference/env_seasonal.md),
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md),
[`landscape_river()`](https://provetelab.org/comdynamite/reference/landscape_river.md),
[`river_covariance()`](https://provetelab.org/comdynamite/reference/river_covariance.md),
[`river_floods()`](https://provetelab.org/comdynamite/reference/river_floods.md),
[`river_network()`](https://provetelab.org/comdynamite/reference/river_network.md)

## Examples

``` r
set.seed(1)
xy <- cbind(x = runif(30, 1, 100), y = runif(30, 1, 100))
e <- env_space_time(xy, n_steps = 600, burn_in = 100)
matplot(t(e[1:5, ]), type = "l", lty = 1, xlab = "step", ylab = "environment")
```
