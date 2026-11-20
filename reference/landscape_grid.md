# Spatially explicit landscape on a square grid

Sites on a `side` x `side` grid, optionally wrapped as a torus so that
no site lies on an edge. The environment is spatially autocorrelated:
independent normal deviates are smoothed by a Gaussian kernel with a
bandwidth of `env_bandwidth` cells, and then given a uniform marginal
distribution on `env_range` by ranks. Emigrants reach the other sites
with probabilities from
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md).

## Usage

``` r
landscape_grid(
  side = 10,
  env_bandwidth = 2,
  disp_scale = 2,
  torus = TRUE,
  env_range = c(0, 100),
  disp_shape = c("exponential", "gaussian"),
  n_axes = 1,
  axis_cor = 0
)
```

## Arguments

- side:

  Number of cells per side.

- env_bandwidth:

  Bandwidth (in cells) of the Gaussian kernel that smooths the
  environment. Larger values give smoother gradients.

- disp_scale:

  Scale (in cells) of the dispersal kernel.

- torus:

  Wrap the grid as a torus?

- env_range:

  Range of the uniform marginal distribution of the environment.

- disp_shape:

  Shape of the dispersal kernel, `"exponential"` (the default) or
  `"gaussian"`; see
  [`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md).

- n_axes:

  Number of environmental axes. With more than one, each axis is a field
  smoothed in the same way, the axes are correlated by `axis_cor` before
  the conversion to ranks, and the landscape also holds `env_axes`, a
  matrix with one column per axis, which
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  uses with a list of traits when `env` is not given. The first axis is
  `env`, the same as with `n_axes = 1` and the same seed.

- axis_cor:

  Correlation between the axes: one number for all pairs, or a
  correlation matrix (`n_axes` x `n_axes`). It applies to the smoothed
  fields, so the correlation of the uniform axes (by ranks) is slightly
  lower.

## Value

An `mc_landscape` object: a list with `env`, `coords` (data frame with
columns `x` and `y`), `dist` (distances among sites, toroidal if
`torus = TRUE`), `dispersal` (row-stochastic matrix with a zero
diagonal) and, with `n_axes > 1`, `env_axes`.

## See also

[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
for site permutations that keep the spatial structure of a torus.

Other landscapes:
[`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md),
[`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md),
[`env_pulse()`](https://provetelab.org/comdynamite/reference/env_pulse.md),
[`env_space_time()`](https://provetelab.org/comdynamite/reference/env_space_time.md),
[`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md),
[`landscape_matrix()`](https://provetelab.org/comdynamite/reference/landscape_matrix.md),
[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md)

## Examples

``` r
set.seed(1)
ls <- landscape_grid(side = 8)
image(matrix(ls$env, 8, 8), asp = 1)
```
