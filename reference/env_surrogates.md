# Spatial surrogates of predictors

Draws predictors with the same spatial autocorrelation as the observed
ones but no relation to anything else, as null models for tests in
landscapes with edges, where
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
do not apply or reject too often. They are the `null` options of
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
and
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
which also take the list returned here (for example drawn with `lonlat`
or `strata`), and also serve as null environments for simulations.

## Usage

``` r
env_surrogates(
  env,
  coords,
  n = 99,
  method = c("variogram", "parametric"),
  lonlat = FALSE,
  strata = NULL
)
```

## Arguments

- env:

  Predictor, or matrix or data frame of predictors, with one row per
  site, such as an `env_time` matrix with one column per step. Or an
  `mc_landscape`, whose environment (or `env_axes`) is replaced, with
  the distances of the landscape.

- coords:

  Matrix or data frame of coordinates with one row per site. Not needed
  for a landscape.

- n:

  Number of surrogates, one whole number, 1 or more.

- method:

  `"variogram"` or `"parametric"`.

- lonlat:

  Are `coords` longitude and latitude in decimal degrees (`TRUE` or
  `FALSE`)? Distances are then great-circle distances in kilometres.

- strata:

  Optional group of each site, such as the city or the basin where it
  lies. Surrogates are drawn within each group on its own, so values
  never move between groups, which have their own autocorrelation. Each
  group needs at least 5 sites; a group where the predictors do not vary
  keeps them.

## Value

List of `n` matrices with the dimensions and column names of `env`, or
of `n` landscapes when `env` is one. For `"parametric"`,
`attr(, "rows")` is a matrix with one row per draw, giving for each site
the row of `env` it received.

## Details

`"variogram"` follows Viladomat et al. (2014): each predictor is
permuted, smoothed with an exponential kernel over a proportion of its
nearest sites, with the proportion chosen so that the variogram of the
smoothed values best matches the observed one up to the 25% quantile of
the distances between sites, and then rescaled, with independent noise
for the nugget, so that the scale and intercept of the variogram match
too. We use binned variograms rather than their kernel-smoothed one.
With several predictors, all of them share the permutation of each draw,
but each one gets its own smoothing, so their correlation is kept only
approximately.

`"parametric"` fits an exponential range by maximum likelihood to the
normal scores of the predictor, draws a Gaussian field with that range
and gives the observed values to the sites in the order of the field, so
each draw is a permutation of sites (in `attr(, "rows")`). With several
predictors, the order is that of the first principal component of the
scaled predictors, and whole rows move together, which keeps their
correlation.

In `data-raw/validation_continuous_nulls.R`, for the site test of
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
with one predictor and no relation, on a grid of 20 x 20 sites, an
L-shaped region of 280 random points and a strip of 50 x 8 sites, with
exponential ranges of 2, 5 and 10 sites, the rate of rejection at the 5%
level was 0.025-0.09 (mean 0.062) for `"variogram"` and 0.03-0.10 (mean
0.066) for `"parametric"`, against 0.24-0.79 for free permutations and
0.055-0.14 for the torus on the grid. Both became liberal (0.08-0.10)
when the range approached the size of the landscape, when few
independent patches of environment are left, so with broad fields p
values near 0.05 deserve caution. These rates are for one predictor;
with several, and in
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
and
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
they were not checked.

For simulations, `env_surrogates()` gives null environments: a
landscape, or an `env_time` matrix, with the same spatial structure but
no relation to the environment that shaped the observed communities. To
ask whether a simulated pattern depends on where the environment is,
rather than only on how much of it is autocorrelated, run the same
simulation on the observed landscape and on its surrogates and compare a
statistic. With `env_time`, `"parametric"` moves the whole series of
each site to another site, so the dynamics of each site are kept and
only their place changes; `"variogram"` smooths each step on its own,
which does not keep the temporal autocorrelation. Landscapes keep their
own distances, such as those of a torus, though the exponential
covariance may then not be positive definite.

## References

Viladomat, J., Mazumder, R., McInturff, A., McCauley, D. J., & Hastie,
T. (2014). Assessing the significance of global and local correlations
under spatial autocorrelation: a nonparametric approach. *Biometrics*,
70, 409–418.
[doi:10.1111/biom.12139](https://doi.org/10.1111/biom.12139)

## See also

[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md).

Other restricted null models:
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md),
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md),
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)

## Examples

``` r
set.seed(1)
coords <- as.matrix(expand.grid(1:12, 1:12))
field <- t(chol(exp(-as.matrix(dist(coords)) / 4))) %*% rnorm(144)
draws <- env_surrogates(field, coords, n = 3)
length(draws)
#> [1] 3
cor(field, draws[[1]])
#>           [,1]
#> [1,] 0.1758885
# null landscapes for the simulator
grid <- landscape_grid(10, env_range = c(0, 100))
nulls <- env_surrogates(grid, n = 2, method = "parametric")
cor(grid$env, nulls[[1]]$env)
#> [1] -0.2500471
```
