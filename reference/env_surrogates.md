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
and also serve as null environments for simulations.

## Usage

``` r
env_surrogates(env, coords, n = 99, method = c("variogram", "parametric"))
```

## Arguments

- env:

  Predictor, or matrix or data frame of predictors, with one row per
  site.

- coords:

  Matrix or data frame of coordinates with one row per site.

- n:

  Number of surrogates.

- method:

  `"variogram"` or `"parametric"`.

## Value

List of `n` matrices with the dimensions and column names of `env`. For
`"parametric"`, `attr(, "rows")` is a matrix with one row per draw,
giving for each site the row of `env` it received.

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
```
