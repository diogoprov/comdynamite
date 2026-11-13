# Genetic and environmental covariance matrices for several traits

Builds one additive genetic covariance matrix \\G\\ and one
environmental covariance matrix \\E\\ per species, for the evolution of
several traits (one per niche axis) in
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
under three families of relationship between them:

## Usage

``` r
genetic_matrices(
  species,
  n_traits = 2,
  family = c("proportional", "independent", "shared"),
  P = NULL,
  correlation = 0,
  heritability = 0.5,
  divergence = 0
)
```

## Arguments

- species:

  Names of the species, or their number.

- n_traits:

  Number of traits (axes). Ignored when `P` is given.

- family:

  One of `"proportional"`, `"independent"` and `"shared"`.

- P:

  Phenotypic covariance matrix, symmetric and positive definite. By
  default, unit variances with a common correlation `correlation`.

- correlation:

  Correlation between all pairs of traits in the default `P`.

- heritability:

  Share \\h^2\\ of the phenotypic variance that is additive genetic,
  between 0 and 1.

- divergence:

  Standard deviation (on the log scale) of the differences among
  species.

## Value

A list with `genetic_var` and `env_var`, each a named list of one matrix
per species, ready for
`evolution_control(genetic_var = m$genetic_var, env_var = m$env_var)`.

## Details

- `"proportional"`: \\G = h^2 P\\ and \\E = (1 - h^2) P\\, both
  proportional to a phenotypic matrix \\P\\ (`P`) that is shared by all
  species up to a scale factor. Under neutral evolution with such
  matrices, the covariance among population means is expected to be
  proportional to that within populations (Lande 1979), the expectation
  tested with phenotypic matrices by Ackermann & Cheverud (2004) and
  whose error rates were studied by Prôa et al. (2013).

- `"independent"`: \\G\\ and \\E\\ have no common structure: each
  species draws two random orientations (eigenvectors) and two sets of
  variances, with total variances \\h^2\\ and \\1 - h^2\\ times that of
  `P`.

- `"shared"`: \\G\\ and \\E\\ share the principal components
  (eigenvectors) of `P` but differ in the variance along them, which
  each species draws around \\h^2\\ and \\1 - h^2\\ times the
  eigenvalues of `P`.

`divergence` sets how much species differ: the standard deviation of the
log of the scale factor of `P` (`"proportional"`) or of the eigenvalues
(`"independent"` and `"shared"`).

## References

Ackermann, R. R., & Cheverud, J. M. (2004). Detecting genetic drift
versus selection in human evolution. *Proceedings of the National
Academy of Sciences*, 101, 17946–17951.
[doi:10.1073/pnas.0405919102](https://doi.org/10.1073/pnas.0405919102)

Lande, R. (1979). Quantitative genetic analysis of multivariate
evolution, applied to brain:body size allometry. *Evolution*, 33,
402–416.
[doi:10.1111/j.1558-5646.1979.tb04694.x](https://doi.org/10.1111/j.1558-5646.1979.tb04694.x)

Prôa, M., O'Higgins, P., & Monteiro, L. R. (2013). Type I error rates
for testing genetic drift with phenotypic covariance matrices: a
simulation study. *Evolution*, 67, 185–195.
[doi:10.1111/j.1558-5646.2012.01746.x](https://doi.org/10.1111/j.1558-5646.2012.01746.x)

## See also

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md),
[`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md),
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

## Examples

``` r
set.seed(1)
m <- genetic_matrices(c("sp1", "sp2"), n_traits = 2, family = "shared",
                      correlation = 0.6, heritability = 0.4)
m$genetic_var$sp1
#>      [,1] [,2]
#> [1,] 0.40 0.24
#> [2,] 0.24 0.40
ev <- evolution_control(genetic_var = m$genetic_var, env_var = m$env_var)
```
