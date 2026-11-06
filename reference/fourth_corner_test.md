# Test the association between a trait and the environment

Weighted correlation between a site predictor and the mean of a species
trait in each site, which is the fourth-corner statistic in the form of
Peres-Neto et al. (2012). Permuting sites tests the link between
composition and the predictor, and permuting species tests the link
between composition and the trait. The association needs both, so the
p-value to report is the larger of the two (ter Braak et al. 2012).

## Usage

``` r
fourth_corner_test(comm, trait, env, n_perm = 999, site_permutations = NULL)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species. Every site
  needs at least one individual.

- trait:

  Numeric vector with the value of each species. Species are matched to
  `comm` by name, in any order, and those in only one of the two are
  dropped with a warning.

- env:

  Numeric vector with the value of each site.

- n_perm:

  Number of permutations of sites and of species.

- site_permutations:

  Optional matrix of site indices with one permutation per row, for
  example from
  [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
  to use instead of free permutations.

## Value

List of class `mc_fourth_corner` with `correlation`, `p_site`, `p_taxa`
and `p_max`.

## Article and original code

The function `TraitEnvCor` of the supplement differs from Appendix A,
and this function from both:

- Appendix A stacks the two null distributions, adds one to the count
  and uses alpha / 2. The code adds two to the count. Here we use the
  larger of the two p-values. The stacked p-value is close to their
  mean, so any test that it rejects at alpha / 2 is also rejected by the
  larger p-value at alpha.

- Appendix A says that the statistic is calculated again for the
  permuted data. The code keeps the observed mean and standard deviation
  of the permuted variable, which then no longer has weighted mean zero.
  Here it is standardised again, as the code of the partition does.

- The code runs `n_perm` - 1 permutations in each null distribution.

- Appendix A allows several traits and predictors. The permutations of
  the code work with one of each, as here.

- The code returns the absolute value, and here the correlation keeps
  its sign.

## References

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

ter Braak, C. J. F., Cormont, A., & Dray, S. (2012). Improved testing of
species traits–environment relationships in the fourth-corner problem.
*Ecology*, 93, 1525–1526.
[doi:10.1890/12-0126.1](https://doi.org/10.1890/12-0126.1)

## Examples

``` r
set.seed(1)
env <- runif(40, 0, 100)
trait <- setNames(runif(15, 0, 100), paste0("sp", 1:15))
comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
colnames(comm) <- names(trait)
fourth_corner_test(comm, trait, env, n_perm = 99)
#> <mc_fourth_corner> trait-environment correlation
#> correlation: 0.7 
#> p (sites): 0.01  p (species): 0.01  p (max): 0.01 
```
