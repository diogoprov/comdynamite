# Partition trait or phylogenetic variation between environment and space

The phylogenetic gradient analysis of Peres-Neto et al. (2012). Trait
values are standardised with species weights, which gives each site a
mean and a standard deviation. The chosen component is related to `env`
and `space` with site weights, and the variation is split into the
fractions of Peres-Neto et al. (2006): `a`, environment alone; `b`,
shared; `c`, space alone; and `d`, unexplained.

## Usage

``` r
trait_partition(
  comm,
  traits,
  env,
  space,
  component = c("mean", "variance"),
  n_perm = 999
)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species. Every site
  needs at least one individual.

- traits:

  Numeric vector or matrix with one row per species, which can hold
  traits or the vectors of
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).
  Species are matched to `comm` by name, in any order, and those in only
  one of the two are dropped with a warning.

- env, space:

  Numeric vector, or matrix or data frame, with one row per site, for
  example the output of
  [`spatial_eigenvectors()`](https://provetelab.org/comdynamite/reference/spatial_eigenvectors.md)
  for `space`.

- component:

  `"mean"` for the composition of sites or `"variance"` for their
  dispersion.

- n_perm:

  Number of permutations of sites and of species.

## Value

List of class `mc_trait_partition` with `fractions`, a data frame with
one row for each of `abc`, `ab`, `bc`, `a`, `b`, `c` and `d`, and the
columns `r2`, `r2_adjusted`, `r2_component`, `p_site`, `p_taxa` and
`p_max`; `component`; and `component_share`, the share of the total
variation in the component.

## Details

The total variation includes both the mean and the dispersion within
sites, so `r2` answers how much of all the variation each fraction
holds, and `r2_component` how much of the component. Fractions `b` and
`d` cannot be tested.

## Original code

Results for the mean component agree with
`PartitionRDATwoDataMatricesPhylogeny` of the authors in the raw and
adjusted fractions of the total. The code departs from the text in the
following:

- The fractions relative to the component mix the two denominators, and
  their adjusted values use the number of predictors of both sets in
  every model, which can give values above one.

- The total variation does not grow with the number of traits, as it
  should according to Appendix A.

- For the variance, the explained sum of squares contains the intercept,
  so `b` absorbs the mean dispersion. Here the intercept is removed.

- The permutations of the variance take the square root twice.

## References

Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006).
Variation partitioning of species data matrices: estimation and
comparison of fractions. *Ecology*, 87, 2614–2625.
[doi:10.1890/0012-9658(2006)87\[2614:VPOSDM\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282006%2987%5B2614%3AVPOSDM%5D2.0.CO%3B2)

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

## Examples

``` r
set.seed(1)
env <- runif(40, 0, 100)
space <- spatial_eigenvectors(cbind(1:40))[, 1:3]
trait <- setNames(runif(15, 0, 100), paste0("sp", 1:15))
comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
colnames(comm) <- names(trait)
trait_partition(comm, trait, env, space, n_perm = 49)
#> <mc_trait_partition> mean component, 0.511 of the total variation
#>   fraction       r2 r2_adjusted r2_component p_site p_taxa p_max
#> 1      abc 0.490557      0.4323      0.96059   0.02   0.02  0.02
#> 2       ab 0.489712      0.4763      0.95894   0.02   0.02  0.02
#> 3       bc 0.020911     -0.0607      0.04095   0.58   0.02  0.58
#> 4        a 0.469647      0.4930      0.91965   0.02   0.02  0.02
#> 5        b 0.020065     -0.0167      0.03929     NA     NA    NA
#> 6        c 0.000846     -0.0439      0.00166   0.98   0.06  0.98
#> 7        d 0.509443      0.5677      0.03941     NA     NA    NA
```
