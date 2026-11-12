# Gradients of the seven scenarios of Peres-Neto et al. (2012)

Two gradients on a transect, as in Appendix C of Peres-Neto et al.
(2012). Unless stated, `x` is evenly spaced on \[-15, 15\] and `y` is
`x` in random order, so that the two are independent.

## Usage

``` r
logistic_gradients(scenario, n_sites = 31)
```

## Arguments

- scenario:

  Number of the scenario, from 1 to 7.

- n_sites:

  Number of sites.

## Value

Data frame with one row per site: `x` and `y`, the gradients that
generate the data, and `x_test` and `y_test`, the ones to analyse.

## Details

1.  The two gradients are equally important and independent.

2.  They are equally important and correlated: `y` is `x` plus a N(0, 1)
    error.

3.  They are independent, and `y` is less important: data come from `y`,
    and `y_test` is `y` with unit variance plus a N(0, 1) error. On the
    scale of `y` this error would have no effect, and with unit variance
    the results agree with Table 5.

4.  Neither is important: data come from scenario 1, and `x_test` and
    `y_test` are N(0, 1) deviates, to estimate the type I error.

5.  `y` has a broad spatial scale, the first Moran's eigenvector map
    (MEM), and `x` is in random order along the transect.

6.  `y` has an intermediate spatial scale, the third MEM, and `x` is as
    in scenario 5.

7.  `x` is the first MEM and `y` the third.

MEMs are rescaled to \[-15, 15\].

## References

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

## See also

Other simulate metacommunities:
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)

## Examples

``` r
set.seed(1)
head(logistic_gradients(6))
#>     x          y x_test     y_test
#> 1 -13  -5.740251    -13  -5.740251
#> 2 -10 -10.606602    -10 -10.606602
#> 3  -6 -13.858193     -6 -13.858193
#> 4  13 -15.000000     13 -15.000000
#> 5  14 -13.858193     14 -13.858193
#> 6  -1 -10.606602     -1 -10.606602
```
