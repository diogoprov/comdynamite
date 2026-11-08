# Likelihood-ratio test between nested DAMOCLES fits

`damocles_lrt()` compares two fits of
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
to the same communities, where `null` is a special case of
`alternative`: constant rates against rates that change among `regimes`
or with an `edge_covariate` (where `beta = 0` gives the constant model),
rates shared by all communities against separate rates, or
`sympatry = 0` against an estimated `sympatry`. Twice the difference of
the total log-likelihoods is compared with a chi-squared distribution
with as many degrees of freedom as extra parameters.

## Usage

``` r
damocles_lrt(null, alternative)
```

## Arguments

- null, alternative:

  Fits of class `mc_damocles_fit` to the same communities.

## Value

A data frame with one row per model: number of free parameters, total
log-likelihood, AIC, and, for the alternative, the statistic, degrees of
freedom and p-value.

## Details

The function counts parameters but cannot check that the models are
nested; that is up to the user. When the null value lies on the boundary
of the parameter space, as `sympatry = 0` or a rate of zero in one
regime, the chi-squared p-value is conservative. The test answers
whether the partition or covariate improves the fit to the local
presences; with partitions or covariates estimated from the same tree by
methods of diversification, a significant result is an association
between those rates and the local dynamics, not evidence of a mechanism.

## See also

[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
[`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md).

## Examples

``` r
set.seed(5)
tr <- treeshapesim::sim_tree_shape(40)
reg <- clade_regimes(tr, list(A = ape::extract.clade(tr, ape::Ntip(tr) + 2)$tip.label))
comm <- simulate_damocles(tr, gamma = 1, mu = 1, n = 5)
constant <- fit_damocles(comm, tr, shared = TRUE)
by_clade <- fit_damocles(comm, tr, shared = TRUE, regimes = reg, vary = "gamma")
damocles_lrt(constant, by_clade)
#>             n_par    loglik      AIC    statistic df   p_value
#> null            2 -91.68836 187.3767           NA NA        NA
#> alternative     3 -91.68833 189.3767 6.611732e-05  1 0.9935123
```
