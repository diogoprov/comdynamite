# Profile-likelihood intervals of DAMOCLES parameters

Confidence intervals for the parameters of
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
from the profile likelihood: the interval holds the values at which the
log-likelihood, maximised over the other parameters, is within
`qchisq(level, 1) / 2` of its maximum. With shared parameters there is
one interval per parameter; otherwise one per community. A bound the
profile does not reach within a factor of `exp(10)` of the estimate
(rates), 10 units (slopes `beta`) or 10 units on the logit scale
(`sympatry`) is reported as the edge of the parameter space (0 or `Inf`
for rates, `-Inf` or `Inf` for slopes, 0 or 1 for `sympatry`): the data
do not bound it.

## Usage

``` r
# S3 method for class 'mc_damocles_fit'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A fit of class `mc_damocles_fit`.

- parm:

  Names of the parameters, among the columns of the fit; by default all
  free parameters.

- level:

  Confidence level.

- ...:

  Not used.

## Value

A data frame with the community (`"shared"` for shared parameters), the
parameter, its estimate and the lower and upper bounds.

## See also

[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md).

## Examples

``` r
set.seed(6)
tr <- treeshapesim::sim_tree_shape(40)
comm <- simulate_damocles(tr, gamma = 1, mu = 1, n = 10)
fit <- fit_damocles(comm, tr, shared = TRUE)
confint(fit)
#>   community parameter     estimate     lower     upper
#> 1    shared     gamma 6.814352e-01 0.5287597 0.8632444
#> 2    shared        mu 2.344120e-09 0.0000000       Inf
```
