# DAMOCLES fits over many configurations of regimes or covariates

`damocles_samples()` fits
[`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
with parameters shared by all communities, once for each configuration
in a list of `regimes` or of `edge_covariate`, and compares each fit
with constant rates by
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md).
The configurations are typically samples from the posterior distribution
of a method that estimates rates of diversification along the tree, so
that the uncertainty in where rates shift, or in branch-specific rates,
carries over to the conclusions about local colonisation and extinction
instead of resting on a single configuration.

## Usage

``` r
damocles_samples(
  comm,
  tree,
  regimes = NULL,
  edge_covariate = NULL,
  vary = c("both", "gamma", "mu"),
  root = c("both", "absent", "present"),
  sympatry = 0,
  cores = 1
)

# S3 method for class 'mc_damocles_samples'
summary(object, probs = c(0.05, 0.5, 0.95), ...)
```

## Arguments

- comm:

  Presences of the species (columns, named as the tips of `tree`) in one
  or more communities (rows): a matrix or data frame, or a named vector
  for one community. Values above zero count as presences. Tips missing
  from `comm` count as absent, with a warning, and species of `comm`
  missing from the tree are dropped.

- tree:

  A rooted, fully dichotomous `phylo` tree with branch lengths.

- regimes, edge_covariate:

  A list of configurations, each a vector with one value per edge of
  `tree` in the order of `tree$edge` (see
  [`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)).
  Give one of the two.

- vary:

  Which rates change among `regimes` or with `edge_covariate`: `"both"`,
  `"gamma"` or `"mu"`; the other one is a single rate for the whole
  tree.

- root:

  Prior on the state of the ancestor at the root: `"both"`, `"absent"`
  or `"present"`.

- sympatry:

  In
  [`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
  the probability that a present ancestor leaves both daughters present
  (see
  [`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)).
  In
  [`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
  a fixed value, or `NULL` to estimate it with the rates (shared or not,
  as they are).

- cores:

  Number of cores, with configurations fitted in parallel. Not available
  on Windows, where one core is used.

- object:

  A result of `damocles_samples()`.

- probs:

  Probabilities of the quantiles.

- ...:

  Not used.

## Value

A list of class `mc_damocles_samples` with `samples`, a data frame with
one row per configuration (`n_par`, `loglik`, `delta_aic` against
constant rates, and the likelihood-ratio `statistic`, `df` and
`p_value`), `constant`, the fit with constant rates, `fits`, the fit of
each configuration, and `gamma` and `mu`, matrices of the rate of each
edge (rows) in each configuration (columns).

## Details

With BAMM output read by `BAMMtools::getEventData()`, for example, the
regimes of sample `i` are `factor(edata$eventVectors[[i]])` and the mean
branch rates of speciation are the columns of
`BAMMtools::getMarginalBranchRateMatrix(edata)$lambda_branch_matrix`,
both in the order of `tree$edge` when `edata` was built on `tree`
(`data-raw/damocles_bamm.R` shows the whole path). With ClaDS
(`RPANDA::fit_ClaDS()`), the branch rates of one posterior sample are
`exp(l_0 + alpha * sampler$alpha_effect + x)`, where `l_0`, `alpha` and
the values `x` of the edges are columns of the chains, transformed as in
`RPANDA::getMAPS_ClaDS()` (`data-raw/damocles_clads.R`). Any other
method works the same way once its output is a vector per edge;
[`dr_rates()`](https://provetelab.org/comdynamite/reference/dr_rates.md)
gives one from the tree alone.

Configurations of regimes and of continuous rates behave differently. A
posterior sample of regimes is one plausible partition, and the spread
of the fits across samples shows how much the conclusions depend on
where the shifts are. A posterior sample of branch rates also carries
the sampling noise of each rate, and a noisy covariate attenuates its
slope (regression dilution): in `data-raw/damocles_clads.R`, with
colonisation simulated from the MAP rates of ClaDS with `beta = 2`, the
MAP rates gave a slope of 1.6 and rejected constant rates, while single
posterior samples gave a median slope of 0.22 and rejected constant
rates in 65 of 100. For continuous rates, fit the posterior mean, median
or MAP of each rate with
[`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
for the estimate, and use the posterior samples to see whether the sign
and the comparison with constant rates hold, knowing that each sample
understates the slope.

[`summary()`](https://rdrr.io/r/base/summary.html) reports, across
configurations, how often the varying rates fit better than constant
rates (p-value below 0.05, and AIC lower by more than 2), quantiles of
the parameters that mean the same in every fit (such as `mu` when only
`gamma` varies, or `beta_gamma`), and quantiles of the rate of each
edge, which are comparable whatever the regimes of each sample. Regimes
differ among posterior samples, so the parameters of each regime are
only comparable within a sample.

The configurations are only as good as the method that produced them.
The reliability of the rate shifts inferred by BAMM has been debated
(Moore et al. 2016; Rabosky et al. 2017); comparing several sources of
regimes, including partitions defined a priori with
[`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md),
and testing a configuration against random partitions with
[`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md)
guard against relying on one of them.

## References

Moore, B. R., Höhna, S., May, M. R., Rannala, B., & Huelsenbeck, J. P.
(2016). Critically evaluating the theory and performance of Bayesian
analysis of macroevolutionary mixtures. *Proceedings of the National
Academy of Sciences*, 113, 9569–9574.
[doi:10.1073/pnas.1518659113](https://doi.org/10.1073/pnas.1518659113)

Rabosky, D. L., Mitchell, J. S., & Chang, J. (2017). Is BAMM flawed?
Theoretical and practical concerns in the analysis of multi-rate
diversification models. *Systematic Biology*, 66, 477–498.
[doi:10.1093/sysbio/syx037](https://doi.org/10.1093/sysbio/syx037)

## See also

[`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md),
[`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md).

## Examples

``` r
set.seed(7)
tr <- treeshapesim::sim_tree_shape(30)
crown <- ape::Ntip(tr) + 2
comm <- simulate_damocles(tr, ifelse(clade_regimes(tr, list(A = crown)) == "A", 3, 0.5),
                          mu = 1, n = 5)
# three configurations, as if sampled from a posterior
configurations <- lapply(c(crown, crown + 1, crown + 2), function(node) {
  clade_regimes(tr, list(A = node))
})
res <- damocles_samples(comm, tr, regimes = configurations, vary = "gamma")
res
#> DAMOCLES fits over 3 configurations, parameters shared by 5 communities
#> Better than constant rates: p < 0.05 in 1 and AIC lower by more than 2 in 1 of 3 
summary(res)
#> DAMOCLES fits over 3 configurations
#> 
#> Proportion of configurations better than constant rates:
#>         measure proportion
#>  p_value < 0.05      0.333
#>  delta_aic < -2      0.333
#> 
#> Likelihood-ratio statistic:
#>   5%  50%  95% 
#> 1.68 3.35 5.72 
#> Difference in AIC (varying - constant):
#>     5%    50%    95% 
#> -3.718 -1.352  0.323 
#> 
#> Parameters shared by every fit:
#>          5%     50%  95%
#> mu 1.88e-09 4.6e-09 3.05
#> 
#> Rates of each edge: see `$edges`
```
