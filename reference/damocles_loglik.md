# Likelihood of local presences under the DAMOCLES model

`damocles_loglik()` computes the exact log-likelihood of the presences
and absences of the species of a clade in a local community, given their
tree and the rates of local colonisation `gamma` and extinction `mu`,
under the model of
[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md).
It follows the pruning algorithm of Pigot & Etienne (2015, Box 1): from
the tips to the root, the probabilities of the observed states below
each branch, given its state, are carried along the branch by the
analytical solution of the backward equations of the two-state process,
and combined at each node by the rule of speciation (allopatric, or with
probability `sympatry` leaving both daughters present). `fit_damocles()`
finds the maximum-likelihood rates, and optionally `sympatry`.

## Usage

``` r
damocles_loglik(
  comm,
  tree,
  gamma,
  mu,
  root = c("both", "absent", "present"),
  sympatry = 0
)

fit_damocles(
  comm,
  tree,
  root = c("both", "absent", "present"),
  start = NULL,
  shared = FALSE,
  sympatry = 0,
  regimes = NULL,
  edge_covariate = NULL,
  vary = c("both", "gamma", "mu")
)
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

- gamma:

  Rate of local colonisation per lineage, in the units of the branch
  lengths: one rate, or one per edge of `tree` in the order of
  `tree$edge`.

- mu:

  Rate of local extinction per lineage, one or one per edge. With rates
  per edge and a root edge, the two edges leaving the root must have the
  same rates, which the root edge takes.

- root:

  Prior on the state of the ancestor at the root: `"both"`, `"absent"`
  or `"present"`.

- sympatry:

  In `damocles_loglik()`, the probability that a present ancestor leaves
  both daughters present (see
  [`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)).
  In `fit_damocles()`, a fixed value, or `NULL` to estimate it with the
  rates (shared or not, as they are).

- start:

  Optional starting rates, `c(gamma = , mu = )`; by default the
  optimisation starts from four combinations of 0.2 and 2 divided by the
  height of the tree.

- shared:

  If `TRUE`, one set of parameters shared by all communities, estimated
  from the sum of their log-likelihoods; if `FALSE` (the default),
  separate parameters for each community.

- regimes:

  Optional regime of each edge of `tree`, in the order of `tree$edge` (a
  factor or character vector), such as from
  [`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md)
  or from a configuration of rate shifts estimated by another method:
  the rates in `vary` are estimated separately for each regime.

- edge_covariate:

  Optional positive number for each edge of `tree`, in the order of
  `tree$edge`, such as branch-specific rates of speciation estimated by
  another method: the rates in `vary` change with it as \\\mathrm{rate}
  \cdot (x / \bar{x}\_g)^\beta\\, where \\x\\ is the covariate and
  \\\bar{x}\_g\\ its geometric mean, with `rate` and \\\beta\\
  estimated. Give `regimes` or `edge_covariate`, not both.

- vary:

  Which rates change among `regimes` or with `edge_covariate`: `"both"`,
  `"gamma"` or `"mu"`; the other one is a single rate for the whole
  tree.

## Value

`damocles_loglik()`: the log-likelihood of each community.
`fit_damocles()`: a data frame of class `mc_damocles_fit` with one row
per community and the columns `gamma`, `mu`, `sympatry` (fixed or
estimated), `loglik`, `richness` and `convergence` (0 when
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) converged). With
`shared = TRUE` every row holds the same rates and convergence code,
`loglik` is the contribution of each community at those rates, and the
attribute `loglik_total` holds their sum.

## Details

At the root, `"both"` averages the likelihoods given an absent and a
present ancestor, as the default of the DAMOCLES package
(`pchoice = 0`); the article sums them, which only adds `log(2)` to the
log-likelihood. The optimisation works on the logarithm of the rates,
with [`stats::optim()`](https://rdrr.io/r/stats/optim.html) from several
starting points, so estimates of zero appear as very small rates. Under
this model the tip states hold information on `gamma` and `mu` only
through their pattern on the tree, so with few species, or with a
community holding nearly all or nearly none of them, the likelihood
surface is flat and the estimates are poor. In
`data-raw/validation_damocles.R`, ten communities simulated with
`gamma = mu = 1` on trees of 60 tips gave estimates of `mu` from zero to
several thousand, with the same likelihoods as
`DAMOCLES::DAMOCLES_ML()`; the ratio of the two rates is better
estimated than each one.

With `shared = TRUE`, `fit_damocles()` estimates one pair of rates for
all communities by maximising the sum of their log-likelihoods, which
treats the communities as independent realisations of the same process
on the same tree. Pooling communities narrows the estimates: in
`data-raw/validation_damocles.R`, with `gamma = mu = 1` on trees of 60
tips and height one (30 fits each), 90% of the estimates of `gamma` fell
between 0.49 and 6.0 from one community and between 0.75 and 1.27 from
twenty, and those of `mu` between zero and 18.9 from one and between
zero and 1.8 from twenty. `mu` stays the harder rate, and some fits
still put it near zero. The assumption of independence matters:
neighbouring sites that share their colonisation history, or communities
nested in one another, are not independent replicates and pooling them
overstates the information on the rates. Twice the difference between
the summed log-likelihoods of separate and shared fits can be compared
with a chi-squared distribution on \\2 (m - 1)\\ degrees of freedom,
with \\m\\ communities, to ask whether the rates differ among
communities.

With `sympatry = NULL`, `fit_damocles()` also estimates the probability
that a present ancestor leaves both daughters present, on the logit
scale. A single community says almost nothing about it: in
`data-raw/validation_damocles.R`, with `sympatry = 0.5` and
`gamma = mu = 1` on trees of 60 tips, 90% of the estimates fell between
0 and 1 from one community, between 0.26 and 0.81 from five and between
0.41 and 0.60 from twenty. Estimate it with shared rates, and compare
the fit with that of `sympatry = 0` by a likelihood ratio (one degree of
freedom, on the boundary of the parameter space, so the chi-squared
p-value is conservative) to ask whether strictly allopatric speciation
is enough.

Rates can change across the tree. `damocles_loglik()` and
[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)
take one rate per edge, and `fit_damocles()` estimates rates that differ
among `regimes` of the edges (a partition defined a priori with
[`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md),
or the configuration of rate shifts of a method such as BAMM) or that
change with an `edge_covariate` (such as branch-specific rates of
speciation from a method such as ClaDS), for colonisation, extinction or
both (`vary`). Either input is a plain vector in the order of
`tree$edge`, so any method can supply it. Rates of diversification
describe speciation and extinction of the whole lineage, not local
colonisation and extinction: the regimes or the covariate are a
hypothesis about the local dynamics, tested by comparing the fit with
the constant model by
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md)
(a single regime, or `beta = 0`, gives the constant model) and by
intervals from
[`confint.mc_damocles_fit()`](https://provetelab.org/comdynamite/reference/confint.mc_damocles_fit.md).
Each regime adds parameters to a likelihood that is already flat, so
these fits need several communities with `shared = TRUE`; the type I
error and power of that comparison are in
`data-raw/validation_damocles.R`: with rates shared by 5 or 20
communities on trees of 60 tips, the test rejected constant rates in
4–7% of 100 replicates without an effect and in 86–99% with a threefold
difference between clades or a slope of one on an autocorrelated
covariate. Those communities were simulated independently; real
communities that share their history will reject more often.

## References

Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for
phylogenetic community structure. *Ecology Letters*, 18, 153–163.
[doi:10.1111/ele.12395](https://doi.org/10.1111/ele.12395)

## See also

[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md),
[`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md).

## Examples

``` r
set.seed(2)
tr <- treeshapesim::sim_tree_shape(60)
comm <- simulate_damocles(tr, gamma = 1, mu = 2, n = 2)
damocles_loglik(comm, tr, gamma = 1, mu = 2)
#>      sim1      sim2 
#> -28.13747 -22.54596 
fit_damocles(comm, tr)
#> DAMOCLES rates by maximum likelihood (Pigot & Etienne 2015)
#>      gamma       mu sympatry loglik richness convergence
#> sim1 1.794 4.61e+00        0  -27.8       11           0
#> sim2 0.476 1.46e-09        0  -21.9        8           0
fit_damocles(comm, tr, shared = TRUE)
#> DAMOCLES rates by maximum likelihood (Pigot & Etienne 2015)
#> Parameters shared by 2 communities; total log-likelihood -50.6593 
#>      gamma   mu sympatry loglik richness convergence
#> sim1  1.12 2.69        0  -28.0       11           0
#> sim2  1.12 2.69        0  -22.6        8           0
```
