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
and combined at each node by the rule of allopatric speciation.
`fit_damocles()` finds the maximum-likelihood rates.

## Usage

``` r
damocles_loglik(comm, tree, gamma, mu, root = c("both", "absent", "present"))

fit_damocles(comm, tree, root = c("both", "absent", "present"), start = NULL)
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
  lengths.

- mu:

  Rate of local extinction per lineage.

- root:

  Prior on the state of the ancestor at the root: `"both"`, `"absent"`
  or `"present"`.

- start:

  Optional starting rates, `c(gamma = , mu = )`; by default the
  optimisation starts from four combinations of 0.2 and 2 divided by the
  height of the tree.

## Value

`damocles_loglik()`: the log-likelihood of each community.
`fit_damocles()`: a data frame of class `mc_damocles_fit` with one row
per community and the columns `gamma`, `mu`, `loglik`, `richness` and
`convergence` (0 when
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) converged).

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
#>      gamma       mu loglik richness convergence
#> sim1 1.794 4.61e+00  -27.8       11           0
#> sim2 0.476 1.46e-09  -21.9        8           0
```
