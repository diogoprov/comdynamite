# Test a partition of the tree against random partitions

`damocles_partition_test()` asks whether rates that change among
`regimes` fit the local presences better than rates that change among
partitions of the same complexity placed at random on the tree
([`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md)).
The statistic is the likelihood ratio of the fit with `regimes` against
constant rates, with parameters shared by all communities
([`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md));
its p-value is the proportion of random partitions with a statistic as
large or larger.

## Usage

``` r
damocles_partition_test(
  comm,
  tree,
  regimes,
  n_rand = 99,
  vary = c("both", "gamma", "mu"),
  root = c("both", "absent", "present"),
  sympatry = 0,
  tolerance = 2,
  max_overlap = 0.25,
  cores = 1
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

- regimes:

  Regime of each edge of `tree`, in the order of `tree$edge`, with two
  or more regimes.

- n_rand:

  Number of random partitions.

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

- tolerance, max_overlap:

  Passed to
  [`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md):
  random clades hold between `1 / tolerance` and `tolerance` times as
  many edges as the clades of `regimes` and share at most `max_overlap`
  of their edges with them (Jaccard index).

- cores:

  Number of cores, with random partitions fitted in parallel. Not
  available on Windows, where one core is used. For reproducible
  p-values with more than one core, call `RNGkind("L'Ecuyer-CMRG")`
  before [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

A list of class `mc_damocles_partition` with the observed `statistic`,
its degrees of freedom `df`, the chi-squared p-value `p_chisq`, the
p-value against random partitions `p_random`, the statistics of the
random partitions `null`, and `fit`, the fit with `regimes`.

## Details

The chi-squared test of
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md)
asks whether rates vary at all along the tree in the way the partition
allows. When they vary somewhere, any partition captures part of that
and the chi-squared test rejects constant rates even for a partition
unrelated to where the rates change. The random partitions keep the
number of parameters and the sizes of the clades, so this test asks
instead whether the partition given is better placed than arbitrary
ones: whether the regimes from a method of diversification, or a
partition defined a priori, carry information on local colonisation and
extinction beyond rate heterogeneity in general.
`data-raw/validation_damocles.R` shows the difference between the two
tests: with 40 replicates of 10 communities on trees of 60 tips and a
clade holding about a quarter of the edges, the chi-squared test
rejected constant rates in 3 of 40 without an effect, in 40 of 40 with
colonisation faster in the clade, and in 29 of 40 with it faster in
another clade of similar size; the test against random partitions
rejected in 4, 37 and 0 of 40. The 4 of 40 without an effect is within
sampling error of 5% but may be slightly liberal.

Random partitions that nearly coincide with the one tested, such as the
clade one edge above it, capture the same rates and leave the test
without power: with `max_overlap = 1` the test against random partitions
rejected in none of 40 replicates with the effect in the clade. The
default `max_overlap = 0.25` keeps only placements elsewhere on the
tree, so the test asks whether the partition beats partitions of the
same complexity placed somewhere else. With the whales and the best BAMM
configuration of BAMMtools (a regime on 69 of 172 edges) and
colonisation three times faster in it, both tests rejected constant
rates (p = 0.007 and 0.01; `data-raw/damocles_bamm.R`).

## See also

[`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md),
[`damocles_samples()`](https://provetelab.org/comdynamite/reference/damocles_samples.md),
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md).

## Examples

``` r
set.seed(9)
tr <- treeshapesim::sim_tree_shape(30)
reg <- clade_regimes(tr, list(A = ape::Ntip(tr) + 2))
comm <- simulate_damocles(tr, ifelse(reg == "A", 3, 0.5), mu = 1, n = 5)
damocles_partition_test(comm, tr, reg, n_rand = 19, vary = "gamma")
#> Partition of the tree against 19 random partitions (DAMOCLES)
#> Likelihood ratio against constant rates: 28.07 on 1 df
#> p (chi-squared): 1.17e-07 
#> p (random partitions): 0.05 
```
