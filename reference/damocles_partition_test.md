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
  n_boot = 0,
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

  Number of random partitions drawn; each distinct one is fitted once.
  One whole number, 0 or more; 0 skips the test against random
  partitions.

- n_boot:

  Number of parametric bootstrap simulations under the constant rates
  estimated. One whole number, 0 or more; 0, the default, skips the
  bootstrap.

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
  many edges as the clades of `regimes`, share at most `max_overlap` of
  their edges with them (Jaccard index), and the whole random partition
  agrees with `regimes` at most as much (adjusted Rand index).
  `tolerance` is one number, 1 or more; `max_overlap` one number between
  0 and 1.

- cores:

  Number of cores, with random partitions and bootstrap simulations
  fitted in parallel. One whole number, 1 or more. Not available on
  Windows, where one core is used. The random partitions and the
  simulated communities are drawn before the fits, so the p-values do
  not depend on the number of cores.

## Value

A list of class `mc_damocles_partition` with the observed `statistic`,
its degrees of freedom `df`, the chi-squared p-value `p_chisq`, the
p-value against random partitions `p_random` (`NA` when no random
partition could be drawn), the statistics of the distinct random
partitions `null`, their number `n_distinct`, `n_rand`, the p-value of
the parametric bootstrap `p_boot` (`NA` without it), the simulated
statistics `null_boot`, and `fit`, the fit with `regimes`.

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

The test needs room on the tree for clades like those of `regimes`
elsewhere. When the partition tested is one of the clades at the root,
or a clade holding most of the tree, the only placements allowed may
split the tree as `regimes` does (its complement) or be few, and the
random statistics then repeat a handful of values: the p-value falls at
its smallest or at one depending on numerical noise. In simulations with
one shift in the rate of speciation and no effect on local rates
(`data-raw/s1_type1_shifts.R`), this made the test reject in 12% to 40%
of replicates before random partitions matching `regimes` were excluded
([`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md))
and before the p-value counted each distinct partition once. The
`n_rand` draws are a way to find distinct placements: the p-value is the
proportion of the distinct random partitions, plus the one tested, with
a statistic as large or larger, so with fewer than 19 distinct ones it
cannot fall below 0.05 (a warning says so). Trees with a strong shift in
diversification often leave that little room. In the same simulations
after both changes (50 replicates per setting, 99 draws), replicates
with one shift had a median of 3 distinct random partitions and 79 of
500 had none; replicates with 19 or more rejected in 8% (37 of 452).

`n_boot` adds a parametric bootstrap of the same likelihood ratio:
communities are simulated with the constant rates estimated from `comm`
([`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)),
as many as in `comm`, and each set is fitted with and without `regimes`;
`p_boot` is the proportion of simulated statistics, plus the observed
one, as large or larger. It asks the question of the chi-squared test
(do rates vary along the tree as the partition allows?) without relying
on the asymptotic distribution of the likelihood ratio, which is
doubtful with few communities and rates near zero, and it needs no room
on the tree, so it is available for any partition. It does not ask
whether the partition is better placed than others: with rate
heterogeneity elsewhere on the tree it rejects, like the chi-squared
test. Each simulation costs two fits, twice a random partition. Set
`n_rand = 0` to run the bootstrap alone.

## See also

[`random_regimes()`](https://provetelab.org/comdynamite/reference/random_regimes.md),
[`damocles_samples()`](https://provetelab.org/comdynamite/reference/damocles_samples.md),
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md).

## Examples

``` r
set.seed(9)
tr <- treeshapesim::sim_tree_shape(30)
reg <- clade_regimes(tr, list(A = ape::Ntip(tr) + 2))
comm <- simulate_damocles(tr, ifelse(reg == "A", 3, 0.5), mu = 1, n = 3)
# few sites and random partitions, to keep the example quick
damocles_partition_test(comm, tr, reg, n_rand = 2, vary = "gamma")
#> Warning: only 2 distinct random partitions in 2 draws: the tree has few places for clades like those of `regimes`, and `p_random` cannot be below 0.33
#> Partition of the tree against constant rates (DAMOCLES)
#> Likelihood ratio against constant rates: 2.857 on 1 df
#> p (chi-squared): 0.091 
#> p (random partitions): 0.333 (2 distinct random partitions)
# \donttest{
# random partitions and the parametric bootstrap refit the model many times
comm <- simulate_damocles(tr, ifelse(reg == "A", 3, 0.5), mu = 1, n = 5)
damocles_partition_test(comm, tr, reg, n_rand = 19, vary = "gamma")
#> Warning: only 18 distinct random partitions in 19 draws: the tree has few places for clades like those of `regimes`, and `p_random` cannot be below 0.053
#> Partition of the tree against constant rates (DAMOCLES)
#> Likelihood ratio against constant rates: 24.81 on 1 df
#> p (chi-squared): 6.33e-07 
#> p (random partitions): 0.0526 (18 distinct random partitions)
damocles_partition_test(comm, tr, reg, n_rand = 0, n_boot = 19, vary = "gamma")
#> Partition of the tree against constant rates (DAMOCLES)
#> Likelihood ratio against constant rates: 24.81 on 1 df
#> p (chi-squared): 6.33e-07 
#> p (parametric bootstrap): 0.05 (19 simulations under constant rates)
# }
```
