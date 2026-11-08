# Local presence along a tree: the dynamic null model of DAMOCLES

`simulate_damocles()` simulates the presence and absence of each lineage
of a clade in one local community, as it diversifies along `tree`, under
the dynamic null model of community assembly of Pigot & Etienne (2015),
DAMOCLES. Along each branch the lineage is either locally absent (state
0) or present (state 1); an absent lineage colonises the community at
rate `gamma` and a present one goes locally extinct at rate `mu`, the
same for every lineage. At each speciation, which is allopatric, a
locally present ancestor leaves one daughter present and the other
absent, chosen at random, and an absent ancestor leaves two absent
daughters. The states at the tips form the community.

## Usage

``` r
simulate_damocles(
  tree,
  gamma,
  mu,
  n = 1,
  root = c("both", "absent", "present"),
  sympatry = 0
)
```

## Arguments

- tree:

  A rooted, fully dichotomous `phylo` tree with branch lengths.

- gamma:

  Rate of local colonisation per lineage, in the units of the branch
  lengths.

- mu:

  Rate of local extinction per lineage.

- n:

  Number of communities to simulate.

- root:

  State of the ancestor at the root, just before the first speciation
  (or at the start of the root edge): `"both"` draws absent or present
  with equal probabilities for each community, as in the simulations of
  Pigot & Etienne (2015); `"absent"` or `"present"` fix it.

- sympatry:

  Probability, between 0 and 1, that a locally present ancestor leaves
  both daughters present at speciation; 0 (the default) is strictly
  allopatric speciation.

## Value

A matrix of zeros and ones with one row per community and one column per
tip, named by the tips, ready to be used as `comm` in the other
functions of the package.

## Details

`sympatry` relaxes the rule of strictly allopatric speciation: with that
probability a present ancestor leaves both daughters present, as after
speciation within the community (in situ or sympatric speciation) or
when one daughter recolonises it before the two can be told apart;
otherwise one daughter is present, as above. `sympatry = 0`, the
default, is the model of Pigot & Etienne (2015). An absent ancestor
always leaves two absent daughters.

The states at the nodes are drawn from the exact transition
probabilities of the two-state process over each branch, which gives the
same distribution of tip states as simulating each colonisation and
extinction event. The model assumes the tree holds the species of the
regional pool, so it is meant for ultrametric trees of extant species.
Branches of length zero and a root edge (`tree$root.edge`) are allowed.

## References

Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for
phylogenetic community structure. *Ecology Letters*, 18, 153–163.
[doi:10.1111/ele.12395](https://doi.org/10.1111/ele.12395)

## See also

[`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
and
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
for the likelihood,
[`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md)
for the null model of phylogenetic community structure.

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(30)
comm <- simulate_damocles(tr, gamma = 0.5, mu = 0.5, n = 5)
rowSums(comm)
#> sim1 sim2 sim3 sim4 sim5 
#>    4    4    2    6    3 
```
