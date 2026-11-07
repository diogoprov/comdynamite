# Check whether a simulation has settled

Compares the first and the second half of the recorded steps. With
phylogenetic competition and high dispersal, regional exclusion can go
on for thousands of steps, and the last step of such a run is not an
equilibrium. A small change does not prove equilibrium either, since
slow exclusion can fall below `tolerance` over a short record: record
steps over a long stretch of the run.

## Usage

``` r
equilibrium_check(x, tolerance = 0.05)
```

## Arguments

- x:

  An `mc_sim` object with `history`, from `record_steps` of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  with at least four recorded steps.

- tolerance:

  Largest relative change in a statistic that counts as settled.

## Value

Data frame with one row for `gamma`, the species in all sites, `alpha`,
the mean species per site, and `abundance`, the mean individuals per
site, with their means in the two halves (`early`, `late`), the relative
`change` and whether it is `settled`, invisibly. It gives a warning if
any statistic has not settled.

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(15)
sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree),
                              competition_matrix(tree, "phylo"), n_steps = 200,
                              seed_steps = 1:5, record_steps = seq(100, 200, by = 10))
suppressWarnings(equilibrium_check(sim))
```
