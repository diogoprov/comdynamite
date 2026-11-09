# Rates of trait evolution from the turnover of partners in a metacommunity

Turns the turnover of interaction partners of each species in a
simulated metacommunity, from
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
into rates of evolution of a continuous trait on the branches of its
tree, for `trait_rates(edge_rates = )`. It implements the expectation
that evolution speeds up with each change in the community around a
lineage and stays slow while the community is stable (McPeek 2017;
scenario 5 of the notes on branch rates).

## Usage

``` r
turnover_rates(
  x,
  tree = x$tree,
  sigma2 = 1,
  k = 1,
  window = NULL,
  partners = c("all", "competitors")
)
```

## Arguments

- x:

  An `mc_sim` object with `history`, that is, run with `record_steps`.

- tree:

  The tree of the species in `x`; by default the one stored in `x`.

- sigma2:

  Rate outside the stretch informed by the metacommunity.

- k:

  Increase of the rate, relative to `sigma2`, per unit of turnover
  (partners gained or lost per site and time step).

- window:

  Length of the end of each tip branch represented by the metacommunity,
  in the units of the branch lengths; `NULL` uses the whole tip branch.
  A window longer than a tip branch is cut to its length.

- partners:

  `"all"` counts every other species as a partner; `"competitors"`
  counts only the species with a non-zero competition coefficient with
  the focal species in `x$alpha` (in either direction).

## Value

Numeric vector with the rate of each branch, in the order of the rows of
`tree$edge`, with the turnover of each species as attribute
`"turnover"`.

## Details

The metacommunity runs in ecological time among the species at the tips,
so its dynamics only inform the tip branches, and how ecological time
relates to branch lengths is a choice of the user, not something the
simulation gives. Here the dynamics are taken to represent the last
`window` units of time of each tip branch: in that stretch the rate is
\\\mathrm{sigma2}\\ (1 + k \cdot \mathrm{turnover})\\, and elsewhere,
including all internal branches, it is `sigma2`. The rate of a tip
branch is the mean over the branch, weighted by time.

## References

McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton
University Press.
[doi:10.23943/princeton/9780691088778.001.0001](https://doi.org/10.23943/princeton/9780691088778.001.0001)

## See also

[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md)

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(10)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tr),
                              competition_matrix(tr, "none"), n_steps = 300,
                              record_steps = seq(100, 300, 20), extirpation = 0.02)
rates <- turnover_rates(sim, k = 20, window = 0.1)
x <- trait_rates(tr, edge_rates = rates)
```
