# History of the geographic ranges of lineages, as a stochastic map

Maps, along the tree of
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
the range of each lineage at each moment, as the combination of areas in
which it occupies at least one site: `"A"`, `"B"`, `"A+B"` and so on,
with the areas in the order of their levels. The result is a tree of
class `simmap`, as from
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
or phytools, with the true history, for
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md)
and
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md),
and for comparison with ranges reconstructed by models such as DEC
(BioGeoBEARS).

## Usage

``` r
range_history(x, areas, tree = c("reconstructed", "full"))
```

## Arguments

- x:

  An `mc_evolving` object from
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md).

- areas:

  Area of each site, in the order of the sites or named by site number;
  a factor gives the order of the areas in the names of the ranges. Area
  names may not contain `"+"`.

- tree:

  `"reconstructed"` (the default) for the tree of extant species,
  `x$tree`, or `"full"` for the tree of all lineages, extinct ones
  included, `x$tree_full`.

## Value

A tree of class `c("simmap", "phylo")` with `maps` (for each edge, the
time spent in each range, from its start; a change in the last
generation of a tip is kept as a piece without duration, so that the tip
has the final range), `mapped.edge` (edges x ranges, the total time) and
`root_range`, or `NULL` when the tree is `NULL`.

## Details

The range of a lineage changes at the end of each generation in which it
colonises the first site, or loses the last site, of an area. A new
lineage starts with the sites it took from its mother. Branches of the
tree follow lineages: a branch that continues a mother after a budding
is the mother, and the stem of a new lineage is the new lineage. In the
reconstructed tree a branch may join several lineages of the full tree
(an extinct mother and the descendant that survived), and its map joins
their ranges in order. Each tip of the reconstructed tree has the range
of the lineage that represents its species (with protracted speciation,
not the union of the incipient lineages of the species).

A node has the range of the mother just before the split (the union of
the ranges that the mother and the new lineage start with, apart from
changes in the same generation), and the tree keeps the range at the
root, which no edge shows, as `root_range`. The stem from the start of
the simulation to the first split (`root.edge`) is not mapped. Nodes of
a starting tree, before the simulation, take the union of the ranges at
the start of the species below them, and the parts of edges before the
start take the range at the start.

## See also

[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md)

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md),
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

## Examples

``` r
set.seed(1)
ls <- landscape_grid(5, disp_scale = 0.5)
areas <- ifelse(ls$coords$x < 2.5, "W", "E")
ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
sim <- simulate_evolving_metacommunity(ls, ancestor, dispersal = 0.01, extirpation = 0.02,
                                       n_steps = 200,
                                       speciation = speciation_control("fission", rate = 0.02))
h <- range_history(sim, areas)
colSums(h$mapped.edge)
#>   W W+E 
#> 132 698 
```
