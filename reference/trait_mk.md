# A discrete trait evolving along a tree, with its whole history

Simulates a discrete character under the Mk model (Lewis 2001): along
each branch the state changes after exponential waiting times, with the
rates of `q`. The result keeps every change, in the format of a
stochastic character map of phytools (class `simmap`), so it is the
known history against which reconstructions can be compared, for example
by
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md).

## Usage

``` r
trait_mk(
  tree,
  q = 1,
  states = c("0", "1"),
  root = NULL,
  alpha = NULL,
  edge_rates = NULL,
  rate_classes = NULL,
  class_switch = 0.1
)
```

## Arguments

- tree:

  A rooted `phylo` tree with branch lengths.

- q:

  One rate for every change between states (the equal-rates model), or a
  square matrix of rates from the state of each row to the state of each
  column, with the states as names; its diagonal is ignored.

- states:

  Names of the states, when `q` is one rate.

- root:

  State at the root; `NULL` draws one with equal probabilities.

- alpha:

  Shape (and rate) of the gamma distribution of the rate multiplier of
  each branch; `NULL` for none. Small values make rates more uneven
  among branches.

- edge_rates:

  Optional multipliers of the rates of each branch, in the order of the
  rows of `tree$edge`, instead of drawing them with `alpha`.

- rate_classes:

  Optional multipliers of the rates in each hidden class, such as
  `c(slow = 0.1, fast = 1)`; classes without names are called `R1`, `R2`
  and so on.

- class_switch:

  Rate of moving from one hidden class to each other one.

## Value

The tree, reordered "cladewise", with class `c("simmap", "phylo")` and,
as in
[`phytools::sim.history()`](https://rdrr.io/pkg/phytools/man/sim.history.html):
`maps`, a list with one element per edge giving the time spent in each
state in the order they occurred; `mapped.edge`, the total time in each
state per edge; `node.states`, the states at the start and end of each
edge; and `states`, the states of the tips. `Q` holds the rate matrix
used. With rates varying among branches, `edge_rates` holds the
multiplier of each edge, in the new order of the edges. With
`rate_classes`, `hidden` holds the same history with the class of each
segment, as states named `"state|class"`, with its own `maps`,
`node.states`, `states` and `Q`, and `classes` holds the class of each
tip.

## Details

Rates can vary among branches in two ways. With `alpha` (or
`edge_rates`), all rates of a branch are multiplied by one factor, drawn
for each branch from a gamma distribution with shape and rate `alpha`
(mean 1), as in the model that
[`phytools::fitgammaMk()`](https://rdrr.io/pkg/phytools/man/fitMk.html)
fits (Revell 2024). With `rate_classes`, the character has hidden rate
classes, as in the hidden-rates models of corHMM (Beaulieu et al. 2013):
in class k, changes between states happen at the rates of `q` times
`rate_classes[k]`, and the lineage moves between classes at rate
`class_switch` without changing its state. The class at the root is
drawn with equal probabilities.

## References

Beaulieu, J. M., O'Meara, B. C., & Donoghue, M. J. (2013). Identifying
hidden rate changes in the evolution of a binary morphological
character: the evolution of plant habit in campanulid angiosperms.
*Systematic Biology*, 62, 725–737.
[doi:10.1093/sysbio/syt034](https://doi.org/10.1093/sysbio/syt034)

Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from
discrete morphological character data. *Systematic Biology*, 50,
913–925.
[doi:10.1080/106351501753462876](https://doi.org/10.1080/106351501753462876)

Revell, L. J. (2024). Testing for heterogeneous rates of discrete
character evolution on phylogenies. *Journal of Evolutionary Biology*,
37, 1591–1602.
[doi:10.1093/jeb/voae119](https://doi.org/10.1093/jeb/voae119)

## See also

[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md)
for continuous traits.

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(10)
h <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
h$states
#>       t1       t2       t3       t4       t5       t6       t7       t8 
#>   "open"   "open" "forest"   "open"   "open" "forest" "forest"   "open" 
#>       t9      t10 
#>   "open"   "open" 
tip_history_metrics(h)
#>     transitions total_nodes prop_transitions last_transition_time stasis_time
#> t1            0           3             0.00            1.0000000   0.7085115
#> t2            0           4             0.00            1.0000000   0.7085115
#> t3            1           4             0.25            0.1435450   0.7085115
#> t4            0           2             0.00            1.0000000   0.7085115
#> t5            0           3             0.00            1.0000000   0.4696704
#> t6            1           4             0.25            0.3535646   0.4696704
#> t7            1           4             0.25            0.3535646   0.4696704
#> t8            0           3             0.00            1.0000000   0.4696704
#> t9            0           4             0.00            1.0000000   0.4696704
#> t10           0           4             0.00            1.0000000   0.4696704
#>     changes time_since_change
#> t1        0        1.00000000
#> t2        0        1.00000000
#> t3        1        0.02466665
#> t4        0        1.00000000
#> t5        0        1.00000000
#> t6        1        0.27902759
#> t7        1        0.27902759
#> t8        0        1.00000000
#> t9        0        1.00000000
#> t10       0        1.00000000
# rates that vary among branches
g <- trait_mk(tr, q = 0.5, alpha = 0.3)
round(g$edge_rates, 2)
#>  [1] 0.14 0.40 0.22 1.87 1.52 0.00 2.22 0.88 0.15 0.03 0.40 1.67 0.29 0.00 0.00
#> [16] 0.46 0.50 0.03
hr <- trait_mk(tr, q = 0.5, rate_classes = c(slow = 0.1, fast = 2), class_switch = 0.2)
hr$classes
#>     t1     t2     t3     t4     t5     t6     t7     t8     t9    t10 
#> "slow" "slow" "slow" "slow" "slow" "slow" "slow" "slow" "slow" "fast" 
```
