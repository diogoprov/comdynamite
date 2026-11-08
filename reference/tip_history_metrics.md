# Tip-based metrics of a known trait history

Computes, for each tip, the metrics of trait evolution of Luza et al.
(2021) from a history of a discrete trait, as
`Herodotools::calc_tip_based_trait_evo()` computes them from a
stochastic character map. Applied to a history from
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
they are the true values that a reconstruction from the tips should
recover. Two metrics in continuous time are added, which count the
changes inside branches too.

## Usage

``` r
tip_history_metrics(history)
```

## Arguments

- history:

  A `simmap` tree, from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  or phytools.

## Value

Data frame with one row per tip and the columns `transitions`,
`total_nodes`, `prop_transitions`, `last_transition_time`,
`stasis_time`, `changes` and `time_since_change`.

## Details

Following the code of `calc_tip_based_trait_evo()`, the state of a node
is the state at the start of the edges that leave it, and the path of a
tip runs from its parent node to the root. `transitions` counts the
changes of state between consecutive nodes of the path, starting from
the state of the tip; `prop_transitions` divides it by the number of
nodes in the path, `total_nodes`. `last_transition_time` sums the length
of the terminal edge and of the edges leading to the nodes of the path,
from the tip towards the root, while their state is the state of the
tip; the root, which has no edge leading to it, adds nothing.
`stasis_time` is the longest time spent in one state within a single
edge, over the terminal edge and the edges leading to the nodes of the
path. Only nodes enter `transitions` and `last_transition_time`, so
changes that revert within a branch are missed; `changes` counts every
change along the path from the root to the tip and `time_since_change`
is the time since the last of them, or the age of the tip from the root
when there was none.

## References

Luza, A. L., Maestri, R., Debastiani, V. J., Patterson, B. D., Hartz, S.
M., & Duarte, L. D. S. (2021). Is evolution faster at ecotones? A test
using rates and tempo of diet transitions in Neotropical Sigmodontinae
(Rodentia, Cricetidae). *Ecology and Evolution*, 11, 18676–18690.
[doi:10.1002/ece3.8476](https://doi.org/10.1002/ece3.8476)

## See also

[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)

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
```
