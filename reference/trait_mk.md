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
trait_mk(tree, q = 1, states = c("0", "1"), root = NULL)
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

## Value

The tree, reordered "cladewise", with class `c("simmap", "phylo")` and,
as in
[`phytools::sim.history()`](https://rdrr.io/pkg/phytools/man/sim.history.html):
`maps`, a list with one element per edge giving the time spent in each
state in the order they occurred; `mapped.edge`, the total time in each
state per edge; `node.states`, the states at the start and end of each
edge; and `states`, the states of the tips. `Q` holds the rate matrix
used.

## References

Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from
discrete morphological character data. *Systematic Biology*, 50,
913–925.
[doi:10.1080/106351501753462876](https://doi.org/10.1080/106351501753462876)

## See also

[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md)

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
