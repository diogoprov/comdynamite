# A continuous trait whose rate of evolution varies among branches

Simulates Brownian motion in which each branch has its own rate
(variance per unit time), and keeps those rates, so that methods which
estimate rates per branch can be compared with the rates that generated
the data. Three models of how rates vary are available:

## Usage

``` r
trait_rates(
  tree,
  model = c("shift", "gbm", "iid"),
  sigma2 = 1,
  shifts = 1,
  shift_size = 4,
  shift_edges = NULL,
  shift_branch = FALSE,
  rate_sd = 0.5,
  root_value = 0
)
```

## Arguments

- tree:

  A rooted `phylo` tree with branch lengths.

- model:

  How rates vary among branches: `"shift"`, `"gbm"` or `"iid"`.

- sigma2:

  Rate at the root (`"shift"`, `"gbm"`) or median rate (`"iid"`).

- shifts:

  Number of shifts, for `"shift"` when `shift_edges` is not given.

- shift_size:

  Multipliers of the rate at each shift, recycled to the number of
  shifts. Values below 1 slow evolution down.

- shift_edges:

  Optional row numbers of `tree$edge` where the shifts happen.

- shift_branch:

  Should a shift change only its own branch, rather than all branches
  descending from it?

- rate_sd:

  Standard deviation of the logarithm of the rate: per square root of
  time for `"gbm"` and among branches for `"iid"`.

- root_value:

  Trait value at the root.

## Value

A list with

- `tips`: the trait at the tips, named by species;

- `nodes`: the trait at the internal nodes, named by their numbers in
  `tree$edge`, the root first;

- `edge_rate`: the rate of each branch, in the order of the rows of
  `tree$edge`;

- `node_rate`: for `"gbm"`, the rate at every node, tips first, named by
  node number; `NULL` otherwise;

- `shifts`: for `"shift"`, a data frame with the `edge` (row of
  `tree$edge`), the `node` at its end and the `multiplier`; `NULL`
  otherwise;

- `tree`: the tree given, and `scaled`: the tree with each branch length
  multiplied by its rate, under which the trait is ordinary Brownian
  motion.

## Details

- `"shift"`: the rate is `sigma2` at the root and is multiplied by
  `shift_size` at a few branches. With `shift_branch = FALSE` the new
  rate is inherited by every branch descending from the shifted one, as
  in the models of Eastman et al. (2011), traitMedusa in motmot (Thomas
  & Freckleton 2012) or bifrost; with `shift_branch = TRUE` only the
  shifted branch changes. The branches are drawn with probabilities
  proportional to their length, or given in `shift_edges`.

- `"gbm"`: the logarithm of the rate evolves itself by Brownian motion
  with standard deviation `rate_sd` per square root of time, starting
  from `log(sigma2)` at the root, as in the variable-rate models of
  Revell (2021) and Martin et al. (2023). The rate of a branch is the
  geometric mean of the rates at its two ends.

- `"iid"`: the rate of each branch is `sigma2` times an independent
  lognormal multiplier with median 1 and `rate_sd` as the standard
  deviation of its logarithm.

The trait is then Brownian motion on the tree whose branch lengths are
multiplied by their rates, drawn at the tips and then at the internal
nodes given the tips, as in
[`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md).

## References

Eastman, J. M., Alfaro, M. E., Joyce, P., Hipp, A. L., & Harmon, L. J.
(2011). A novel comparative method for identifying shifts in the rate of
character evolution on trees. *Evolution*, 65, 3578–3589.
[doi:10.1111/j.1558-5646.2011.01401.x](https://doi.org/10.1111/j.1558-5646.2011.01401.x)

Martin, B. S., Bradburd, G. S., Harmon, L. J., & Weber, M. G. (2023).
Modeling the evolution of rates of continuous trait evolution.
*Systematic Biology*, 72, 590–605.
[doi:10.1093/sysbio/syac068](https://doi.org/10.1093/sysbio/syac068)

Revell, L. J. (2021). A variable-rate quantitative trait evolution model
using penalized-likelihood. *PeerJ*, 9, e11997.
[doi:10.7717/peerj.11997](https://doi.org/10.7717/peerj.11997)

Thomas, G. H., & Freckleton, R. P. (2012). MOTMOT: models of trait
macroevolution on trees. *Methods in Ecology and Evolution*, 3, 145–151.
[doi:10.1111/j.2041-210X.2011.00132.x](https://doi.org/10.1111/j.2041-210X.2011.00132.x)

## See also

[`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md)
for a single rate,
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
for discrete traits with rates that vary among branches.

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(20)
x <- trait_rates(tr, "shift", shifts = 1, shift_size = 8)
x$shifts
#>   edge node multiplier
#> 1   38   20          8
head(x$tips)
#>         t1         t2         t3         t4         t5         t6 
#>  0.3876716  0.3725845  0.1821453 -0.3773427 -0.5333873 -0.5354236 
g <- trait_rates(tr, "gbm", rate_sd = 0.5)
summary(g$edge_rate)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5037  0.7585  0.8401  0.8745  0.9986  1.1787 
# as the environmental optima of a metacommunity
traits <- species_traits(tr, optima = g$tips)
```
