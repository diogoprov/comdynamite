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
  shift_rate = NULL,
  shift_branch = FALSE,
  relative = TRUE,
  rate_sd = 0.5,
  root = 0,
  edge_rates = NULL,
  root_value = NULL
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

  Number of shifts, for `"shift"` when neither `shift_edges` nor
  `shift_rate` is given.

- shift_size:

  Multipliers of the rate at each shift, recycled to the number of
  shifts, or a function of the number of shifts that gives them. Values
  below 1 slow evolution down.

- shift_edges:

  Optional row numbers of `tree$edge` at whose start the shifts happen.

- shift_rate:

  Optional rate of shifts per unit of branch length, for shifts as a
  Poisson process along the branches.

- shift_branch:

  Should a shift last only to the end of its branch, rather than pass to
  all branches descending from it?

- relative:

  Multiply the rate before each shift by `shift_size` (`TRUE`), or set
  the rate after it to `sigma2` times `shift_size` (`FALSE`)?

- rate_sd:

  Standard deviation of the logarithm of the rate: per square root of
  time for `"gbm"` and among branches for `"iid"`.

- root:

  Trait value at the root.

- edge_rates:

  Optional rate of each branch, in the order of the rows of `tree$edge`,
  such as those of
  [`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md).
  When given, `model` and its arguments are ignored.

- root_value:

  Deprecated; use `root`.

## Value

A list with

- `tips`: the trait at the tips, named by species;

- `nodes`: the trait at the internal nodes, named by their numbers in
  `tree$edge`, the root first;

- `edge_rates`: the rate of each branch, in the order of the rows of
  `tree$edge`, as in
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  (also as `edge_rate`, its name before 0.7.0, until 1.0.0);

- `node_rate`: for `"gbm"`, the rate at every node, tips first, named by
  node number; `NULL` otherwise;

- `shifts`: for `"shift"`, a data frame with the `edge` (row of
  `tree$edge`), the `node` at its end, the `position` of the shift as
  time from the start of the edge, the `multiplier` and the `rate` after
  the shift; `NULL` otherwise;

- `tree`: the tree given, and `scaled`: the tree with each branch length
  multiplied by its rate, under which the trait is ordinary Brownian
  motion.

## Details

- `"shift"`: the rate is `sigma2` at the root and changes at a few
  points of the tree. With `shift_branch = FALSE` the new rate is
  inherited by every branch descending from the shift, as in the models
  of Eastman et al. (2011), traitMedusa in motmot (Thomas & Freckleton
  2012), BAMM (Rabosky 2014) or bifrost; with `shift_branch = TRUE` it
  lasts only to the end of the branch. By default `shifts` shifts happen
  at the start of branches drawn with probabilities proportional to
  their length, or of the branches in `shift_edges`. With `shift_rate`,
  shifts instead happen as a Poisson process along the branches, as
  under the prior of BAMM: their number is Poisson with mean
  `shift_rate` times the total length of the tree, and they fall at
  uniform positions over it, so also within branches. At each shift the
  rate is multiplied by `shift_size` (`relative = TRUE`), or becomes
  `sigma2` times `shift_size` whatever it was before
  (`relative = FALSE`), as each new regime of BAMM or bifrost has a rate
  of its own. `shift_size` can be a function of the number of shifts
  that draws them, such as `function(n) stats::rlnorm(n, 0, 1)`. The
  rate of a branch is its mean over the branch, weighted by time.

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

Rabosky, D. L. (2014). Automatic detection of key innovations, rate
shifts, and diversity-dependence on phylogenetic trees. *PLoS ONE*, 9,
e89543.
[doi:10.1371/journal.pone.0089543](https://doi.org/10.1371/journal.pone.0089543)

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

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md),
[`plot_rates()`](https://provetelab.org/comdynamite/reference/plot_rates.md),
[`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md),
[`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md),
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md),
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md),
[`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md),
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(20)
x <- trait_rates(tr, "shift", shifts = 1, shift_size = 8)
x$shifts
#>   edge node position multiplier rate
#> 1   38   20        0          8    8
head(x$tips)
#>         t1         t2         t3         t4         t5         t6 
#>  0.3876716  0.3725845  0.1821453 -0.3773427 -0.5333873 -0.5354236 
# shifts as a Poisson process, each new regime with a rate of its own, as in BAMM
b <- trait_rates(tr, "shift", shift_rate = 0.1, relative = FALSE,
                 shift_size = function(n) stats::rexp(n))
b$shifts
#> [1] edge       node       position   multiplier rate      
#> <0 rows> (or 0-length row.names)
g <- trait_rates(tr, "gbm", rate_sd = 0.5)
summary(g$edge_rates)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.3990  0.7932  0.8618  0.8986  0.9967  1.7629 
# as the environmental optima of a metacommunity
traits <- species_traits(tr, optima = g$tips)
```
