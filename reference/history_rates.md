# Rates of trait evolution from a history along the tree

Turns the history of a discrete character along a tree, such as the
habitat or region of each lineage from
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
or a stochastic map of phytools, into the rate of evolution of a
continuous trait on each branch, for `trait_rates(edge_rates = )`. The
rate is a function of the history at each moment, and the rate of a
branch is its mean over the branch, weighted by time. Three models link
the history to the rate:

## Usage

``` r
history_rates(
  history,
  model = c("state", "burst", "sympatry"),
  sigma2 = 1,
  state_rates = NULL,
  burst = 4,
  decay = 0.1,
  r = -0.1
)
```

## Arguments

- history:

  A tree with the history of a discrete character, of class `simmap`,
  from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  or phytools.

- model:

  `"state"`, `"burst"` or `"sympatry"`.

- sigma2:

  Rate before any change of state (`"burst"`) or scale of the rate
  (`"sympatry"`).

- state_rates:

  For `"state"`, the rate in each state, named by state.

- burst, decay:

  For `"burst"`, the relative increase of the rate at each change of
  state and the time scale of its decay.

- r:

  For `"sympatry"`, the change in the logarithm of the rate per lineage
  in the same state.

## Value

Numeric vector with the rate of each branch, in the order of the rows of
`history$edge`.

## Details

- `"state"`: each state has its own rate, `state_rates`, as in Brownian
  motion with several regimes (O'Meara et al. 2006).

- `"burst"`: each change of state, such as the colonisation of a new
  habitat or region, multiplies the rate by `1 + burst`, and the excess
  decays exponentially with time scale `decay`, also along the
  descending branches, until the next change. The rate at time \\u\\
  after the last change is \\\mathrm{sigma2}\\ (1 + \mathrm{burst}\\
  e^{-u / \mathrm{decay}})\\; before any change it is `sigma2`. This
  follows the expectation of McPeek (2017) that selection, and so
  evolution, is strongest when a lineage first adapts to a new
  community.

- `"sympatry"`: the rate depends on the number of lineages that share
  the state at each moment, counting the lineage itself, as
  \\\mathrm{sigma2}\\ e^{r n}\\, with \\n\\ that number; a negative
  \\r\\ slows evolution down as lineages accumulate in the same region,
  as in the diversity-dependent model with biogeography of Drury et al.
  (2016).

## References

Drury, J., Clavel, J., Manceau, M., & Morlon, H. (2016). Estimating the
effect of competition on trait evolution using maximum likelihood
inference. *Systematic Biology*, 65, 700–710.
[doi:10.1093/sysbio/syw020](https://doi.org/10.1093/sysbio/syw020)

McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton
University Press.
[doi:10.23943/princeton/9780691088778.001.0001](https://doi.org/10.23943/princeton/9780691088778.001.0001)

O'Meara, B. C., Ané, C., Sanderson, M. J., & Wainwright, P. C. (2006).
Testing for different rates of continuous trait evolution using
likelihood. *Evolution*, 60, 922–933.
[doi:10.1111/j.0014-3820.2006.tb01171.x](https://doi.org/10.1111/j.0014-3820.2006.tb01171.x)

## See also

[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md)

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(30)
habitat <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
rates <- history_rates(habitat, "state", state_rates = c(forest = 1, open = 5))
x <- trait_rates(habitat, edge_rates = rates)
# a burst of evolution after each change of habitat
b <- history_rates(habitat, "burst", burst = 9, decay = 0.1)
summary(b)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.000   1.000   1.000   1.393   1.000   5.078 
# slower evolution as lineages accumulate in the same habitat
s <- history_rates(habitat, "sympatry", r = -0.1)
```
