# Traits evolving by an Ornstein-Uhlenbeck process with one or several optima

Simulates traits along a tree under the Ornstein-Uhlenbeck (OU) model of
Hansen (1997), with the multi-optimum ("multi-peak") extension of Butler
& King (2004). Along each branch the trait is pulled towards the optimum
`theta` of the selective regime of that branch with strength `alpha`,
while drifting with rate `sigma2`. Values are drawn exactly from the
transition distribution of the process: over a branch of length t, a
trait with value x at the start has, at the end, a normal distribution
with mean `theta + (x - theta) * exp(-alpha * t)` and variance
`sigma2 / (2 * alpha) * (1 - exp(-2 * alpha * t))`.

## Usage

``` r
traits_ou(
  tree,
  n_traits = 1,
  alpha = 0,
  sigma2 = 1,
  theta = 0,
  regimes = NULL,
  root = 0
)
```

## Arguments

- tree:

  A rooted `phylo` tree with branch lengths.

- n_traits:

  Number of independent traits.

- alpha:

  Strength of selection towards the optimum (0 = Brownian motion).

- sigma2:

  Rate of the stochastic (drift) component.

- theta:

  Optimum of each regime: one value for a single optimum, or a named
  vector with one value per regime (names as in `regimes`).

- regimes:

  Selective regime of each branch, in the order of the edges of `tree`
  reordered "cladewise" (`ape::reorder.phylo(tree, "cladewise")`), as
  returned by
  [`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md).
  `NULL` puts every branch in one regime.

- root:

  Trait value at the root.

## Value

Matrix with one row per tip (tip labels as row names) and one column per
trait.

## Details

With `alpha = 0` the process is Brownian motion with rate `sigma2`, and
the optima play no role. Each branch belongs to a single regime, so a
shift in optimum happens at the start of the branch where the regime
changes.

## References

Butler, M. A., & King, A. A. (2004). Phylogenetic comparative analysis:
a modeling approach for adaptive evolution. *The American Naturalist*,
164, 683–695. [doi:10.1086/426002](https://doi.org/10.1086/426002)

Hansen, T. F. (1997). Stabilizing selection and the comparative analysis
of adaptation. *Evolution*, 51, 1341–1351.
[doi:10.1111/j.1558-5646.1997.tb01457.x](https://doi.org/10.1111/j.1558-5646.1997.tb01457.x)

## See also

[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md)
to assign regimes to clades.

## Examples

``` r
set.seed(1)
tr <- ape::rcoal(20)
# Brownian motion
traits_ou(tr, alpha = 0)[1:3, ]
#>        t8        t7       t11 
#> 0.7324355 1.0476049 0.6082487 
# two optima: the clade of node 25 evolves towards 10, the rest towards 0
reg <- paint_regimes(tr, nodes = 25)
x <- traits_ou(tr, alpha = 2, theta = c(`0` = 0, `1` = 10), regimes = reg)
round(x[, 1], 1)
#>   t8   t7  t11  t17   t4  t13  t20   t9  t15  t19  t16  t18   t5  t14  t10  t12 
#>  0.5  0.5  0.6  0.7  0.6  0.4  6.8  6.9  0.8  0.1  0.5 -0.3 -0.3 -0.3  0.0  0.2 
#>   t6   t1   t3   t2 
#>  0.2  0.3 -0.1 -0.2 
```
