# Trait evolution by Brownian motion with Pagel's lambda

Draws from a multivariate normal distribution with covariance
`lambda * C + (1 - lambda) * diag(C)`, where `C` is the phylogenetic
covariance matrix of `tree`.

## Usage

``` r
trait_bm(tree, lambda = 1)
```

## Arguments

- tree:

  A `phylo` tree with branch lengths.

- lambda:

  Pagel's lambda, in \[0, 1\].

## Value

Named numeric vector with the simulated trait value of each tip.

## Examples

``` r
set.seed(1)
trait_bm(treeshapesim::sim_tree_shape(5), lambda = 0.5)
#>           t1           t2           t3           t4           t5 
#> -0.005767173  2.275670676  0.763593461 -0.554890164 -1.062009857 
```
