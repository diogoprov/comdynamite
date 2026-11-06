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
trait_bm(ape::rcoal(5), lambda = 0.5)
#>         t2         t5         t4         t3         t1 
#>  1.0262565  0.6694432 -0.4217225 -1.5578819  0.4996051 
```
