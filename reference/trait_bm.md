# Trait evolution by Brownian motion with Pagel's lambda

Draws from a multivariate normal distribution with covariance \\\lambda
C + (1 - \lambda)\\ \mathrm{diag}(C)\\, where \\C\\ is the phylogenetic
covariance matrix of `tree`.

## Usage

``` r
trait_bm(tree, lambda = 1, nodes = FALSE)
```

## Arguments

- tree:

  A `phylo` tree with branch lengths.

- lambda:

  Pagel's lambda, in \[0, 1\].

- nodes:

  Also draw the values at the internal nodes? They are drawn after the
  tips, given them, so the tips are the same as with `nodes = FALSE` and
  the same seed. Under lambda, the nodes are those of the tree whose
  internal branches are multiplied by lambda and whose terminal branches
  are lengthened to keep the tips where they were.

## Value

Named numeric vector with the simulated trait value of each tip. With
`nodes = TRUE`, the attribute `nodes` holds the values at the internal
nodes, named by their numbers in `tree$edge`, the root first.

## Examples

``` r
set.seed(1)
trait_bm(treeshapesim::sim_tree_shape(5), lambda = 0.5)
#>           t1           t2           t3           t4           t5 
#> -0.005767173  2.275670676  0.763593461 -0.554890164 -1.062009857 
```
