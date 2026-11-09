# Phylogenetic beta diversity between communities

Mean phylogenetic distance between the species of each pair of
communities (COMDIST), the metric used by Peres-Neto & Kembel (2015).
Values in `comm` are the weights of the species, so incidence data give
the unweighted mean and abundances give the abundance-weighted mean.

## Usage

``` r
phylo_beta(comm, tree)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species, named as the
  tips of `tree`. Every site needs at least one individual. Or an
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  object, whose parts fill the arguments left out.

- tree:

  A `phylo` object with the same species as `comm`.

## Value

A `dist` object.

## Details

Two communities with the same species have a distance equal to their
mean pairwise distance, not zero, whereas the distance of a community to
itself is set to zero. Peres-Neto & Kembel (2015) describe the matrix as
fully Euclidean, but this does not always hold. If sites a and b both
have species x and y, site c has only x and site d has only y, then a
and b lie halfway between c and d and should coincide, yet their
distance is not zero. The pseudo-F of
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
does not require a Euclidean matrix (McArdle & Anderson 2001).

## References

McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to
community data: a comment on distance-based redundancy analysis.
*Ecology*, 82, 290–297.
[doi:10.1890/0012-9658(2001)082\[0290:FMMTCD\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282001%29082%5B0290%3AFMMTCD%5D2.0.CO%3B2)

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

## See also

Other the whole tree:
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
[`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md),
[`plot.mc_pcps()`](https://provetelab.org/comdynamite/reference/plot.mc_pcps.md),
[`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(5)
comm <- matrix(rbinom(20, 1, 0.7) + diag(1, 4, 5), 4,
               dimnames = list(NULL, tree$tip.label))
phylo_beta(comm > 0, tree)
#>          1        2        3
#> 2 1.168432                  
#> 3 1.097033 1.108714         
#> 4 1.251324 1.261214 1.251324
```
