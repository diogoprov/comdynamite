# Phylogenetic eigenvectors

Principal coordinates of the patristic distances among species, the
first step of the phylogenetic community composition of Peres-Neto &
Kembel (2015).

## Usage

``` r
phylo_eigenvectors(tree, squared = FALSE)
```

## Arguments

- tree:

  A `phylo` object.

- squared:

  Use \\-\mathrm{PD}^2 / 2\\, as in a classic principal coordinate
  analysis, instead of \\-\mathrm{PD} / 2\\, the form printed by
  Peres-Neto & Kembel (2015)?

## Value

Matrix with one row per tip and one column per positive eigenvalue.

## References

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

## Examples

``` r
set.seed(1)
dim(phylo_eigenvectors(treeshapesim::sim_tree_shape(10)))
#> [1] 10  9
```
