# Moran's eigenvector maps

Spatial predictors of Dray et al. (2006). Distances above the largest
edge of the minimum spanning tree are truncated, and the eigenvectors of
the centred similarity matrix with positive Moran's I are kept, in order
of decreasing spatial scale.

## Usage

``` r
spatial_eigenvectors(coords)
```

## Arguments

- coords:

  Matrix of coordinates with one row per site, or a `dist` object.

## Value

Matrix with one row per site and one column per map, with their Moran's
I in the attribute `moran`.

## Article and original code

The function `MEM` of the supplement differs from Appendix A of
Peres-Neto et al. (2012) in two points, and here we follow the code:

- Appendix A recommends centring with site weights so that the maps stay
  orthogonal in a weighted regression, and the code does not.
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  standardises the maps with site weights, and orthogonality is not
  needed to partition variation.

- Appendix A sets the distances between sites that are not connected to
  zero, and the code to four times the threshold.

## References

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

Dray, S., Legendre, P., & Peres-Neto, P. R. (2006). Spatial modelling: a
comprehensive framework for principal coordinate analysis of neighbour
matrices (PCNM). *Ecological Modelling*, 196, 483–493.
[doi:10.1016/j.ecolmodel.2006.02.015](https://doi.org/10.1016/j.ecolmodel.2006.02.015)

## Examples

``` r
dim(spatial_eigenvectors(cbind(1:20)))
#> [1] 20  9
```
