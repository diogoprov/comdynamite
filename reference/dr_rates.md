# Tip diversification rates (DR) of Jetz et al. (2012)

`dr_rates()` computes the DR statistic of Jetz et al. (2012), a measure
of the recent rate of speciation of each tip from the tree alone: the
inverse of its equal-splits measure of evolutionary distinctiveness,
`DR_i = 1 / sum_j l_j / 2^(j - 1)`, where `l_j` is the length of the
`j`-th edge on the path from tip `i` to the root, the terminal edge
first. Tips in recently and rapidly diversifying clades have short,
heavily shared paths and a high DR.

## Usage

``` r
dr_rates(tree, per = c("tip", "edge"))
```

## Arguments

- tree:

  A rooted `phylo` tree with branch lengths.

- per:

  `"tip"` for one value per tip, named by the tips, or `"edge"` for one
  value per edge, in the order of `tree$edge`.

## Value

A numeric vector of DR values.

## Details

With `per = "edge"`, each edge gets the geometric mean of the DR of the
tips below it (the DR of its tip for a terminal edge), a value per edge
in the order of `tree$edge` that can be given as the `edge_covariate` of
[`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
or
[`damocles_samples()`](https://provetelab.org/comdynamite/reference/damocles_samples.md).
DR is defined for tips; the value of an internal edge summarises its
descendants and is not an estimate of the rate along that edge, unlike
the branch-specific rates of model-based methods such as ClaDS.
`data-raw/damocles_clads.R` checks the tip values against the
equal-splits measure of
[`picante::evol.distinct()`](https://rdrr.io/pkg/picante/man/evol.distinct.html)
and compares DR with ClaDS as covariates.

## References

Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O.
(2012). The global diversity of birds in space and time. *Nature*, 491,
444–448. [doi:10.1038/nature11631](https://doi.org/10.1038/nature11631)

## See also

[`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
[`damocles_samples()`](https://provetelab.org/comdynamite/reference/damocles_samples.md).

## Examples

``` r
set.seed(10)
tr <- treeshapesim::sim_tree_shape(20)
head(dr_rates(tr))
#>       t1       t2       t3       t4       t5       t6 
#> 2.227964 2.227964 1.436217 2.056434 2.056434 1.163367 
head(dr_rates(tr, per = "edge"))
#> [1] 1.924622 2.227964 2.227964 2.227964 1.436217 2.768076
```
