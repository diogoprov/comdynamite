# Validating the methods of Leibold et al. (2010) and Peres-Neto et al. (2012)

The functions
[`fourth_corner_test()`](https://provetelab.org/comdynamite/reference/fourth_corner_test.md),
[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md),
[`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md),
[`spatial_eigenvectors()`](https://provetelab.org/comdynamite/reference/spatial_eigenvectors.md)
and
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
are translations of MATLAB code from the supplements of these articles.
We checked them in three ways. Each help page lists the differences
between the text of the article and the code.

## Comparison with the original code

We ran the MATLAB code in Octave on the numerical example of Appendix A
of Peres-Neto et al. (2012) and on simulated data. The supplement lacks
the functions `eigen` and `WeightedMeanVar`, which we wrote. Weighted
coordinates, correlations, Moran’s eigenvector maps and the raw and
adjusted fractions agreed to within 1e-13. Correlations also agreed with
`ade4::fourthcorner()`, and maps with `adespatial::dbmem()`.

## Type I error and power

The script `data-raw/validation_gradient.R` follows the simulations of
the article, with 1000 data sets and 999 permutations (100 sites and 50
species, one environmental variable). The article stacks the two null
distributions and uses alpha / 2, shown here as “stacked”. Rates near
0.05 are expected for a true null hypothesis, with a standard error near
0.007.

Trait and environment (Table 4 of the article):

| Case | Tolerance | Row, article | Row, here | Column, article | Column, here | Stacked, article | Stacked, here | Larger p, here |
|----|----|----|----|----|----|----|----|----|
| Random trait | 1 | 0.836 | 0.812 | 0.056 | 0.050 | 0.030 | 0.048 | 0.050 |
| Random trait | 20 | 0.916 | 0.923 | 0.052 | 0.042 | 0.046 | 0.042 | 0.042 |
| Random trait | 60 | 0.984 | 0.985 | 0.058 | 0.048 | 0.048 | 0.047 | 0.048 |
| Random environment | 1 | 0.052 | 0.053 | 0.678 | 0.695 | 0.052 | 0.052 | 0.053 |
| Random environment | 20 | 0.056 | 0.052 | 0.832 | 0.866 | 0.038 | 0.051 | 0.052 |
| Random environment | 60 | 0.058 | 0.046 | 0.954 | 0.962 | 0.036 | 0.046 | 0.046 |

Global test of
[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md),
with one Brownian trait:

| Case | Tolerance | Stacked, article | Stacked, here | Larger p, here |
|----|----|----|----|----|
| Random environment (Table 1) | 1 | 0.028 | 0.039 | 0.051 |
| Random phylogeny (Table 2) | 1 | 0.042 | 0.040 | 0.040 |
| Power (Table 3) | 1 | 0.880 | 0.875 | 0.877 |
| Power (Table 3) | 20 | 0.948 | 0.934 | 0.936 |
| Power (Table 3) | 60 | 0.948 | 0.940 | 0.941 |

With tolerances of 20 and 60, which the article does not report for
Tables 1 and 2, the larger p-value gave rates from 0.037 to 0.055. The
test of the final vector after a significant global test kept these
rates, because the global test comes first. A test of that vector alone
is not valid, as it was chosen to fit the predictor.

### Weights of the eigenvectors

The correlation between the vector and the trait that generated the data
depends on how the eigenvectors are weighted. The text uses correlations
and the code their absolute values.

| Tolerance | PEV with absolute values, here | PEV in Table 3 | PEV with signs, here | PEV-B, here | PEV-B in Table 3 |
|----|----|----|----|----|----|
| 1 | 0.480 | 0.514 | 0.685 | 0.698 | 0.676 |
| 20 | 0.455 | 0.526 | 0.722 | 0.707 | 0.703 |
| 60 | 0.473 | 0.508 | 0.740 | 0.715 | 0.688 |

Absolute values give about what the article reports for PEV, which
suggests that its Table 3 used them (this is our interpretation). With
signs, PEV matches PEV-B, so the gain of the Brownian vectors reported
in the article may come from this difference.

### Partition

We did not reproduce Table 5, which depends on details of Appendix C.
Following the demonstration script of the authors instead, with a vector
fixed in advance and predictors without relation to the data, rejection
rates of the larger p-value were 0.023 to 0.028, because both null
hypotheses are true. With environment and space acting, the global test
for environment passed in 37% of the data sets. In those, the rejection
rate was 0.99 for fraction `a` and 0.03 for fraction `c`, as expected
for a vector estimated for environment alone.

## The daphniids of Leibold et al. (2010)

The data set `daphniids` rebuilds the data of the article from the
public survey. With the larger p-value, no predictor is significant at
any node, as species permutations have little power with 13 species. The
article permutes rows and columns together. With sites alone, at a level
of 0.05 without correction and 999 permutations, environmental results
compare with Figure 3a as follows:

| Node | Figure 3a | Here |
|----|----|----|
| *retrocurva*, *parvula*, *catawba* | COND, NTL | COND, NTL |
| Same plus *pulicaria*, *pulex* | AV_DEP, PHEQ | AV_DEP, PHEQ |
| *longispina*, *galeata*, *rosea* | AV_DEP, ELEV | AV_DEP, ELEV |
| All species | TN:TP, AV_DEP, ELEV, DOC, COLOR | COLOR, DOC, AV_DEP |
| All *Daphnia* except *laevis*, *dubia* | AV_DEP, PHEQ, ELEV, DOC | AV_DEP, PHEQ, ELEV, DOC, NTL, TN_TP |
| Clade with *ambigua* | ELEV | ELEV, TN_TP |
| *retrocurva*, *parvula* | PHEQ, AV_DEP | AV_DEP |
| *pulicaria*, *pulex* | ELEV, NTL | none |
| *galeata*, *rosea* | COLOR, DOC | none |
| *laevis*, *dubia* | COND | NTL, PHEQ, PTL, TN_TP |

Nodes with rare species differ most: *D. laevis* and *D. rosea* occur in
one and two lakes of the rebuilt set. Biogeographic results differ more
(no match at most nodes), which we expected, because the article does
not give the area of each lake and we approximated rivers by longitude.

## References

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
<https://doi.org/10.1111/j.1461-0248.2010.01523.x>

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
<https://doi.org/10.1890/11-0494.1>
