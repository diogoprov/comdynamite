# Classic metrics of trait spacing and spread

Metrics of the distribution of trait values within a community used in
tests of community assembly: nearest-neighbour distances (Weiher et al.
1998; Stubbs & Wilson 2004), all pairwise distances, the minimum
spanning tree, and the range and kurtosis (Cornwell & Ackerly 2009).
Aiba et al. (2013) compared the statistical power of many of these
metrics. Each function returns one number, so that it can be given as
`index` to
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md);
with `multivariate = TRUE` there, the metrics of distances receive the
matrix of traits of each site.

## Usage

``` r
nn_distances(x, na.rm = FALSE)

mnnd(x, relative = FALSE, na.rm = FALSE)

sdnnd(x, relative = FALSE, na.rm = FALSE)

cvnnd(x, na.rm = FALSE)

min_nnd(x, relative = FALSE, na.rm = FALSE)

max_nnd(x, relative = FALSE, na.rm = FALSE)

mnd(x, relative = FALSE, na.rm = FALSE)

sdnd(x, relative = FALSE, na.rm = FALSE)

cvnd(x, na.rm = FALSE)

mst_length(x, relative = FALSE, na.rm = FALSE)

mst_ratio(x, na.rm = FALSE)

trait_range(x, na.rm = FALSE)

trait_kurtosis(x, na.rm = FALSE)
```

## Arguments

- x:

  Numeric vector (one trait), matrix or data frame (traits in columns,
  individuals or species in rows) or a `dist` object.

- na.rm:

  Drop missing values (rows with any missing value of a matrix) before
  the computation; with `FALSE` (the default) missing values give `NA`.

- relative:

  Divide by the range of the trait (one trait) or by the largest
  pairwise distance (several traits), as the `div_range` option of cati.

## Value

A single number, or `NA` when there are too few values: two for the
means, minima, maxima, the minimum spanning tree and the range, three
for the standard deviations and coefficients of variation, four for the
kurtosis. `nn_distances()` returns a vector with one distance per value
(all `NA` when a value is missing and `na.rm = FALSE`).

## Details

- `nn_distances()`: distance from each value to its nearest neighbour (a
  vector).

- `mnnd()`, `sdnnd()`, `cvnnd()`, `min_nnd()`, `max_nnd()`: mean,
  standard deviation, coefficient of variation (sd / mean), minimum and
  maximum of the nearest-neighbour distances.

- `mnd()`, `sdnd()`, `cvnd()`: mean, standard deviation and coefficient
  of variation of all pairwise distances.

- `mst_length()`: total length of the minimum spanning tree;
  `mst_ratio()`: shortest over longest edge of the tree, as in
  `MinMaxMST()` of the package cati.

- `trait_range()`: range of a single trait; `trait_kurtosis()`: excess
  kurtosis of a single trait (fourth standardised moment, with
  population moments, minus 3).

Distances are Euclidean among the rows of a matrix, or the absolute
differences of a vector; a `dist` object (for instance Gower distances
from [`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html))
is used as it is. Identical values are at distance zero and count as
such: an individual whose nearest neighbour has the same value has a
nearest-neighbour distance of zero (cati discards zero distances). For a
single trait the nearest neighbours and the minimum spanning tree are
found from the sorted values, without the matrix of distances, so that
large communities are fast.

Traits on different scales should be standardised over all communities
before the metrics are computed (for instance with
[`scale()`](https://rdrr.io/r/base/scale.html) on the whole data).
Standardising within each community, as cati does by default, changes
the units from one community to the next and between observed and null
communities.

## References

Aiba, M., Katabuchi, M., Takafumi, H., Matsuzaki, S. S., Sasaki, T., &
Hiura, T. (2013). Robustness of trait distribution metrics for community
assembly studies under the uncertainties of assembly processes.
*Ecology*, 94, 2873–2885.
[doi:10.1890/13-0269.1](https://doi.org/10.1890/13-0269.1)

Cornwell, W. K., & Ackerly, D. D. (2009). Community assembly and shifts
in plant trait distributions across an environmental gradient in coastal
California. *Ecological Monographs*, 79, 109–126.
[doi:10.1890/07-1134.1](https://doi.org/10.1890/07-1134.1)

Stubbs, W. J., & Wilson, J. B. (2004). Evidence for limiting similarity
in a sand dune community. *Journal of Ecology*, 92, 557–567.
[doi:10.1111/j.0022-0477.2004.00898.x](https://doi.org/10.1111/j.0022-0477.2004.00898.x)

Weiher, E., Clarke, G. D. P., & Keddy, P. A. (1998). Community assembly
rules, morphological dispersion, and the coexistence of plant species.
*Oikos*, 81, 309–322.
[doi:10.2307/3547051](https://doi.org/10.2307/3547051)

## See also

[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md).

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
x <- c(1, 2, 2, 5, 9)
nn_distances(x)
#> [1] 1 0 0 3 4
mnnd(x)
#> [1] 1.6
mst_length(x)
#> [1] 8
mnnd(cbind(x, rev(x)))
#> [1] 3.449242
# as an index of communities against a null model
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:5), each = 15),
                  species = sample(letters[1:4], 75, replace = TRUE), trait = rnorm(75))
community_index(ind, mnnd, null = "regional.ind", n_perm = 99)
#>   site trait index   observed        ses    p
#> 1   s1 trait  mnnd 0.21951480  1.5305257 0.12
#> 2   s2 trait  mnnd 0.11803015 -0.6224601 0.58
#> 3   s3 trait  mnnd 0.10781498 -0.9020233 0.42
#> 4   s4 trait  mnnd 0.06145019 -1.9225412 0.04
#> 5   s5 trait  mnnd 0.12445061 -0.4838972 0.70
```
