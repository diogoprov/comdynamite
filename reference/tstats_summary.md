# Summary of the null-model tests of T-statistics

Summarises the result of
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
with null models: for each trait and statistic, the share of sites below
and above the null distribution, the mean standardised effect size and a
test of the mean over sites, and the sites that deviate.

## Usage

``` r
tstats_summary(x, level = 0.05)
```

## Arguments

- x:

  Data frame from
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  with `nperm > 0`.

- level:

  Size of the two-sided test of each site and of the mean. A site can
  deviate only if `1 / (1 + nperm) < level / 2`, so `nperm` must be at
  least `2 / level`; a warning is given otherwise.

## Value

A list with two data frames. `traits` has one row per trait and
statistic, with the columns `trait`, `statistic`, `n_sites` (testable
sites), `percent_low`, `percent_high`, `mean_ses`, `p_mean_low`,
`p_mean_high`, `sites_low` and `sites_high` (site names separated by
spaces). `sites` has one row per site, trait and statistic, with
`observed`, `ses`, `p_low`, `p_high` and `deviation` (`"low"`, `"high"`,
`"none"` or `NA`).

## Details

Each site is tested in both directions from the null values kept by
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md):
the p-value of a low value is `(1 + r) / (1 + n)`, where `r` is the
number of null values at or below the observed one and `n` the number of
null values, and that of a high value counts null values at or above it.
A site deviates when one of the two is below `level / 2`, so that the
two tests together have size `level`. A site cannot be tested when its
observed value is missing or when every null value equals the observed
one, as for `T_IP.IC` at a site with one species, where the statistic is
1 whatever the permutation; such sites get `NA` and are left out of the
shares, instead of a small p-value in both directions.

The test of the mean compares the mean of the observed statistic over
the testable sites with the mean over the same sites of each
randomisation, so that the null distribution of the mean keeps the
dependence among sites within a randomisation.

## See also

[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md).

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
res <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 99)
tstats_summary(res)$traits
#>   trait statistic n_sites percent_low percent_high    mean_ses p_mean_low
#> 1 trait   T_IP.IC      15   100.00000     0.000000 -14.3546962       0.01
#> 2 trait   T_IC.IR      15    26.66667     6.666667  -0.8778137       0.04
#> 3 trait   T_PC.PR      15     0.00000     0.000000   0.3734805       0.89
#>   p_mean_high
#> 1        1.00
#> 2        0.97
#> 3        0.12
#>                                                                                         sites_low
#> 1 site1 site2 site3 site4 site5 site6 site7 site8 site9 site10 site11 site12 site13 site14 site15
#> 2                                                                        site2 site4 site5 site12
#> 3                                                                                                
#>   sites_high
#> 1           
#> 2      site6
#> 3           
```
