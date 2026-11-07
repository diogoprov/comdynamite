# Summary statistics of a metacommunity

The same statistics for simulated and observed data, so that they can be
compared or used for inference. Each statistic that lacks its inputs is
`NA`, so the result always has the same names and results can be bound
by rows.

## Usage

``` r
summary_stats(
  x,
  tree = NULL,
  trait = NULL,
  env = NULL,
  coords = NULL,
  n_null = 99
)
```

## Arguments

- x:

  An
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  object, an `mc_sim` object, a matrix with one row per site and one
  column per species, or an array of sites x species x steps. In an
  array, the last step gives the statistics that are not temporal. From
  an `mc_sim`, these come from `abundance` and the temporal ones from
  `history`, and `env`, `coords` and `trait` default to those of the
  simulation, the trait being the environmental optimum. From an
  `mc_data` object, they default to its parts: the trait is the column
  `optimum`, or the only numeric column of `traits`, and `env` the only
  numeric column of `env`.

- tree:

  Optional `phylo` object, matched to species by name.

- trait:

  Optional numeric vector with one value per species, named by species.

- env:

  Optional numeric vector with one value per site.

- coords:

  Optional matrix or data frame with one row per site and the
  coordinates in columns.

- n_null:

  Number of shuffles of species for `ses_mpd` and `ses_mntd`.

## Value

Named numeric vector.

## Details

`alpha`, `gamma` and `beta` are the mean number of species per site, the
number of species in all sites and `gamma / alpha`. `occupancy` and
`occupancy_sd` are the mean and the standard deviation of the share of
sites each species occupies, over the species present. `abundance` is
the mean number of individuals per site, `evenness` the mean Pielou
evenness of sites with two or more species, `evenness_regional` the same
for all sites together and `sad_sd` the standard deviation of the
logarithm of regional abundances, which describes the species abundance
distribution. With 0/1 data, these last four describe frequencies of
occurrence.

`ses_mpd` and `ses_mntd` need `tree`. They are the standardised effect
sizes of the mean pairwise distance and of the mean nearest taxon
distance, from presences, averaged over sites with two or more species.
The null model shuffles species across the tips, as
`null.model = "taxa.labels"` of
[`picante::ses.mpd()`](https://rdrr.io/pkg/picante/man/ses.mpd.html), so
the pool is all the species shared by `x` and `tree`, including those
without individuals. Compare results only between data with the same
pool. This is the only part that uses random numbers.

`decay` needs `coords`. It is the correlation between Bray-Curtis
dissimilarity and distance among occupied sites. `trait_env` needs
`trait` and `env`: the correlation between the mean trait of each site,
weighted by abundance, and the environment. `trait_sd` is the standard
deviation of the trait among all individuals. `trait_k` needs `trait`
and `tree`: Blomberg's K (Blomberg et al. 2003) of the trait over all
the species shared by `trait` and `tree`, including those without
individuals, as
[`picante::Kcalc()`](https://rdrr.io/pkg/picante/man/Kcalc.html) gives
it. For a simulation the trait is the optimum, which is not observed in
the field; the K of a measured trait carries less of the signal of the
optima.

The temporal statistics need at least two steps, from `history` or an
array. `cv_population`, `cv_site` and `cv_regional` are mean
coefficients of variation over time: of each species in each site, of
the individuals of each site, and of each species in all sites, leaving
out those that are always absent. `turnover` is the mean Bray-Curtis
dissimilarity of a site between consecutive steps. Guzman et al. (2022)
found that temporal statistics added explanatory power, and that
dispersal depended on them.

## References

Blomberg, S. P., Garland, T., & Ives, A. R. (2003). Testing for
phylogenetic signal in comparative data: behavioral traits are more
labile. *Evolution*, 57, 717–745.
[doi:10.1111/j.0014-3820.2003.tb00285.x](https://doi.org/10.1111/j.0014-3820.2003.tb00285.x)

Guzman, L. M., Thompson, P. L., Viana, D. S., Vanschoenwinkel, B.,
Horváth, Z., Ptacnik, R., Jeliazkov, A., Gascón, S., Lemmens, P.,
Anton-Pardo, M., Langenheder, S., De Meester, L., & Chase, J. M. (2022).
Accounting for temporal change in multiple biodiversity patterns
improves the inference of metacommunity processes. *Ecology*, 103,
e3683. [doi:10.1002/ecy.3683](https://doi.org/10.1002/ecy.3683)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(15)
sim <- simulate_metacommunity(landscape_grid(side = 5), species_traits(tree),
                              competition_matrix(tree, "phylo"), n_steps = 60,
                              seed_steps = 1:5, record_steps = seq(40, 60, by = 5))
round(summary_stats(sim, tree, n_null = 19), 3)
#>             alpha             gamma              beta         occupancy 
#>             5.720            13.000             2.273             0.440 
#>      occupancy_sd         abundance          evenness evenness_regional 
#>             0.334           239.720             0.766             0.841 
#>            sad_sd           ses_mpd          ses_mntd             decay 
#>             1.185             0.185            -0.251             0.163 
#>         trait_env          trait_sd           trait_k     cv_population 
#>             0.951            19.708             0.720             0.771 
#>           cv_site       cv_regional          turnover 
#>             0.082             0.275             0.107 
```
