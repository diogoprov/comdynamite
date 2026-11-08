# Known biogeographic histories for Herodotools

Herodotools (Nakamura et al. 2024) estimates how assemblages were built
over evolutionary time: when lineages arrived in each region, from
where, and how traits changed along the tree. Its inputs are a tree, the
occurrence of species in sites, the region of each site and the area of
each ancestor. A simulation can give all of these with the history
known, so the estimates can be compared with the truth. This article
shows the steps and summarises `data-raw/validation_herodotools.R`.
Herodotools is not on CRAN, so its calls are shown but not run here;
install it from GitHub (`GabrielNakamura/Herodotools`) to run them.

## A history of areas and a metacommunity that follows it

The areas evolve along the tree as a discrete trait, with
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
whose states are the regions. The result keeps every change, in the
`simmap` format of phytools.

``` r

library(comdynamite)
set.seed(4)
tree <- treeshapesim::sim_tree_shape(30)
areas <- trait_mk(tree, q = 1, states = c("A", "B"))
table(areas$states)
#> 
#>  A  B 
#> 13 17
```

Each species can then occur only in the sites of its region, through
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md)
and the `pool` argument of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
Here 20 sites are split between the two regions.

``` r

regions <- rep(c("A", "B"), each = 10)
sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "phylo"), n_steps = 300,
                              pool = region_pool(areas, regions))
input <- as_herodotools(sim, regions, history = areas)
str(input, max.level = 1)
#> List of 4
#>  $ W             : num [1:20, 1:30] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ tree          :List of 4
#>   ..- attr(*, "class")= chr "phylo"
#>  $ biogeo        :'data.frame':  20 obs. of  1 variable:
#>  $ ancestral.area:'data.frame':  29 obs. of  1 variable:
```

[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md)
names its parts as the arguments of Herodotools, and `ancestral.area`
holds the true area of each node.

``` r

library(Herodotools)
arrival <- do.call(calc_age_arrival, input)
dispersal <- do.call(calc_dispersal_from, input)
insitu <- do.call(calc_insitu_diversification, input)
```

The true ages at which each lineage arrived in the region of each site
come from the history:

``` r

true_age <- arrival_ages(areas, sim, regions)
round(true_age[c(1, 11), 1:8], 2)
#>          t1 t2   t3  t4   t5   t6 t7 t8
#> site1    NA NA   NA 0.2 0.18 0.18 NA NA
#> site11 0.85 NA 0.85  NA   NA   NA NA NA
```

Herodotools dates an arrival by the age of the oldest ancestor node in
the region, so it can only see arrivals at nodes.

## Tip-based metrics of trait evolution

`calc_tip_based_trait_evo()` computes the metrics of Luza et al. (2021)
from stochastic character maps estimated from the tips.
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md)
computes them from the true history, and adds the number of changes
along the path to each tip, which includes changes that the nodes do not
show.

``` r

truth <- tip_history_metrics(areas)
head(truth)
#>    transitions total_nodes prop_transitions last_transition_time stasis_time
#> t1           1           4        0.2500000            1.0000000   0.3810036
#> t2           1           4        0.2500000            1.0000000   0.3810036
#> t3           1           3        0.3333333            1.0000000   0.3976755
#> t4           2           3        0.6666667            0.2484184   0.4023181
#> t5           2           5        0.4000000            0.2484184   0.4023181
#> t6           2           5        0.4000000            0.2484184   0.4023181
#>    changes time_since_change
#> t1       1         0.8512062
#> t2       1         0.8512062
#> t3       1         0.8512062
#> t4       2         0.1991536
#> t5       2         0.1828201
#> t6       2         0.1828201
```

``` r

estimate <- calc_tip_based_trait_evo(tree, areas$states, nsim = 20)
```

## Validation

**Same definitions.** We ran `calc_tip_based_trait_evo()` with the true
history in place of its stochastic maps, on 30 histories (25 species,
three states). Transitions, number of nodes, last transition time and
stasis time were identical to
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md)
for every tip (Herodotools 2.0.1).

**Recovery from the tips.** Trees of 40 species with height 1, two
states, 20 histories per rate of change, 20 stochastic maps averaged per
tip (symmetric model, rates by MCMC, as in
`calc_tip_based_trait_evo()`). Spearman correlation between the estimate
and the truth, and mean bias:

| Changes per tree height | Transitions  | Last transition time | Stasis time  |
|-------------------------|--------------|----------------------|--------------|
| 0.3                     | 0.68 (+0.09) | 0.72 (-0.01)         | 0.97 (-0.00) |
| 1                       | 0.45 (+0.09) | 0.67 (-0.04)         | 0.87 (-0.00) |
| 3                       | 0.38 (-0.48) | 0.44 (+0.06)         | 0.69 (+0.02) |

Even with the true history, transitions counted at nodes miss changes
inside branches: per tip, 0.28 against 0.31 changes at the lowest rate,
0.77 against 0.91 at the middle one and 1.81 against 3.01 at the
highest.

**Age of arrival.** The same histories, with the metacommunity simulated
as above (20 sites, 10 per region). `calc_age_arrival()` was given the
true ancestral areas or areas reconstructed by maximum likelihood
([`ape::ace()`](https://rdrr.io/pkg/ape/man/ace.html), equal rates, most
likely state at each node). Spearman correlation with the true age of
arrival, over species and sites, and mean bias:

| Changes per tree height | Nodes reconstructed correctly | True areas | Reconstructed areas |
|----|----|----|----|
| 0.3 | 99% | 0.88 (-0.03) | 0.58 (-0.02) |
| 1 | 91% | 0.89 (-0.07) | 0.53 (+0.03) |
| 3 | 69% | 0.70 (-0.02) | 0.24 (+0.22) |

With true areas, the loss comes from dating arrivals at nodes. With
reconstructed areas, one wrong node deep in the tree changes the age of
arrival of a whole clade in every site: in one history with 97% of the
nodes right, 99% of the ages changed. Mean ages of arrival were 0.86,
0.64 and 0.31.

These numbers depend on the size of the tree, the number of regions and
the method of reconstruction; the script can be changed to test others,
such as BioGeoBEARS.

## References

Luza, A. L., Maestri, R., Debastiani, V. J., Patterson, B. D., Hartz, S.
M., & Duarte, L. D. S. (2021). Is evolution faster at ecotones? A test
using rates and tempo of diet transitions in Neotropical Sigmodontinae
(Rodentia, Cricetidae). *Ecology and Evolution*, 11, 18676–18690.
<https://doi.org/10.1002/ece3.8476>

Nakamura, G., Rodrigues, A. V., Luza, A. L., Maestri, R., Debastiani, V.
J., & Duarte, L. D. S. (2024). Herodotools: An R package to integrate
macroevolution, biogeography and community ecology. *Journal of
Biogeography*, 51, 533–543. <https://doi.org/10.1111/jbi.14610>
