# Inference by ABC with random forests

The simulations of comdynamite can be turned around: instead of asking
what patterns a set of processes produces, we can ask which values of
the parameters of a process make simulated metacommunities look like an
observed one. Approximate Bayesian computation (ABC) does this by
simulating many data sets with parameters drawn from a prior, describing
each with summary statistics, and keeping the parameters whose
statistics resemble the observed ones. With random forests (Pudlo et
al. 2016; Raynal et al. 2019), there is no tolerance or distance to
choose: a forest learns how the parameters relate to the statistics in
the reference table and weighs its rows by their similarity to the
observed data.

This vignette goes through the whole workflow with the anurans of 40
urban ponds in two Brazilian cities (`anurans`):

1.  build a landscape from the coordinates of the ponds and write the
    simulator;
2.  simulate a reference table with
    [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md);
3.  check that the parameters can be recovered at all, with
    [`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md);
4.  choose between models with
    [`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md);
5.  estimate the parameters for the observed data with
    [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
    and check the fit with
    [`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
    [`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md)
    and cross-validation over spatial blocks with
    [`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md).

The reference tables (about 20,000 simulations) take about a quarter of
an hour on six cores, and they and the forests were computed in advance
with
[`data-raw/abc_anurans.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/abc_anurans.R);
the rest runs here.

``` r

library(comdynamite)
library(ggplot2)
res <- readRDS("abc.rds")
prob <- function(p) if (round(p, 2) >= 1) "> 0.99" else round(p, 2)
```

## The data and the simulator

The ponds of each city are a few kilometres apart, and the two cities
are about 900 km from each other, so we treat them as two
metacommunities that share the values of the parameters. Each pond has a
hydrology score, from small, shallow and warm ponds (0) to large and
deep ones (100), and the urban cover within 500 m. Each city has its own
regional pool: the anurans that breed in lentic water and were recorded
by GBIF in its ecoregion and within 200 km of its ponds (see
[`?anurans`](https://provetelab.org/comdynamite/reference/anurans.md)).
The pools are much larger than the sets of species found in the ponds,
so the simulations must explain which species of the region reach and
persist in the ponds.

``` r

table(anurans$city)
#> 
#>        Campo Grande São José dos Campos 
#>                  20                  20
c(recorded = ncol(anurans$comm), pool = ncol(anurans$pool))
#> recorded     pool 
#>       44      155
tapply(rowSums(anurans$pool), anurans$city, max)
#>        Campo Grande São José dos Campos 
#>                  53                 146
city <- ifelse(anurans$city == "Campo Grande", "cg", "sjc")
pool <- anurans$pool
rownames(pool) <- paste0("site", seq_along(city))
```

[`landscape_points()`](https://provetelab.org/comdynamite/reference/landscape_points.md)
keeps the real layout of the ponds. With `lonlat = TRUE`, distances are
great-circle distances in kilometres, and `disp_scale` is the scale of
the exponential dispersal kernel in kilometres.

``` r

ponds <- landscape_points(anurans$coords, anurans$env[, "hydrology"], disp_scale = 2,
                          lonlat = TRUE)
ponds
#> <mc_landscape> 40 sites; points, exponential dispersal kernel 
#> environment:    0 to 100
```

The model has six parameters: the fraction of individuals that emigrate
at each step (`log_dispersal`, on a log10 scale), the scale of the
dispersal kernel (`log_kernel`, log10 km), the strength of competition
between the most closely related species (`inter_max`, with competition
decreasing with phylogenetic distance), the mean niche breadth along the
hydrology axis (`sigma_mean`), the phylogenetic signal of the niche
optima (`lambda`, Pagel’s lambda), and the share of species that
tolerate any urban cover (`q_tol`). Following the chapter from which the
data come, each of the other species gets a threshold drawn from U(0,
100)% and is kept out of ponds whose urban cover exceeds it; the filter
acts on the pool of each pond, before the dynamics. The simulator takes
one draw of the parameters (`scenario`) and the objects shared by all
simulations (`shared`):

``` r

shared <- list(tree = anurans$tree, coords = anurans$coords, env = anurans$env[, "hydrology"],
               urban = anurans$env[, "urban_cover"], pool = pool, city = city)
prior <- list(log_dispersal = c(-4, log10(0.3)), log_kernel = log10(c(0.5, 20)),
              inter_max = c(0, 1), sigma_mean = c(5, 100), lambda = c(0, 1), q_tol = c(0, 1))
urban_pool <- function(pool, urban, q_tol) {
  threshold <- ifelse(runif(ncol(pool)) < q_tol, Inf, runif(ncol(pool), 0, 100))
  filtered <- pool * outer(urban, threshold, "<=")
  lost <- colSums(filtered) == 0  # a species excluded everywhere keeps its pool
  filtered[, lost] <- pool[, lost]
  filtered
}
simulate <- function(scenario, shared) {
  landscape <- landscape_points(shared$coords, shared$env, 10^scenario$log_kernel, lonlat = TRUE)
  simulate_metacommunity(landscape,
                         species_traits(shared$tree, scenario$lambda,
                                        sigma_mean = scenario$sigma_mean,
                                        sigma_sd = scenario$sigma_mean / 2),
                         competition_matrix(shared$tree, "phylo",
                                            inter_max = scenario$inter_max),
                         dispersal = 10^scenario$log_dispersal, n_steps = 500,
                         pool = urban_pool(shared$pool, shared$urban, scenario$q_tol))
}
set.seed(1)
sim <- simulate(list(log_dispersal = -2, log_kernel = 0.3, inter_max = 0.5, sigma_mean = 20,
                     lambda = 0.5, q_tol = 0.3), shared)
#> simulate_metacommunity() implements the method of:
#>   Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik, R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A process-based metacommunity framework linking local and regional scale community ecology. Ecology Letters, 23, 1314-1329. https://doi.org/10.1111/ele.13568
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
rowSums(sim$abundance > 0)[1:10]
#>  site1  site2  site3  site4  site5  site6  site7  site8  site9 site10 
#>     13     11     11      7     19     15     11      6      5     14
```

## Summary statistics

The observed data are counts from surveys, whose totals depend on effort
and on how easily each species is detected, so we describe both the
observed and the simulated metacommunities by incidence alone: the
statistics of
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
that need no abundances (mean and total richness, beta diversity,
occupancy, the standardised effect sizes of MPD and MNTD, and distance
decay), plus the correlations of richness and of compositional
dissimilarity with hydrology. They are computed separately for each
city, with the prefixes `cg_` and `sjc_`, over the species of its
regional pool, which is also the pool of the null model of the
standardised effect sizes. See `data-raw/abc_anurans.R` for the
function; the observed statistics are

``` r

round(res$observed[grepl("^(cg|sjc)_", names(res$observed))], 3)
#>            cg_alpha            cg_gamma             cg_beta        cg_occupancy 
#>               9.450              27.000               2.857               0.350 
#>     cg_occupancy_sd          cg_ses_mpd         cg_ses_mntd            cg_decay 
#>               0.258              -0.758              -0.129               0.154 
#>     cg_richness_env  cg_composition_env           sjc_alpha           sjc_gamma 
#>               0.231              -0.052               8.000              26.000 
#>            sjc_beta       sjc_occupancy    sjc_occupancy_sd         sjc_ses_mpd 
#>               3.250               0.308               0.276              -0.235 
#>        sjc_ses_mntd           sjc_decay    sjc_richness_env sjc_composition_env 
#>               0.105               0.234               0.113              -0.030
```

With the simulator and the statistics in hand,
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
draws parameters uniformly from the prior and runs one simulation per
draw, in parallel with future:

``` r

future::plan(future::multicore, workers = 6)
table <- abc_table(prior, simulate, n = 10200, setup = function() shared,
                   summarise = summarise, seed = 2040)
```

## Can the parameters be recovered?

Before looking at the observed data, we ask whether the statistics carry
information about each parameter.
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
leaves 200 simulations out of the reference table, treats each as if it
were observed, and compares the estimates with the true values. A
parameter is worth estimating when the estimates follow the truth (R²
well above zero) and the 90% intervals contain the truth at least about
90% of the time. The intervals of random forests are approximate: with
many statistics they tend to be somewhat wider than the exact ones
(Raynal et al. 2019), so a coverage a little above 90% is expected,
while one well below it would say that the intervals are too narrow to
be trusted.

``` r

print(res$recovery, digits = 2)
#>       parameter   n      r2    bias coverage coverage_se width
#> 1 log_dispersal 200 0.73298  0.0157     0.95       0.015  0.49
#> 2    log_kernel 200 0.24312 -0.0043     0.91       0.021  0.80
#> 3     inter_max 200 0.94467  0.0029     0.94       0.016  0.23
#> 4    sigma_mean 200 0.58988  0.0030     0.94       0.016  0.61
#> 5        lambda 200 0.00012 -0.0224     0.91       0.020  0.89
#> 6         q_tol 200 0.68877  0.0080     0.91       0.020  0.59
```

``` r

ggplot(res$estimates, aes(truth, estimate)) +
  geom_abline(colour = "grey60") +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "True value", y = "Posterior median") +
  theme_minimal()
```

![Estimated against true values of the five parameters for the
pseudo-observations, with the line of
equality.](abc_files/figure-html/recovery-plot-1.png)

Competition (`inter_max`) and emigration are well recovered, and niche
breadth and tolerance to urban cover reasonably so, while the scale of
the kernel is only weakly identified and the phylogenetic signal of the
niche (`lambda`) not at all: with incidence in 40 ponds, its posterior
will be its prior. The coverage is at or a little above 90% for all six
(0.9 to 0.96, with a standard error of about 0.02), so the intervals are
not too narrow and the posterior of `lambda` is honest about this, just
uninformative. This recovery holds for the tree and the ponds of the
anurans only. In simulations on small grids
(`data-raw/abc_variable_design.R`), the recovery of emigration,
competition, niche breadth and detection changed little when every
simulation drew its own tree and landscape, but that of `lambda`,
already near the limit with one tree (R² 0.43 to 0.53), fell to 0.32.
What the forests learn about `lambda` therefore seems to depend on the
tree used, and a recovery study with the observed tree, as here, is the
one to trust.

## Which model?

The same statistics may come from metacommunities without competition.
We simulated a second table with `inter_max` fixed at zero (the other
five parameters from the same priors) and let
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md)
tell the two models apart:

``` r

res$choice$error
#> [1] 0.005462564
res$choice$confusion
#>              phylogenetic none
#> phylogenetic         9869   87
#> none                   22 9976
print(res$chosen, digits = 2)
#>   observation        model probability votes_phylogenetic votes_none
#> 1           1 phylogenetic           1                491          9
```

The out-of-bag error shows that the two models are easy to tell apart
when one of them is true. For the observed data, 491 of the 500 trees
vote for the model with phylogenetic competition, and the posterior
probability of that model is \> 0.99.

## The posterior

[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
grows one regression forest per parameter, and
[`predict()`](https://rdrr.io/r/stats/predict.html) gives the posterior
median and a 90% interval for the observed statistics:

``` r

print(res$posterior, digits = 2)
#>   observation     parameter estimate    low  high  mean
#> 1           1 log_dispersal   -1.049 -3.545 -0.62 -1.41
#> 2           1    log_kernel    0.017 -0.270  0.63  0.10
#> 3           1     inter_max    0.573  0.354  0.95  0.59
#> 4           1    sigma_mean   31.786  6.611 98.36 44.87
#> 5           1        lambda    0.724  0.059  0.97  0.62
#> 6           1         q_tol    0.310  0.035  0.81  0.35
```

The rows of the reference table that the forests weigh give an
approximate sample from the posterior, which we can set against the
uniform prior:

``` r

draws <- stack(res$draws[names(res$prior)])
limits <- data.frame(ind = names(res$prior), low = sapply(res$prior, `[`, 1),
                     high = sapply(res$prior, `[`, 2))
ggplot(draws, aes(values)) +
  geom_histogram(bins = 25, fill = "grey40") +
  geom_blank(aes(x = low), limits) + geom_blank(aes(x = high), limits) +
  facet_wrap(~ind, scales = "free") +
  labs(x = "Value (panel width = prior range)", y = "Draws") +
  theme_minimal()
```

![Histograms of posterior draws of the five parameters, with the prior
range shown by the horizontal extent of each
panel.](abc_files/figure-html/posterior-plot-1.png)

These draws use the mean of the weights of the six forests, so they are
broader than the intervals of each parameter.

## Checking the fit

A posterior is only as good as the model behind it.
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md)
asks whether the observed statistics lie within the range of the
reference table and how many rows carry the weight (effective sample
size, ESS):

``` r

print(res$diagnose, digits = 2)
#>   observation     parameter n_positive ess khat khat_threshold n_outside
#> 1           1 log_dispersal        319  81 0.49           0.60         1
#> 2           1    log_kernel        448 131 0.34           0.62         1
#> 3           1     inter_max        497 155 0.45           0.63         1
#> 4           1    sigma_mean        393  96 0.45           0.61         1
#> 5           1        lambda        377 146 0.17           0.61         1
#> 6           1         q_tol        362 101 0.46           0.61         1
#>      outside
#> 1 cg_ses_mpd
#> 2 cg_ses_mpd
#> 3 cg_ses_mpd
#> 4 cg_ses_mpd
#> 5 cg_ses_mpd
#> 6 cg_ses_mpd
```

The weight is spread over enough rows for every parameter, but one
observed statistic, the standardised MPD of Campo Grande, lies outside
the central 99% of the reference table: few simulations under the priors
produced ponds with species as closely related as the observed ones. The
posterior predictive check below shows the same misfit from another
angle.

A posterior predictive check compares the observed statistics with those
of the simulations drawn from the posterior. It is most telling for
statistics left out of the fit, so we fitted a second set of forests
without the phylogenetic statistics and checked those:

``` r

print(res$ppc_phylo, digits = 2)
#>      statistic observed median   low high     p   n
#> 1   cg_ses_mpd    -0.76   1.02  0.16  1.7 0.000 500
#> 2  cg_ses_mntd    -0.13   0.35 -0.74  1.4 0.400 500
#> 3  sjc_ses_mpd    -0.23   1.47  0.50  2.1 0.028 500
#> 4 sjc_ses_mntd     0.11   0.11 -1.05  1.4 0.996 500
```

The model chosen with a posterior probability of \> 0.99 fails here.
With competition among close relatives, the simulated ponds hold species
more distantly related than expected from the regional pool
(standardised MPD above zero), whereas in the observed ponds the species
are as related as, or more related than, random draws from the pool. The
observed MPD is in the lower tail in both cities (p \< 0.002 in Campo
Grande and p = 0.028 in São José dos Campos). With all the statistics in
the fit, the forests pull the posterior towards the observed MPD and the
misfit is milder:

``` r

print(res$ppc_all[c("statistic", "observed", "median", "low", "high", "p")], digits = 2)
#>              statistic observed median    low  high     p
#> 1             cg_alpha    9.450  7.500  4.200 14.87 0.408
#> 2             cg_gamma   27.000 19.000 11.000 34.00 0.264
#> 3              cg_beta    2.857  2.566  1.538  3.88 0.580
#> 4         cg_occupancy    0.350  0.390  0.258  0.65 0.580
#> 5      cg_occupancy_sd    0.258  0.252  0.166  0.31 0.872
#> 6           cg_ses_mpd   -0.758  0.757 -0.704  1.58 0.080
#> 7          cg_ses_mntd   -0.129  0.151 -0.723  1.09 0.648
#> 8             cg_decay    0.154  0.062 -0.102  0.23 0.304
#> 9      cg_richness_env    0.231  0.235 -0.056  0.52 0.976
#> 10  cg_composition_env   -0.052  0.053 -0.092  0.62 0.268
#> 11           sjc_alpha    8.000 10.550  6.195 27.40 0.440
#> 12           sjc_gamma   26.000 28.000 17.000 61.10 0.852
#> 13            sjc_beta    3.250  2.738  1.501  4.03 0.308
#> 14       sjc_occupancy    0.308  0.365  0.248  0.67 0.308
#> 15    sjc_occupancy_sd    0.276  0.264  0.170  0.33 0.776
#> 16         sjc_ses_mpd   -0.235  1.232 -0.085  1.97 0.068
#> 17        sjc_ses_mntd    0.105  0.081 -0.916  1.03 0.964
#> 18           sjc_decay    0.234  0.226  0.055  0.41 0.952
#> 19    sjc_richness_env    0.113 -0.211 -0.570  0.11 0.100
#> 20 sjc_composition_env   -0.030  0.209  0.030  0.66 0.020
```

Why, then, does model choice favour competition so strongly? The
statistics that weigh most in the forest of `inter_max` are the mean and
total richness of each city:

``` r

head(res$importance[res$importance$parameter == "inter_max", ], 4)
#>   parameter statistic importance
#> 1 inter_max sjc_gamma  0.3331401
#> 2 inter_max sjc_alpha  0.2184167
#> 3 inter_max  cg_gamma  0.1932001
#> 4 inter_max  cg_alpha  0.1246896
```

With pools of 50 to 150 species, the observed ponds hold few of them,
and among the processes in the simulator, competition is the one that
can keep richness that low. Any process that limits richness but is
missing from the model, such as colonisation from a region where most
species are rare, may therefore be read as competition. The phylogenetic
statistics, which measure who coexists rather than how many, seem to
point instead to a filter that favours some clades, for instance a
tolerance to urban cover shared by related species, whereas the filter
of the simulator acts at random with respect to phylogeny.

We tested the first of these hypotheses with simulations in
[`data-raw/abc_conserved_filter.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/abc_conserved_filter.R).
In a second version of the simulator, a trait drawn on the tree by
Brownian motion with Pagel’s lambda (`lambda_filter`, uniform between 0
and 1) ranks the species, so that related species have similar
thresholds of tolerance to urban cover; with `lambda_filter = 0` it
reduces to the random filter. With and without competition, this gives
four models, each with 10,000 simulations. The conserved filter does not
explain the misfit. For the observed ponds, the models with competition
received 497 of the 500 votes, 427 of them for competition with the
conserved filter (posterior probability 0.82), and under this model the
observed MPD stays in the lower tail (p = 0.016 in Campo Grande and p =
0.06 in São José dos Campos). The two filters are also hard to tell
apart: with competition, about 30% of the simulations of each filter
were assigned to the other, and `lambda_filter` was not recovered (R²
0.06), so the votes for the conserved filter say little about
conservatism itself.

The same simulations measured the power of the check. We fitted the
model of this vignette to 200 pseudo-observations from a conserved
filter without competition and checked their phylogenetic statistics.
The standardised MPD was flagged (p \< 0.05) in 39% of them in São José
dos Campos and in 23% in Campo Grande, against 4% when the model was
true; the MNTD was flagged no more often than under the true model.
Interestingly, the power hardly changed with `lambda_filter`, so the
check seems to detect the absence of competition, which the statistics
of richness in the fit do not reveal, more than the conservatism of the
filter. Model choice between the two models of this vignette assigned
199 of these 200 pseudo-observations to the model without competition,
so a conserved filter alone does not produce the pattern read here as
competition. The check says what the model misses, not why, and the
misfit of the MPD remains unexplained by a conserved filter.

The second hypothesis, limited colonisation from the region, was tested
in the same way in
[`data-raw/abc_colonisation.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/abc_colonisation.R).
In the simulator above, every species of the pool receives the same
number of seeds in every pond. In the alternative, the mean number of
seeds of each species (`seed_mean`, which
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
accepts per species) is `10^log_seed` times a regional abundance drawn
from a lognormal distribution with standard deviation `sigma_region` on
the log scale, so that a few species are common in the region and most
are rare. The log abundances were drawn on the tree with Pagel’s lambda
(`lambda_region`), so that related species may be similarly common. With
`sigma_region = 0` and `log_seed = log10(0.5)`, this is the model above.
Because the standardised MPD and MNTD of the observed ponds change with
the draws of the null model (standard deviation 0.10 to 0.19 among sets
of 19 draws), their observed values in this study were averaged over 100
sets of draws.

Unlike the conserved filter, limited colonisation is read as
competition. Among 199 pseudo-observations from limited colonisation
without competition, model choice between the two models of this
vignette assigned 103 to competition, and this share rose from 13% when
the regional abundances were nearly even (`sigma_region` below 1) to 89%
when they were most uneven (above 2). The posterior predictive check
flagged the standardised MPD of these pseudo-observations in 43% of them
in São José dos Campos and in 25% in Campo Grande, against 2% to 4% when
the model was true. With all four models, however, competition was
rarely chosen for these pseudo-observations (10 of 199), so the
alternative model, once included, removes most of the confusion. For the
observed ponds, the models with competition still received 454 of the
500 votes, and competition with limited colonisation was chosen with a
posterior probability of 0.60. Under this model, the estimate of
`inter_max` was lower than above (0.35, 90% interval 0.11 to 0.72), and
the observed MPD stayed in the lower tail, but less so (p = 0.04 in
Campo Grande and p = 0.08 in São José dos Campos). Limited colonisation
therefore seems to account for part of the low richness that the model
above attributed to competition, and for part of the misfit of the MPD,
but not for all of it. Its parameters were poorly recovered (R² 0.12 for
`log_seed`, 0.37 for `sigma_region` and none for `lambda_region` under
competition), so the data say little about how uneven the region is.

Leaving out total richness and beta diversity (total over mean richness)
changes little, so the choice does not rest on them alone:

``` r

print(res$reduced$chosen, digits = 2)
#>   observation        model probability votes_phylogenetic votes_none
#> 1           1 phylogenetic        0.99                489         11
print(res$reduced$posterior, digits = 2)
#>   observation     parameter estimate    low  high  mean
#> 1           1 log_dispersal   -0.998 -3.545 -0.61 -1.34
#> 2           1    log_kernel    0.017 -0.270  0.83  0.13
#> 3           1     inter_max    0.687  0.369  0.97  0.69
#> 4           1    sigma_mean   27.007  6.154 97.03 41.91
#> 5           1        lambda    0.660  0.059  0.97  0.61
#> 6           1         q_tol    0.328  0.016  0.81  0.36
```

## Cross-validation over space

The last check asks whether the model fitted to some ponds predicts the
others. With spatial blocks, the ponds left out are not neighbours of
those kept (Roberts et al. 2017). In simulations with a wrong model, in
which species tracked a hidden, spatially autocorrelated environment,
blocks rejected it only slightly more often than random folds: by about
two percentage points in a grid of 100 sites, and not at all in a grid
of 25 (`data-raw/abc_cv.R`). We therefore use blocks because they keep
the test from resting on neighbours, not because they add much power. We
split each city into two blocks with
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md),
so that each of the four folds leaves out half of one city, and
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md)
computed the statistics of the ponds kept and of those left out in every
simulation:

``` r

ggplot(data.frame(anurans$coords, fold = factor(res$folds), city = anurans$city),
       aes(lon, lat, colour = fold)) +
  geom_point() +
  facet_wrap(~city, scales = "free") +
  theme_minimal()
```

![Ponds of the two cities coloured by
fold.](abc_files/figure-html/folds-1.png)

[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md)
fits the forests on the statistics of the ponds kept and predicts those
of the ponds left out.
[`summary()`](https://rdrr.io/r/base/summary.html) gives, for each
statistic, the fraction of folds in which the observed value fell inside
the 90% predictive interval. These intervals are conservative (they hold
the observed value more often than 90% when the model is true), so a
statistic outside them counts against the model more than a statistic
inside them counts for it:

``` r

print(res$cv_summary, digits = 2)
#>          statistic n coverage z_mean z_rms rejected
#> 1            alpha 4     1.00   0.15  0.60     0.00
#> 2             beta 4     1.00   1.00  1.11     0.00
#> 3  composition_env 4     0.50  -1.20  1.28     0.50
#> 4            decay 4     0.50   0.17  2.47     0.50
#> 5            gamma 4     1.00   0.56  1.09     0.00
#> 6        occupancy 4     1.00  -0.84  0.92     0.00
#> 7     occupancy_sd 4     1.00  -0.13  0.36     0.00
#> 8     richness_env 4     0.75   0.58  1.14     0.00
#> 9         ses_mntd 4     1.00  -0.45  0.55     0.00
#> 10         ses_mpd 4     0.25  -2.57  2.63     0.75
```

With only four folds, these fractions are coarse, but most statistics of
the ponds left out are predicted. The exceptions repeat the posterior
predictive check: the standardised MPD of the ponds left out falls below
the prediction in most folds, so the misfit is not confined to a part of
either city. The link between composition and the environment is also
below the prediction in two folds. Distance decay misses in two folds,
once above and once below the prediction; it rests on few pairs of ponds
within half a city, so these misses may be noise, but they may also mean
that a single exponential kernel for both cities is too simple a
description of space.

## What this example shows

The workflow ran from a simulator to checks of the fit, and each step
said something the previous one could not. Recovery showed, before we
looked at the data, that `lambda` cannot be learned from them. Model
choice preferred competition among close relatives with a probability
near one, and the posterior predictive checks and cross-validation then
showed that this model does not reproduce the phylogenetic structure of
the ponds. Of the two hypotheses this suggested, a filter conserved in
the phylogeny did not remove the misfit in simulations, and limited
colonisation from the region reduced both the misfit and the estimate of
competition without removing either. An earlier version of this
analysis, with the pool of each city restricted to the species recorded
there, made the observed total richness the largest value a simulation
could reach, and that bound alone pushed the estimate of competition
towards zero; the choice of the regional pool matters as much as the
choice of the model. The analysis is meant to illustrate the functions,
not to test hypotheses about these frogs, which would also need a model
of detection if abundances were used (see
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md)).

## Planning a field study of trait variation

Incidence and abundances say little about the genetic parameters of the
evolutionary engine
([`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)),
so a study that wants to estimate them needs measured individuals. A
simulation study in
[`data-raw/abc_identifiability.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/abc_identifiability.R)
asked which statistics of individuals, and which sampling designs, make
the genetic variance G, the ratio of effective to census size Ne/N and
the environmental variance V_E recoverable by
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
with the same criterion as above (R² of at least 0.5 and coverage of the
90% interval of at least 0.85). It used eight species on a 10 x 10 grid
and 10,000 simulations per case, in eight cases: a static or a warming
environment; one G for every population, or G of each population set by
mutation, drift and selection (`mutation_var`); and ancestral optima the
same in every run, or drawn from priors as they would be unknown for
real species. The individuals of each design were summarised with
[`population_stats()`](https://provetelab.org/comdynamite/reference/population_stats.md),
which takes the output of
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
or field measurements:

``` r

ind <- sample_individuals(sim, n = 10, source = "evolution", step = c(350, 400))
population_stats(ind, abundance = "abundance", step = "step", coords = coords)
```

With the means of the populations alone (their regression on the
environment), G was well recovered (R² 0.89 to 0.92), but V_E (0.15 to
0.45) and Ne/N (0.08 to 0.33) were not. The variance within populations
is what separates G from V_E: with it, V_E reached R² 0.72 to 0.78 with
one G and 0.94 to 0.96 with mutation, and the mutational variance, when
it was the parameter, rose from 0.64 to 0.71 to 0.82 to 0.91. Ne/N
stayed the hardest parameter. With the QST-like divergence beyond the
environment (Spitze 1993), the semivariance of neighbouring populations,
the relation of the residuals to population size and, with two surveys,
the change of the means between them (the basis of the temporal method,
Waples 1989), its R² rose to 0.34 to 0.61, and it met the criterion in
four of the eight cases with the richest design; it fell short under
warming with mutation in particular. Unknown optima did not lower the
recovery of G or Ne/N; their centre was recovered with R² 0.40 to 0.76
and their spread with 0.23 to 0.35 (one G) or 0.67 to 0.73 (mutation).

Among the designs, the number of populations mattered more than the
number of individuals in each. With 10 individuals at 20 of the 100
sites instead of at every occupied site, the R² of Ne/N fell from 0.29
to 0.52 to 0.18 to 0.38, and that of the mutational variance by 0.07 to
0.13; with 5 or 30 individuals per population instead of 10, the R²
changed by at most 0.08. A second survey 50 steps before the last raised
the R² of Ne/N by 0.10 and 0.11 in the static environment with one G,
and by 0.03 or less in the other cases.

Two further studies qualify the role of a second survey. When the
parameters were G and the strength of drift, G / (Ne/N), instead of G
and Ne/N
([`data-raw/abc_g_ne.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/abc_g_ne.R)),
the strength of drift was recovered with R² 0.92 to 0.96 in every case,
but Ne/N computed from the two estimates was never better than its own
forest. Ten consecutive surveys, summarised by the change of population
means per generation (`temporal_step` in
[`population_stats()`](https://provetelab.org/comdynamite/reference/population_stats.md);
[`data-raw/abc_temporal_drift.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/abc_temporal_drift.R)),
raised the R² of Ne/N by 0.07 at most over two surveys 50 steps apart.
Since drift per generation is proportional to G / Ne (Lande 1976), trait
means, even repeated in time, seem to inform the strength of drift
rather than Ne itself, which would need independent information on G,
such as a heritability estimated in a common garden.

With 200 pseudo-observations, differences of about 0.05 in R² are close
to the sampling error, so these are tendencies rather than rules, and
they hold for this simulator and these priors. A field study can use
them as a starting point: measure about ten individuals in as many
populations as possible along the environmental gradient, survey more
than once if the strength of drift matters, and run
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
with its own simulator and design before going to the field.

## References

Lande, R. (1976). Natural selection and random genetic drift in
phenotypic evolution. *Evolution*, 30, 314–334.
<https://doi.org/10.1111/j.1558-5646.1976.tb00911.x>

Pudlo, P., Marin, J.-M., Estoup, A., Cornuet, J.-M., Gautier, M., &
Robert, C. P. (2016). Reliable ABC model choice via random forests.
*Bioinformatics*, 32, 859–866.
<https://doi.org/10.1093/bioinformatics/btv684>

Raynal, L., Marin, J.-M., Pudlo, P., Ribatet, M., Robert, C. P., &
Estoup, A. (2019). ABC random forests for Bayesian parameter inference.
*Bioinformatics*, 35, 1720–1728.
<https://doi.org/10.1093/bioinformatics/bty867>

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J.,
Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder,
B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann,
C. F. (2017). Cross-validation strategies for data with temporal,
spatial, hierarchical, or phylogenetic structure. *Ecography*, 40,
913–929. <https://doi.org/10.1111/ecog.02881>

Spitze, K. (1993). Population structure in *Daphnia obtusa*:
quantitative genetic and allozymic variation. *Genetics*, 135, 367–374.
<https://doi.org/10.1093/genetics/135.2.367>

Waples, R. S. (1989). A generalized approach for estimating effective
population size from temporal changes in allele frequency. *Genetics*,
121, 379–391. <https://doi.org/10.1093/genetics/121.2.379>
