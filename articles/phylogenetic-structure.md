# Dispersal, niche conservatism and phylogenetic community structure

This vignette uses comdynamite, which implements the model of Thompson
et al. (2020), to ask an ecological question: **how do dispersal,
phylogenetic signal in niches and competition among relatives combine to
shape the phylogenetic structure of local communities?** Phylogenetic
clustering is often read as evidence of environmental filtering on
conserved niches, and overdispersion as evidence of competition (Webb et
al. 2002). Because a simulation fixes the processes, it shows how far
those readings hold once dispersal also operates.

The simulations take a few minutes, so the results were computed in
advance and are loaded below. The complete script is
[`data-raw/vignette_ecology.R`](https://github.com/diogoprov/comdynamite/blob/main/data-raw/vignette_ecology.R).

## Design

We cross three factors:

- **dispersal**, the emigration probability, from 0.001 to 0.5;
- **niche conservatism**: optima without phylogenetic signal
  (`lambda = 0`) or evolving under Brownian motion (`lambda = 1`);
- **competition**: none among species, or stronger between close
  relatives.

Each replicate draws one tree of 30 species and one landscape of 25
sites, and uses them in all scenarios, so that scenarios are compared on
the same regional pool. Traits are drawn once per value of `lambda` and
shared across dispersal rates and competition types.

``` r

library(comdynamite)

tree <- ape::rcoal(30)
land <- landscape_implicit(25)
traits <- species_traits(tree, lambda = 1, sigma_mean = 8, sigma_sd = 2)
alpha <- competition_matrix(tree, "phylo", inter_max = 0.5)
sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.01, n_steps = 1000)
```

For each simulation we record the mean local richness (alpha), the
number of species in the metacommunity (gamma), the mean pairwise
Jaccard dissimilarity among sites (beta), and the standardised effect
size of the mean pairwise distance (SES MPD; computed with picante,
Kembel et al. 2010), averaged over sites. SES MPD is negative under
phylogenetic clustering and positive under overdispersion.

``` r

pa <- (sim$abundance > 0) * 1
pa <- pa[rowSums(pa) > 1, , drop = FALSE]
z <- picante::ses.mpd(pa, cophenetic(tree), null.model = "taxa.labels",
                      runs = 99)$mpd.obs.z
mean(z[is.finite(z)])
```

The null model draws from all the species in the tree, which is the
regional pool. Using only the species that persisted would remove the
effect of regional extinction from the comparison.

## Results

``` r

results <- readRDS("phylogenetic-structure.rds")
means <- aggregate(cbind(alpha, gamma, beta, ses_mpd) ~ dispersal + lambda + competition,
                   data = results, FUN = mean)
knitr::kable(means, digits = 2)
```

| dispersal | lambda | competition | alpha | gamma | beta | ses_mpd |
|----------:|-------:|:------------|------:|------:|-----:|--------:|
|      0.00 |      0 | none        |  7.43 | 29.37 | 0.80 |    0.06 |
|      0.01 |      0 | none        | 10.08 | 29.43 | 0.74 |    0.05 |
|      0.05 |      0 | none        | 17.08 | 29.37 | 0.51 |   -0.06 |
|      0.10 |      0 | none        | 20.95 | 29.23 | 0.36 |    0.01 |
|      0.50 |      0 | none        | 24.01 | 25.07 | 0.07 |   -0.07 |
|      0.00 |      1 | none        |  7.33 | 29.47 | 0.80 |   -1.65 |
|      0.01 |      1 | none        |  9.93 | 29.43 | 0.74 |   -0.97 |
|      0.05 |      1 | none        | 16.81 | 29.43 | 0.51 |   -0.41 |
|      0.10 |      1 | none        | 20.78 | 29.30 | 0.36 |   -0.25 |
|      0.50 |      1 | none        | 23.73 | 24.83 | 0.07 |    0.00 |
|      0.00 |      0 | phylo       |  3.35 | 18.80 | 0.85 |    0.61 |
|      0.01 |      0 | phylo       |  4.71 | 18.83 | 0.80 |    0.49 |
|      0.05 |      0 | phylo       |  8.71 | 18.60 | 0.61 |    0.52 |
|      0.10 |      0 | phylo       | 11.08 | 17.93 | 0.46 |    0.58 |
|      0.50 |      0 | phylo       | 10.79 | 11.37 | 0.08 |    0.72 |
|      0.00 |      1 | phylo       |  3.01 | 16.23 | 0.84 |   -0.29 |
|      0.01 |      1 | phylo       |  4.13 | 16.63 | 0.80 |   -0.15 |
|      0.05 |      1 | phylo       |  7.46 | 16.00 | 0.61 |    0.29 |
|      0.10 |      1 | phylo       |  9.73 | 15.63 | 0.46 |    0.41 |
|      0.50 |      1 | phylo       | 10.03 | 10.60 | 0.08 |    0.79 |

``` r

groups <- unique(means[c("lambda", "competition")])
cols <- c("#E69F00", "#0072B2", "#009E73", "#D55E00")
old <- par(mfrow = c(1, 3), mar = c(4, 4, 1, 1), mgp = c(2.4, 0.7, 0))
for (v in c("ses_mpd", "alpha", "beta")) {
  plot(means$dispersal, means[[v]], log = "x", type = "n", xlab = "dispersal",
       ylab = c(ses_mpd = "SES MPD", alpha = "local richness", beta = "beta (Jaccard)")[v])
  if (v == "ses_mpd") abline(h = 0, col = "grey70")
  for (g in seq_len(nrow(groups))) {
    d <- means[means$lambda == groups$lambda[g] &
                 means$competition == groups$competition[g], ]
    lines(d$dispersal, d[[v]], col = cols[g], lwd = 2, type = "b", pch = 15 + g)
  }
}
legend("bottomleft", paste0("lambda = ", groups$lambda, ", ", groups$competition),
       col = cols, pch = 15 + seq_len(nrow(groups)), lwd = 2, bty = "n", cex = 0.8)
```

![SES MPD, local richness and beta diversity against dispersal, for the
four combinations of phylogenetic signal and
competition.](phylogenetic-structure_files/figure-html/unnamed-chunk-5-1.png)

``` r

par(old)
```

Four patterns stand out.

**Clustering needs both conserved niches and limited dispersal.**
Without competition and with `lambda = 1`, communities are clustered at
low dispersal (SES MPD near -1.7) and the signal fades as dispersal
increases, reaching zero at 0.5. Niche conservatism alone does not
produce clustering: immigrants keep species in sites that do not suit
them, so composition no longer reflects the match between niches and
environment.

**Without phylogenetic signal, filtering leaves no phylogenetic trace.**
With `lambda = 0` and no competition, SES MPD stays near zero at all
dispersal rates, even though species are strongly sorted by the
environment at low dispersal (beta near 0.8).

**Competition among relatives can hide filtering.** With `lambda = 1`,
competition among relatives moves SES MPD from about -1.7 to about -0.3
at the lowest dispersal rate, which is close to the random expectation
although both processes are acting strongly. At higher dispersal,
communities become overdispersed. A value near zero therefore does not
imply that assembly is random.

**Taxonomic beta diversity does not see the phylogenetic signal.** Beta
diversity depends on dispersal and competition but is the same for
`lambda = 0` and `lambda = 1`. In this model, the signal changes which
relatives co-occur, not how different sites are in composition.

Competition also reduces the size of the metacommunity: it halves local
richness and, at the highest dispersal rate, gamma falls to about 11 of
the 30 species. A likely explanation, which these simulations do not
test directly, is that high dispersal exposes all species to the same
competitors in every site.

## Limitations

These results hold for the values used here. They come from one tree
size, one landscape size and one strength of competition, and the
magnitudes will change with those choices. The means are also averages
over 30 replicates with a standard error of up to 0.2 for SES MPD, so
small differences between scenarios should not be interpreted.

With competition among relatives and high dispersal, 1000 steps is not
an equilibrium: in a test run, regional exclusion continued slowly up to
3000 steps. Without competition, richness was stable by 300 steps.
Before drawing conclusions from a new parameter set, check that the
statistics of interest have stabilised, for example with `record_steps`
in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

SES MPD here uses presence and absence, so a single immigrant counts as
much as a large population. Sink populations therefore have their full
weight, and an abundance-weighted statistic may show a weaker effect of
dispersal.

## References

Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H.,
Ackerly, D. D., Blomberg, S. P., & Webb, C. O. (2010). Picante: R tools
for integrating phylogenies and ecology. *Bioinformatics*, 26,
1463–1464. <https://doi.org/10.1093/bioinformatics/btq166>

Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik,
R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A
process-based metacommunity framework linking local and regional scale
community ecology. *Ecology Letters*, 23, 1314–1329.
<https://doi.org/10.1111/ele.13568>

Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475–505.
<https://doi.org/10.1146/annurev.ecolsys.33.010802.150448>
