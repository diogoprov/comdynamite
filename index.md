# comdynamite

Process-based metacommunity simulations in R, with phylogenetic
structure.

comdynamite implements the model of Thompson et al. (2020):
density-independent growth set by the match between each species’
environmental optimum and the local environment, Beverton–Holt
competition, Poisson demographic stochasticity, and dispersal among
sites. On top of the original model it offers:

- **phylogenetic structure**: optima, niche breadths and maximum growth
  rates with phylogenetic signal (Pagel’s lambda), traits under single-
  or multi-optimum Ornstein-Uhlenbeck processes, and competition that is
  stronger between close relatives;
- **spatially implicit or explicit landscapes**: uniform dispersal among
  exchangeable sites, or a grid (optionally a torus) with an exponential
  dispersal kernel and a spatially autocorrelated environment;
- **restricted null models** for grids: site permutations of a torus
  that keep its spatial structure;
- an environment that can change over time, and abundances recorded at
  chosen steps;
- **rule-based assembly** without population dynamics, as in Kraft et
  al. (2007): limiting similarity (the GAUSE algorithm of Colwell &
  Winkler 1984), habitat filtering and neutral assembly, on one or
  several traits.

## Installation

``` r

remotes::install_github("diogoprov/comdynamite")
```

## Example

``` r

library(comdynamite)
set.seed(1)
tr <- ape::rcoal(40)
sim <- simulate_metacommunity(
  landscape  = landscape_implicit(50),
  traits     = species_traits(tr, lambda = 1),
  alpha      = competition_matrix(tr, type = "phylo"),
  dispersal  = 0.01,
  n_steps    = 1000
)
sim$abundance  # sites x species
```

See
[`vignette("comdynamite")`](https://provetelab.org/comdynamite/articles/comdynamite.md)
for species sorting vs. mass effects, grids and time-varying
environments,
[`vignette("phylogenetic-structure")`](https://provetelab.org/comdynamite/articles/phylogenetic-structure.md)
for an ecological question worked through simulations, and
[`vignette("evaluating-methods")`](https://provetelab.org/comdynamite/articles/evaluating-methods.md)
for the evaluation of a statistical method. The
[website](https://provetelab.org/comdynamite/) has the articles and the
function reference.

## Related packages

- **[mcomsimr](https://github.com/plthompson/mcomsimr)**, by Patrick
  Thompson, is the original R implementation of the model of Thompson et
  al. (2020). comdynamite was written independently from the published
  equations. It differs in its conventions (`sigma` is the standard
  deviation of the Gaussian response, so it equals `sqrt(2)` times the
  niche breadth of Thompson et al.; `alpha[j, k]` is the effect of
  species k on species j) and adds phylogenetic traits and competition,
  grid landscapes on a torus with their restricted permutations, and
  abundances recorded over time.

- **[gen3sis](https://github.com/project-gen3sis/R-package)** (Hagen et
  al.

  2021. is a general engine for eco-evolutionary simulations on
        landscapes that change over geological time. Speciation,
        extinction, dispersal and trait evolution are simulated, so the
        **phylogeny is an outcome** of the simulation, and the
        ecological rules are functions that the user writes. In
        comdynamite the **phylogeny is an input**: the tree is fixed,
        traits are placed on it with a chosen phylogenetic signal, and a
        fixed population model runs over ecological time (no
        speciation). This is what simulation studies of methods need,
        because tree shape and phylogenetic signal can be set
        independently and replicates can be paired by tree across
        scenarios. The two can be combined: a phylogeny produced by
        gen3sis can be used as the regional pool of a comdynamite
        simulation.

- **[treeshapesim](https://github.com/diogoprov/treeshapesim)**
  simulates trees with separately controlled imbalance and branching
  times, to use as regional pools in comdynamite. Tree shape can then be
  a factor in the design of a simulation study:

  ``` r

  tr <- treeshapesim::sim_tree_shape(40, beta = -1, c = -0.5)
  sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tr, lambda = 1),
                                competition_matrix(tr, type = "phylo"),
                                dispersal = 0.01, n_steps = 1000)
  ```

## Citation

Provete, D. B. (2026). comdynamite: Process-based metacommunity
simulations with phylogenetic structure. R package.
<https://doi.org/10.5281/zenodo.23993533>

This DOI always points to the latest version; each release also has its
own DOI on Zenodo. In R, `citation("comdynamite")` gives the entry with
the installed version.

## References

Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik,
R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A
process-based metacommunity framework linking local and regional scale
community ecology. *Ecology Letters*, 23, 1314–1329.
<https://doi.org/10.1111/ele.13568>

Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M.,
Rangel, T. F., & Pellissier, L. (2021). gen3sis: A general engine for
eco-evolutionary simulations of the processes that shape Earth’s
biodiversity. *PLoS Biology*, 19, e3001340.
<https://doi.org/10.1371/journal.pbio.3001340>

Colwell, R. K., & Winkler, D. W. (1984). A null model for null models in
biogeography. In D. R. Strong, D. Simberloff, L. G. Abele, & A. B.
Thistle (Eds.), *Ecological Communities: Conceptual Issues and the
Evidence* (pp. 344–359). Princeton University Press.
<https://doi.org/10.1515/9781400857081.344>

Kraft, N. J. B., Cornwell, W. K., Webb, C. O., & Ackerly, D. D. (2007).
Trait evolution, community assembly, and the phylogenetic structure of
ecological communities. *The American Naturalist*, 170, 271–283.
<https://doi.org/10.1086/519400>
