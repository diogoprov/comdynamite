# comdynamite

Methods that link traits, phylogenies, biogeography and metacommunities
in R, to analyse data and to simulate it.

The simulations work at the level of populations in ecological time:
abundances in every site grow, compete and disperse, and, when asked,
each population carries the genetic means of one or more traits that
evolve by selection and drift under quantitative genetics (Lande 1979),
with a covariance matrix *G* per species. Lineages can split into new
species inside the metacommunity, so the same run gives communities,
trait distributions of individuals and the true phylogeny. Because the
truth is known, the analysis functions of the package can be tested on
these simulations before they are applied to field data.

## What do you want to do?

**Analyse data**

| Question | Functions | Article |
|----|----|----|
| Which clades respond to the environment or to biogeographic events? | [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md), [`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md), [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md) | [Analysing clades](https://provetelab.org/comdynamite/articles/analysing-clades.html) |
| Does phylogenetic composition follow the environment across the whole tree? | [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md), [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md), [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md) | [Testing the whole tree](https://provetelab.org/comdynamite/articles/whole-tree.html) |
| Which phylogenetic axis relates to a predictor? | [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md) | [Analysing clades](https://provetelab.org/comdynamite/articles/analysing-clades.html) |
| How much of trait or phylogenetic variation is due to environment and space? | [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md), [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md) | [Analysing clades](https://provetelab.org/comdynamite/articles/analysing-clades.html) |
| My sites are neighbours: which null model? | [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md), [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md), [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md) | [Analysing clades](https://provetelab.org/comdynamite/articles/analysing-clades.html) |
| How do I format my own data? |  | [Preparing your data](https://provetelab.org/comdynamite/articles/preparing-data.html) |

**Simulate data**

| Question | Functions | Article |
|----|----|----|
| Dynamics with growth, competition and dispersal, as in Thompson et al. (2020) | [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md) | [Get started](https://provetelab.org/comdynamite/articles/comdynamite.html) |
| Static communities along gradients | [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md), [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md) | [Testing the whole tree](https://provetelab.org/comdynamite/articles/whole-tree.html) |
| Assembly by rules, without dynamics, as in Kraft et al. (2007) | [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md) | [Watching assembly](https://provetelab.org/comdynamite/articles/watching-assembly.html) |
| Traits with phylogenetic signal or under Ornstein-Uhlenbeck regimes | [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md), [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md) | [Get started](https://provetelab.org/comdynamite/articles/comdynamite.html) |
| Local adaptation, migration load and evolutionary rescue, with one or several traits | [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md), [`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md), [`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md) | [Quantitative genetics in the metacommunity](https://provetelab.org/comdynamite/articles/quantitative-genetics.html) |
| Speciation inside the metacommunity, with the full and reconstructed trees | [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md), [`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md) |  |
| Individuals sampled from simulated populations, and tests of drift | [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md), [`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md) |  |
| Does a method work? Type I error and power by simulation | all of the above | [Evaluating methods](https://provetelab.org/comdynamite/articles/evaluating-methods.html) |

In R,
[`?comdynamite`](https://provetelab.org/comdynamite/reference/comdynamite-package.md)
gives the same map.

## Installation

``` r

install.packages("comdynamite", repos = c("https://diogoprov.r-universe.dev",
                                          "https://cloud.r-project.org"))
# or
remotes::install_github("diogoprov/comdynamite")
```

## Large analyses

The permutations of
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
and
[`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
run in parallel when you choose a plan with the future package, and
results can be reproduced with
[`set.seed()`](https://rdrr.io/r/base/Random.html) for any number of
workers.
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
has its own `cores` argument: it took about 4 minutes with 600 species,
100 sites and 999 permutations on one core.

``` r

future::plan(future::multisession, workers = 4)
```

Inputs can be tibbles of numeric columns, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) or
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
turn results into tables.

## Example

``` r

library(comdynamite)
set.seed(1)
tr <- treeshapesim::sim_tree_shape(40)
sim <- simulate_metacommunity(
  landscape  = landscape_implicit(50),
  traits     = species_traits(tr, lambda = 1),
  alpha      = competition_matrix(tr, type = "phylo"),
  dispersal  = 0.01,
  n_steps    = 1000
)
sim$abundance  # sites x species
```

The [website](https://provetelab.org/comdynamite/) has the articles and
the function reference. To read the vignettes in R, install with
`build_vignettes = TRUE`.

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

Provete, D. B. (2026). comdynamite: Linking traits, phylogenies and
biogeography in metacommunities. R package.
<https://doi.org/10.5281/zenodo.23993533>

This DOI always points to the latest version; each release also has its
own DOI on Zenodo. In R, `citation("comdynamite")` gives the entry with
the installed version.

## References

Lande, R. (1979). Quantitative genetic analysis of multivariate
evolution, applied to brain:body size allometry. *Evolution*, 33,
402–416. <https://doi.org/10.1111/j.1558-5646.1979.tb04694.x>

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
