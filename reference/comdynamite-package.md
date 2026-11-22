# comdynamite: Linking Traits, Phylogenies and Biogeography in Metacommunities

Methods that link traits, phylogenies, biogeography and metacommunities,
to analyse data and to simulate it. Analyses include the node regression
of Leibold et al. (2010)
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)
, with biogeographic events and restricted null models, the phylogenetic
gradient analysis of Peres-Neto et al. (2012)
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1) , and the
phylogenetic community composition and beta diversity tests of
Peres-Neto and Kembel (2015)
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)
. Simulations follow the process-based model of Thompson et al. (2020)
[doi:10.1111/ele.13568](https://doi.org/10.1111/ele.13568) , with
environmental filtering, competition, demographic stochasticity and
dispersal in spatially implicit or grid landscapes. Species traits can
carry phylogenetic signal or evolve under an Ornstein-Uhlenbeck process,
and competition can be structured by relatedness. Rule-based assembly
follows Kraft et al. (2007)
[doi:10.1086/519400](https://doi.org/10.1086/519400) . The node
regression and the phylogenetic gradient analysis are translations of
the MATLAB code in the supplements of Leibold et al. (2010) and
Peres-Neto et al. (2012); the simulation models were implemented from
the published equations.

## Analyse data

- Which clades respond to the environment or to biogeographic events?
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
  codes the nodes,
  [`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
  the events and
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  tests them.
  [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)
  and
  [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)
  show the results.

- Does phylogenetic composition follow the environment across the whole
  tree?
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  and
  [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  work on phylogenetic eigenvectors and
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  on the phylogenetic distances between communities, from
  [`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md),
  shown by
  [`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md).

- Which phylogenetic axis relates to a predictor?
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).

- How much of trait or phylogenetic variation is due to environment and
  space?
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  and
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md).

- Are sites neighbours?
  [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
  gives restricted null models for a torus, and
  [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
  and
  [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
  for the events of a landscape with edges.

## Simulate data

- Dynamics with growth, competition and dispersal:
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  from
  [`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md)
  or
  [`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)
  and
  [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md).

- Static communities along gradients:
  [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)
  and
  [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md).

- Assembly by rules, without dynamics:
  [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md).

- Traits on a tree:
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  [`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md)
  and
  [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md),
  with correlated traits from
  [`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md).

## Arguments

Arguments that choose among options, such as `method`, `null` or
`model`, are matched by
[`match.arg()`](https://rdrr.io/r/base/match.arg.html) and are
case-sensitive: `null = "variogram"` works and `null = "Variogram"`
stops with the list of valid options. Each scalar argument is checked
for its length, type and range before any computation, and missing or
infinite values in `comm`, `env` and the other data inputs stop the
function with a message that names the argument; remove or impute
missing values before the analysis. Factors used as labels (the site and
species columns of the data on individuals, `group` of
[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
`regimes` of the DAMOCLES functions, `strata` of
[`env_surrogates()`](https://provetelab.org/comdynamite/reference/env_surrogates.md))
are taken as unordered: the order of the levels of an ordered factor is
ignored. Factors and character columns are refused as predictors in
`env`; code them as dummy variables.

Inputs are converted once, at the start: community tables become double
matrices (integer counts included, so the `comm` of
[`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
is double), and site and species labels in results are character, even
when the input columns were factors. Attributes and classes of numeric
columns, such as units of the units package, are dropped without a
message: values are taken in the units given, and results are in those
units. The analyses take either an
[`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
object, which keeps the parts matched by site and species names, or the
parts as separate arguments (`comm`, `tree`, `traits`, `env`), which are
matched by the same names inside each function.

## Citation

Cite the package with `citation("comdynamite")` and, with it, the
original papers of the methods you use. Each function that implements a
published method shows these papers the first time it runs in a session,
and
[`method_references()`](https://provetelab.org/comdynamite/reference/method_references.md)
lists them by function; see
[`method_references()`](https://provetelab.org/comdynamite/reference/method_references.md)
for how to turn the message off.

## Audience

The package is written for community ecologists, biogeographers and
evolutionary biologists who analyse species composition, traits and
phylogenies across sites, and for methodologists who want to check, by
simulation, whether a method recovers a known process. It assumes
familiarity with R and with the methods it implements, whose original
papers each help page cites.

## Where to start

The articles at <https://provetelab.org/comdynamite/articles/> go
through each question. "Preparing your data" shows the format that the
analyses expect.

## See also

Useful links:

- <https://provetelab.org/comdynamite/>

- <https://github.com/diogoprov/comdynamite>

- Report bugs at <https://github.com/diogoprov/comdynamite/issues>

## Author

**Maintainer**: Diogo B. Provete <dbprovete@gmail.com>
([ORCID](https://orcid.org/0000-0002-0097-0651))

Authors:

- Diogo B. Provete <dbprovete@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-0097-0651))
