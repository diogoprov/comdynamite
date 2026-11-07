# Package index

## Simulate metacommunities

The process-based model, from landscape and traits to abundances.

- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  : Simulate a metacommunity
- [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)
  : Simulate a static metacommunity
- [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md)
  : Simulate presences and absences from logistic responses
- [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md)
  : Gradients of the seven scenarios of Peres-Neto et al. (2012)
- [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md)
  : Metacommunity dynamics from growth rates
- [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md)
  : Growth rates of each species at each site

## Landscapes

Spatially implicit and grid landscapes, dispersal and environment.

- [`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md)
  : Spatially implicit landscape
- [`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md)
  : Spatially explicit landscape on a square grid
- [`dispersal_kernel()`](https://provetelab.org/comdynamite/reference/dispersal_kernel.md)
  : Exponential dispersal kernel
- [`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md)
  : Environment that changes over time

## Species traits

Traits with phylogenetic signal and evolutionary regimes.

- [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)
  : Species traits, with or without phylogenetic signal
- [`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md)
  : Trait evolution by Brownian motion with Pagel's lambda
- [`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md)
  : Traits evolving by an Ornstein-Uhlenbeck process with one or several
  optima
- [`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md)
  : Assign selective regimes to clades

## Competition

Interspecific competition, optionally structured by relatedness.

- [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md)
  : Competition matrix

## Clades and biogeography

Which clades respond to the environment or to biogeographic events?

- [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
  : Node contrasts of a phylogeny
- [`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
  : Matrix of biogeographic events
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  : Test clades against the environment, node by node
- [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)
  : Plot the results of a node test on the phylogeny
- [`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)
  : Plot the semipartial correlations of a node test

## The whole tree

Does phylogenetic composition follow the environment?

- [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  : Test phylogenetic community composition against the environment
- [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
  : Principal coordinates of phylogenetic structure
- [`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md)
  : Phylogenetic eigenvectors
- [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  : Test phylogenetic beta diversity against the environment
- [`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md)
  : Phylogenetic beta diversity between communities
- [`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md)
  : Community matrix ordered by phylogenetic beta diversity

## Environment and space

Phylogenetic vectors and the partition of trait variation.

- [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  : Represent a phylogeny by one vector related to a predictor
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
  : Partition trait or phylogenetic variation between environment and
  space
- [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md)
  : Plot the fractions of a trait partition

## Restricted null models

Permutations for sites that are neighbours.

- [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
  : Site permutations that keep the structure of a torus
- [`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
  : Permutations that move biogeographic events as rigid bodies
- [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
  : Permutations that move the patches of biogeographic events

## Results

Tables, summaries and plots of results.

- [`as.data.frame(`*`<mc_node_test>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`summary(`*`<mc_node_test>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`plot(`*`<mc_node_test>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`as.data.frame(`*`<mc_pcc>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`summary(`*`<mc_pcc>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`as.data.frame(`*`<mc_pcps>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`summary(`*`<mc_pcps>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`as.data.frame(`*`<mc_phylo_beta>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`summary(`*`<mc_phylo_beta>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`plot(`*`<mc_phylo_beta>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`as.data.frame(`*`<mc_phylo_vector>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`summary(`*`<mc_phylo_vector>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`as.data.frame(`*`<mc_trait_partition>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`summary(`*`<mc_trait_partition>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  [`plot(`*`<mc_trait_partition>`*`)`](https://provetelab.org/comdynamite/reference/result-methods.md)
  : Results as data frames, summaries and plots

## Data

- [`daphniids`](https://provetelab.org/comdynamite/reference/daphniids.md)
  : Daphniids of the northeastern United States

## Rule-based assembly

Community assembly without population dynamics.

- [`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)
  : Assemble a local community from a regional pool by rules
