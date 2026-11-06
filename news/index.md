# Changelog

## comdynamite 0.5.0

- [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md):
  static metacommunities from the Gaussian response of Dray & Legendre
  (2008), also used by Peres-Neto et al. (2012), as a baseline without
  population dynamics. Data can be Poisson counts, binomial presences or
  the expected values. With several gradients, for example environment
  and space, the expected abundance is the product of the responses.
- [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md)
  and
  [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md):
  presences and absences from logistic responses to two gradients, with
  the seven scenarios of Appendix C of Peres-Neto et al. (2012).
- [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  and
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  phylogenetic gradient analysis of Peres-Neto et al. (2012), tested
  against the original code in Octave. The test between a trait and the
  environment stays in
  [`ade4::fourthcorner()`](https://adeverse.github.io/ade4/reference/fourthcorner.html).
- Functions that take a community and a phylogeny, traits or contrasts
  now match species by name and, as
  [`picante::match.phylo.comm()`](https://rdrr.io/pkg/picante/man/data.checking.html),
  drop with a warning the species found in only one of them. Repeated
  names are an error.
- [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md):
  semipartial correlations of each predictor, with their tests,
  following Appendix A of Peres-Neto et al. (2012).
- [`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md):
  results of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  on the tree, as in Figure 3 of Leibold et al. (2010), as a phylogram
  or a fan.
  [`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md):
  fractions of
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md).
- `daphniids`: data of Leibold et al. (2010), rebuilt from the public
  survey.
- [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
  now keeps nodes that join two tips, as the article does.
  `cherries = FALSE` gives the original code.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md):
  `semipartial` replaces `partial_r2`, following the article, and there
  are `p_adjusted` (false discovery rate), biogeographic events, the
  joint permutation of the article (`p_joint`), 9999 permutations by
  default and `cores`.
- [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
  and
  [`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md):
  phylogenetic community composition of Peres-Neto & Kembel (2015), with
  site and species permutations.
- [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
  and
  [`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md):
  phylogenetic beta-diversity framework of the same article, a
  distance-based redundancy analysis of the mean phylogenetic distance
  between communities.
- [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  and
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md):
  regression of clades on the environment, node by node, from Leibold et
  al. (2010), with site and species permutations.
- [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  and
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md):
  stochastic extirpation of populations (`extirpation`), used by
  Thompson et al. (2020) in their competition-colonisation scenario.
- New article on the validation of the model, with the differences
  between the article of Thompson et al. (2020) and the code archived
  with it.
- New article with animations of metacommunity dynamics.

## comdynamite 0.1.0

- First version: the model of Thompson et al. (2020) with phylogenetic
  traits and competition, spatially implicit and grid landscapes, torus
  permutations and time-varying environments. Code from the pcps-sim
  project.

- Rule-based assembly
  ([`assemble_community()`](https://provetelab.org/comdynamite/reference/assemble_community.md)):
  limiting similarity (GAUSE algorithm of Colwell & Winkler 1984),
  habitat filtering and neutral assembly, as in Kraft et al. (2007).

- [`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md):
  traits under an Ornstein-Uhlenbeck process with one or several optima
  (multi-peak OU; Brownian motion when `alpha = 0`), with
  [`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md)
  to assign regimes to clades.

- New vignettes:
  [`vignette("phylogenetic-structure")`](https://provetelab.org/comdynamite/articles/phylogenetic-structure.md),
  on how dispersal, niche conservatism and competition shape
  phylogenetic community structure, and
  [`vignette("evaluating-methods")`](https://provetelab.org/comdynamite/articles/evaluating-methods.md),
  on type I error and power with simulations.

- [`vignette("evaluating-methods")`](https://provetelab.org/comdynamite/articles/evaluating-methods.md)
  shows tree shape as a factor of the design, with trees from
  treeshapesim (in Suggests).
