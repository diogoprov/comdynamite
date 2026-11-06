# Changelog

## comdynamite (development version)

- [`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md):
  static metacommunities from unimodal species responses with Poisson
  counts, the Gaussian response of Peres-Neto et al. (2012), as a
  baseline without population dynamics. With several gradients, for
  example environment and space, the expected abundance is the product
  of the responses.
- [`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md)
  and
  [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md):
  presences and absences from logistic responses to two gradients, with
  the seven scenarios of Appendix C of Peres-Neto et al. (2012).

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
