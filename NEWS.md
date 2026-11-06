# comdynamite (development version)

- `simulate_static()`: static metacommunities from the Gaussian response of Dray
  & Legendre (2008), also used by Peres-Neto et al. (2012), as a baseline without
  population dynamics. Data can be Poisson counts, binomial presences or the
  expected values.
  With several gradients, for example environment and space, the expected
  abundance is the product of the responses.
- `simulate_logistic()` and `logistic_gradients()`: presences and absences from
  logistic responses to two gradients, with the seven scenarios of Appendix C of
  Peres-Neto et al. (2012).
- `pcc_test()` and `phylo_eigenvectors()`: phylogenetic community composition of
  Peres-Neto & Kembel (2015), with site and species permutations.
- New article with animations of metacommunity dynamics.

# comdynamite 0.1.0

- First version: the model of Thompson et al. (2020) with phylogenetic traits
  and competition, spatially implicit and grid landscapes, torus permutations
  and time-varying environments. Code from the pcps-sim project.

- Rule-based assembly (`assemble_community()`): limiting similarity (GAUSE
  algorithm of Colwell & Winkler 1984), habitat filtering and neutral assembly,
  as in Kraft et al. (2007).
- `traits_ou()`: traits under an Ornstein-Uhlenbeck process with one or several
  optima (multi-peak OU; Brownian motion when `alpha = 0`), with
  `paint_regimes()` to assign regimes to clades.
- New vignettes: `vignette("phylogenetic-structure")`, on how dispersal, niche
  conservatism and competition shape phylogenetic community structure, and
  `vignette("evaluating-methods")`, on type I error and power with simulations.
- `vignette("evaluating-methods")` shows tree shape as a factor of the design,
  with trees from treeshapesim (in Suggests).
