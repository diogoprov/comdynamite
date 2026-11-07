# comdynamite (development version)

- `patch_permutations()`: restricted null model for the events of `node_test()` in landscapes
  with edges, which moves the patches of each area and keeps their number and size. It rejected
  3.0% to 6.3% of tests at the 5% level in `data-raw/validation_edges.R`, against 44% to 55%
  with free permutations.
- `node_test()`: `area_permutations` takes one matrix of permutations per event.
- New title and description, since the package now covers analyses as well as simulations.
- `?comdynamite` and the README map questions to functions, and the reference and articles on
  the website are grouped by question.
- New articles: "Preparing your data" and "Testing the whole tree", on `pcc_test()` and
  `phylo_beta_test()`.

# comdynamite 0.5.2

- `biogeography_matrix()`: events of `node_test()` from the area of each site, as pairs of
  areas stated in advance or as each area against the others.
- `node_test()`: `site_permutations` gives a restricted null model. For events, it moves the
  areas across the sites.
- `plot_node_test()`: biogeographic events are triangles, as in Leibold & Chase (2018).
- `plot_trait_partition()`: two overlapping bars along the total variation, as in Figure 1 of
  Desdevises et al. (2003), so that the shared fraction is seen as shared. New `labels`.
- Progress bars with an estimate of the time left in `node_test()`, `pcc_test()`,
  `phylo_beta_test()`, `phylo_vector()` and `trait_partition()`, in interactive sessions.
- Examples, README and vignettes use trees of treeshapesim, in Suggests. The results of
  `vignette("phylogenetic-structure")` and `vignette("evaluating-methods")` were run again.
- The help of functions adapted from the MATLAB code of the supplements says so.

Corrections that change results:

- `node_test()`: draws are now ranked on F statistics. The statistic of each predictor, taken
  from the original code, depended on the scale of the response, and neither statistic took
  into account that draws of species occupy different numbers of sites. `p_partial`,
  `p_partial_joint`, `p_taxa`, `p_max` and `p_joint` differ from those of 0.5.1, and so do the
  same p-values of events. `statistic`, `p_site` and `p_partial_site` of predictors do not.
- `node_test()`: a node with a clade without individuals is no longer tested, as it was a
  repeat of another node.
- `node_test()`: with `cores` above 1, each call now uses its own random streams.
- `pcc_test()` and `phylo_beta_test()`: species without individuals leave the analysis. In
  `pcc_test()` they could turn numerical noise into an axis.
- `trait_partition()`: species without individuals leave the analysis, which could stop it.

Checks of data:

- All analysis functions match sites by name when the data have row names, and stop if the
  names differ. Before, misordered predictors gave a wrong result without warning.
- Missing values, predictors without variance and columns that are not numeric now stop with
  a message that names them. `phylo_vector()` gave p = 1 in these cases.
- `node_test()` checks `biogeography` before the permutations and explains a tree in `nodes`.
- Results keep the numbers of permutations and the null model, and `print()` shows them
  with the smallest possible p-value. Nodes that were not tested are reported, and
  `phylo_vector()` says when no vector was retained.

# comdynamite 0.5.1

- `node_contrasts()`: `branch_lengths = TRUE` gives the coding with branch lengths of the
  supplement of Leibold et al. (2010).
- `plot_node_test()`: with `type = "fan"`, the names of the variables follow the radii and
  no longer overlap the names of the species.
- `plot_node_heatmap()`: semipartial correlations of `node_test()` as a table of nodes by
  variables.
- `simulate_logistic()` now uses the values of the code of Ives & Helmus (2011): the default
  of `intercept_sd` is 0.5, and slopes have a standard deviation of `slope_sd` among species.
- `logistic_gradients()`: in scenario 3, the error is added to `y` with unit variance and only
  in `y_test`, and in scenarios 5 and 6 `x` is in random order along the transect.
- New article on analysing clades and traits. The article on validation now compares Table 5
  of Peres-Neto et al. (2012).

# comdynamite 0.5.0

- `simulate_static()`: static metacommunities from the Gaussian response of Dray
  & Legendre (2008), also used by Peres-Neto et al. (2012), as a baseline without
  population dynamics. Data can be Poisson counts, binomial presences or the
  expected values.
  With several gradients, for example environment and space, the expected
  abundance is the product of the responses.
- `simulate_logistic()` and `logistic_gradients()`: presences and absences from
  logistic responses to two gradients, with the seven scenarios of Appendix C of
  Peres-Neto et al. (2012).
- `phylo_vector()` and `trait_partition()`: phylogenetic gradient analysis of Peres-Neto et al.
  (2012), tested against the original code in Octave.
  The test between a trait and the environment stays in `ade4::fourthcorner()`.
- Functions that take a community and a phylogeny, traits or contrasts now match
  species by name and, as `picante::match.phylo.comm()`, drop with a warning the
  species found in only one of them. Repeated names are an error.
- `trait_partition()`: semipartial correlations of each predictor, with their tests,
  following Appendix A of Peres-Neto et al. (2012).
- `plot_node_test()`: results of `node_test()` on the tree, as in Figure 3 of Leibold
  et al. (2010), as a phylogram or a fan. `plot_trait_partition()`: fractions of
  `trait_partition()`.
- `daphniids`: data of Leibold et al. (2010), rebuilt from the public survey.
- `node_contrasts()` now keeps nodes that join two tips, as the article does.
  `cherries = FALSE` gives the original code.
- `node_test()`: `semipartial` replaces `partial_r2`, following the article, and
  there are `p_adjusted` (false discovery rate), biogeographic events, the joint
  permutation of the article (`p_joint`), 9999 permutations by default and `cores`.
- `pcc_test()` and `phylo_eigenvectors()`: phylogenetic community composition of
  Peres-Neto & Kembel (2015), with site and species permutations.
- `phylo_beta_test()` and `phylo_beta()`: phylogenetic beta-diversity framework of
  the same article, a distance-based redundancy analysis of the mean phylogenetic
  distance between communities.
- `node_test()` and `node_contrasts()`: regression of clades on the environment,
  node by node, from Leibold et al. (2010), with site and species permutations.
- `simulate_metacommunity()` and `mc_dynamics()`: stochastic extirpation of
  populations (`extirpation`), used by Thompson et al. (2020) in their
  competition-colonisation scenario.
- New article on the validation of the model, with the differences between the
  article of Thompson et al. (2020) and the code archived with it.
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
