# comdynamite (development version)

- `trait_partition()` with `site_permutations`: each set of predictors now moves alone in its
  tests, with the response and the other set in place: `env` for `ab`, `a` and its
  semipartials, `space` for `bc`, `c` and its semipartials; `abc` still moves both. Up to
  0.6.1 the moves took both sets together in every test, which for `c` compared space with an
  environment that no longer explained anything: on a grid it rejected 22% of data sets with a
  strong environment without structure and never with a structured one. In
  `data-raw/validation_fractions_null.R` (choice on one generator, check on another), moving
  `env` alone rejected `a` in 3.7% to 7.1% of data sets where the environment had no effect, where
  free permutations reached 20% with an environment structured outside `space`, and moving `space`
  alone rejected `c` in 3.0% to 8.0% of null data sets, where free permutations reached 24% with
  a strong, structured environment. When the environment acts on a trait other than the one
  analysed, the test of sites of `a` still rejects, in 45% to 75%, as with free permutations:
  read `p_max`, which stayed at 2% to 4%.
  The power of that test of `c` was complete for a field in the span of `space` and 25% to 28%
  for a field that `space` holds only in part, since the maps moved around the grid explain
  such a field about as well: it asks whether the maps in `space` matter beyond the same maps
  elsewhere, not whether the remainder has any spatial structure. The result records
  `permutations$moved` and the print says which set moves. Free permutations are unchanged.
- Permutations of `pcc_test()`, `pcps_test()`, `phylo_beta_test()` and `trait_partition()` run in
  parallel under a plan of the future package, which is optional.
- Install from r-universe, as the README now shows.
- `summary_stats()`: the same statistics for simulated and observed data, covering diversity,
  occupancy, abundance distribution, SES MPD and MNTD, distance decay, trait-environment relation
  and, with more than one step, temporal variation.
- `pcps_test()`: principal coordinates of phylogenetic structure of Duarte et al. (2016), the
  method of the package PCPS, with tests of all the dissimilarity and of each axis, and
  `site_permutations` for a restricted null model, which that package does not have. Its
  statistics match those of PCPS 1.0.9 (`data-raw/validation_pcps.R`). With autocorrelated
  fields unrelated to each other, free permutations of sites rejected 29-80% at the 5% level;
  the torus rejected 4-7% when fields wrap around the grid and 6-25% when they stop at its edges.
- `torus_permutations()`: rectangular grids, with `side` of length two, and new `symmetries`
  and `n_perm`.

# comdynamite 0.6.1

- `trait_partition()`: the tests of fractions `a` and `c`, and of the semipartials, now use the
  share of the component that the other set leaves unexplained as statistic, which is monotone
  in the partial F, instead of the fraction itself, which the original code tests. The
  permutations are the same. The old test lost its power when the other set explained much:
  in `data-raw/validation_fractions.R` it rejected `a` in none of 500 data sets where the
  environment had no effect and space a strong one (new: 5.6% to 5.8%), and in 15% where the
  environment had a weak effect beside space (new: 100%). `fractions` has a new column
  `statistic`, and p-values are `NA` when a fraction cannot be estimated. See the validation
  article for the limits: a free permutation of sites rejects `a` too often when `env` has
  structure at scales outside `space`, and `site_permutations` are nearly degenerate for `c`.
- `trait_partition()`: test of the shared fraction `b` with `shared_null`, by Moran spectral
  randomisation of `env` from `coords` (adespatial) or by moves of `env` alone with
  `site_permutations`, as in `adespatial::envspace.test()`. `b` is tested only when `env` is
  related to `space` (`level`), and the result has `shared`. `plot_trait_partition()` colours
  `b` by its test when there is one.
- `as.data.frame()`, `summary()` and `plot()` methods for the results of the analyses.
- `rigid_permutations()`: restricted null model for the events of `node_test()` that moves the
  whole event by a rotation or a reflection and a translation, so it keeps the shape of the
  sides and the distance between them. It is better than `patch_permutations()` for sides with
  room to move, and gives few or no placements for events that span the landscape. Its
  `independent` argument covers those at the price of an assumption on the range of
  autocorrelation, and remains liberal in some layouts. See its help for the rates.
- `rigid_permutations()` off a regular lattice now enumerates the moves, on a grid of
  rotations and reflections with a translation to each site, instead of drawing them at
  random, so which placements exist no longer changes between runs. In
  `data-raw/validation_offgrid.R` it rejected 4.7% and 5.9% of the tests at the 5% level for
  adjacent areas off a lattice, and 6.8% to 8.0% for two blocks placed at random, where it is
  liberal.
- `rigid_permutations()` and `patch_permutations()`: the default `gap` follows the fourth
  nearest site, the step of a full grid and the spacing of sparse sites, instead of the
  nearest one. `rigid_permutations()` thus no longer stops on a grid with holes, and the holes
  between sites count as inside the landscape while the rim beyond an edge does not. The
  default of `patch_permutations()` is unchanged on a full grid where most sites are off the edge, and
  larger with sparse or irregular sites, where it changes the patches and the results: with
  25% to 40% of the sites of a square and areas in opposite corners it rejected 5.6% to 5.9%,
  against 15% to 17% before (`data-raw/validation_gap.R`).
- `rigid_permutations()` draws the placements at random when an event has more than `n_perm`
  of them. Before, with fewer than twice `n_perm`, it took the first ones in the order of
  enumeration.
- `node_test()`: `permutations$areas` now gives the number of permutations of each event, which
  may differ, and `print()` shows them.
- `plot_phylo_beta()`: community matrix with sites ordered by `phylo_beta_test()` and species
  by the tree, as Peres-Neto & Kembel (2015) suggest.
- `node_test()`: `r2_weighted`, which stays between 0 and 1, unlike `r2`, inherited from the
  original code; `n_draws`, the valid draws of sites. Errors in parallel workers now stop the
  function. The default of `n_perm` is 999, as in the other functions.
- `trait_partition()`: `site_permutations`.
- `pcc_test()`, `phylo_beta_test()`, `phylo_vector()` and `trait_partition()` accept empty sites
  with `site_permutations`, since these sites are part of the landscape.
- The significance level is now `level` in `phylo_vector()`, `plot_node_test()`,
  `plot_node_heatmap()` and `plot_trait_partition()`, because `alpha` also names the
  competition matrix. `alpha` still works, with a warning.
- Names of results, which break code that reads them: `p_taxa` is now `p_species`,
  `p_partial_taxa` is `p_partial_species`, and `p_partial`, the larger of two p-values, is
  `p_partial_max`, as in `p_max`. `r2` of `pcc_test()` is `r2_mean` and that of
  `phylo_beta_test()` is `r2_adjusted`, which say what they hold. In `node_test()` and
  `trait_partition()`, `r2` is unchanged. `statistic` stays, as the test statistic of each
  function.
- `phylo_vector()`: `predictor` is now `env`, as in the other functions, and
  `plot_trait_partition()`: `test = "taxa"` is now `"species"`. The old forms work, with a
  warning.
- `patch_permutations()`: evaluated with autocorrelation up to the scale of the landscape, and
  advice on `gap` for sparse sites. The rates given in 0.6.0 averaged over random places of
  areas. With fixed areas the test is liberal for sides far apart (6% to 10% at the 5% level)
  and conservative for neighbouring sides, and the help now says so.
- `simulate_metacommunity()` matches `alpha` to `traits` by name. Before, a matrix with
  species in another order gave wrong results without warning.
- `growth_rates()`, `simulate_static()` and `simulate_metacommunity()` accept traits as a matrix
  or a tibble, and species without names are called sp1, sp2, ... `simulate_static()` matches
  sets of traits by name. `assemble_community()` warns when it invents names.
- Predictors must vary among occupied sites with `site_permutations`.
- Help: weights of `node_test()`, and what the tests of fractions of `trait_partition()` ask.
- `pcc_test()` gives an `r2` of 0, not `NaN`, when no adjusted R2 is positive.
- `node_contrasts()` accepts clades with branches of length zero, and `plot_node_heatmap()`
  accepts a predictor and an event with the same name.
- Clearer errors in `simulate_metacommunity()`, `mc_dynamics()` with one site,
  `logistic_gradients()` and `plot_node_heatmap()` without tested nodes.

# comdynamite 0.6.0

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
