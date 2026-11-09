# comdynamite (development version)

- Plots for several traits, in base graphics: `plot_traits()` draws traits from `traits_ou()`
  with their tree, as a phylomorphospace (two traits) or a traitgram (one trait), with branches
  coloured by regime and the colour changing where a mapped history changes;
  `plot()` for `rate_matrix()` draws each regime's rate matrix as an ellipse; `plot_assembly()`
  draws a community from `assemble_community()` in trait space, with the boundary of the filter
  or the contours of limiting similarity in the metric of the assembly. New article,
  "Correlated traits and assembly in trait space".
- `traits_ou()`: `theta` can be a matrix with one column per trait, so that a regime shifts
  each trait's optimum by its own amount; `root` can have one value per trait and is checked;
  a mapped history is matched to the tree by its edges, so it may come in another order of
  edges. A named `theta` or `sigma2` of length one is no longer recycled to every regime: it
  must name every regime, as before rate matrices were added.
- `rate_matrix()` returns an object of class `rate_matrix`, with `print()` and `plot()`
  methods. `assemble_community()` returns an object of class `mc_assembly`, which keeps the
  rule, distance, optimum, traits and the matrix of the metric, with a `print()` method; it
  stops on missing or non-numeric traits instead of keeping such species or failing later.
- The limiting-similarity rule of `assemble_community()` keeps the nearest neighbour of each
  species up to date instead of scanning the whole distance matrix at each removal: 23 times
  faster for a pool of 2000 species reduced to 100 (from 10.6 to 0.46 s), with the same species
  removed and the same random draws.
- Simulations with a niche of several axes keep their tree. `as_mc_data()` keeps one column of
  optima per axis, `summary_stats()` uses the first axis, and `sample_individuals()` stops with
  a clear message; both use the `env` that set growth when `simulate_metacommunity(env = )` was
  given one. Axes of `env` without names take the names of the list of traits, or `axis1`,
  `axis2`, ...

- Plots of null-model tests: `plot()` for `tstats()` and `community_index()` draws, site by
  site, the null interval (hatched by default) with the observed value coloured by the test, or
  the standardised effect sizes; `null_intervals()` returns the table behind it.
- `trait_densities()` gives trait probability densities of populations and communities
  (Carmona et al. 2016, 2019) with one bandwidth for all populations, and
  `plot_population_traits()` draws them with adjustable transparency.
- Classic metrics of trait spacing and spread, for use with `community_index()`: `mnnd()`,
  `sdnnd()`, `cvnnd()`, `min_nnd()`, `max_nnd()`, `mnd()`, `sdnd()`, `cvnd()`, `mst_length()`,
  `mst_ratio()`, `trait_range()`, `trait_kurtosis()` and `nn_distances()`. Zero distances count
  (cati discards them), and one-trait metrics use sorted values instead of distance matrices.
- `community_index()` accepts a named list of indices, computed on the same randomised
  communities, and returns them in one table with an `index` column.
- Sampling effort: `sample_individuals()` gains `cover` (measure only the most abundant species
  of each site, as in Pakeman & Quested 2007) and `min_n`; `subsample_individuals()` applies
  the same schemes to real data, and `trait_rarefaction()` follows any statistic over the
  number of individuals per population or over `cover`, with a `plot()` method that marks the
  80% rule of Pakeman & Quested (2007) with a dashed line.
- `tstats()`: the null model of `T_IC.IR` now keeps the observed regional variance as the
  denominator. It was computed from the values drawn, which repeat individuals across sites;
  this biased the null upward and made the tests anti-conservative (about 6% of sites
  rejected per tail at a nominal 4% in simulations with random traits).
- New functions for analyses of intraspecific trait variation, written from the published
  definitions to replace functions of cati 0.99.6 that give wrong results (an audit found
  bugs reported as adrientaudiere/cati#13 to #21): `tstats_summary()` summarises the tests of
  `tstats()` (shares of sites below and above the null, test of the mean over sites, the sites
  that deviate), with sites that cannot be tested left out instead of given small p-values;
  `trait_flex()` splits community-weighted means into species turnover, intraspecific
  variability and covariation (Lepš et al. 2011), with environmental variables matched to sites
  by name, and a `plot()` method (ggplot2) with explained and residual sums of squares; `nested_variance()` gives REML variance components across nested levels (nlme in
  Suggests) and stops when the factors are not ordered from the broadest to the finest level,
  with a `plot()` method that draws stacked bars of the shares in Okabe-Ito colours (ggplot2
  in Suggests);
  `community_index()` compares any trait index of communities with four null models, with all
  traits of an individual kept together in multivariate indices. The standardised effect size
  is now `NA` when the null values do not vary, as for an index that a randomisation cannot
  change.

- `population_shifts()` estimates local adaptation (`h`) and character displacement (`c`)
  together from the population means of individuals, by a regression within species on the
  environment of the site and on the push predicted from co-occurring competitors, with the
  individuals of each population split in halves against the bias of sampling error and the
  kernel width profiled. `tstats(null_ip = "species")` adds a null model for T_IP.IC that
  exchanges population means among the populations of each species, so that its effect size
  is near zero when the site does not matter, positive under convergence and negative under
  divergence; the null model of Taudiere & Violle (2016) is far below zero for any set of
  distinct species. In `data-raw/pilot_tstats_identifiability.R` (ABC, `h ~ U(0, 1)`,
  `c ~ U(0, 10)`, ten individuals per population, metacommunities of the pseudo-observations
  left out), the T-statistics alone recover `h` (R2 0.87 to 0.89 over tolerances) but
  hardly `c` (0.56 to 0.58), while the estimates of `population_shifts()` recover both (0.99,
  and 0.95 to 0.96; recoverable at tolerance 0.05). `group_sum()`, used inside `tstats()`, now gives zero to empty groups, and `tstats()`
  no longer fails on a trait without values.

- `sample_individuals(displacement = , displacement_width = )` adds character displacement:
  the populations of each site sit at the equilibrium between selection toward their mean
  without competitors and a push away from close competitors through a Gaussian competition
  kernel, weighted by their abundance. `tstats()` now has the three null models of Taudiere &
  Violle (2016) in base R (`nperm`, `alternative`), with standardised effect sizes and
  p-values, so cati is needed only for comparison. The null model for T_IP.IC in
  `cati::Tstats()` 0.99.6 puts together permutations of different sites, so its effect sizes
  for T_IP.IC are wrong (adrientaudiere/cati#12); its regular-expression match of site names
  was also reported upstream (adrientaudiere/cati#11). In `data-raw/demo_tstats.R`, character
  displacement lowers T_IP.IC and raises T_PC.PR and offsets local adaptation: with
  `local_adaptation = 0.9`, `displacement = 10` takes T_IP.IC from 0.61 to 0.011 and T_PC.PR
  from 0.010 to 0.59.

- New vignette "Rates of trait evolution that vary among branches" (`branch-rates`): rates as a
  model (`trait_rates()`), from a history along the tree (`history_rates()`) and from
  ecological dynamics (`turnover_rates()`, `food_web_rates()`), and how well two estimators
  recover them, computed in advance by `data-raw/vignette_rates.R`. RRphylo follows the change
  that happened on each branch (mean Spearman rho 0.35 to 0.37) but hardly the true rate (0.06
  to 0.08, against a ceiling of 0.24 to 0.35 set by one realisation per branch); Brownian
  motion with regimes in mvMORPH recovers the ratio between the rates of two habitats (medians
  1.02, 2.12 and 5.46 for true ratios 1, 2 and 5); and the RRphylo rates of tip branches
  follow the turnover of partners weakly (median rho 0.10).

- `fit_damocles(shared = TRUE)` and `damocles_test(shared = TRUE)` estimate one pair of
  DAMOCLES rates for several communities on the same tree from the sum of their
  log-likelihoods. A single community rarely holds enough information on colonisation and
  extinction; pooling communities narrows the estimates (`data-raw/validation_damocles.R`),
  provided they can be treated as independent realisations of the process.
- DAMOCLES rates can change across the tree. `simulate_damocles()` and `damocles_loglik()`
  take one rate per edge; `fit_damocles()` and `damocles_test()` estimate rates that differ
  among `regimes` of the edges or change with an `edge_covariate`, for colonisation,
  extinction or both (`vary`). Regimes can be clades defined a priori (`clade_regimes()`) or
  rate shifts from methods of diversification (BAMM, ClaDS and similar), given as plain
  vectors so that no such package is needed. `damocles_lrt()` compares nested fits (constant
  against varying rates, shared against separate, allopatric against estimated `sympatry`)
  and `confint()` gives profile-likelihood intervals. Type I error and power of the
  comparison with constant rates are in `data-raw/validation_damocles.R`.
- `simulate_damocles()`, `damocles_loglik()`, `fit_damocles()` and `damocles_test()` take
  `sympatry`, the probability that a locally present ancestor leaves both daughters present
  at speciation (in situ speciation, or recolonisation before the daughters diverge),
  relaxing the strictly allopatric rule of Pigot & Etienne (2015), which stays the default
  (`sympatry = 0`). The likelihood remains exact; `fit_damocles(sympatry = NULL)` estimates it,
  which needs several communities with shared rates (`data-raw/validation_damocles.R`).

- `food_web_rates()` turns the rates of trait change of the resource, consumer and predator
  from `food_web_dynamics()` into rates per branch for `trait_rates(edge_rates = )`: each tip
  belongs to a trophic level, and over the last `window` units of time of its branch the rate
  is `sigma2 * (1 + k * speed)`, where speed is the mean rate of change of the mean trait of
  that level while present. As in `turnover_rates()`, only tip branches are informed; the
  speeds are directional changes under selection, so they give relative rates among levels,
  not Brownian variances.

- Community assembly with several traits. `assemble_community(distance = )` measures trait
  space by Euclidean distance (`"euclidean"`, the default and the earlier behaviour),
  Euclidean distance on traits standardised by their standard deviation in the pool
  (`"standardized"`), or Mahalanobis distance with the covariance of the pool or a given
  `trait_cov` (`"mahalanobis"`), in both limiting similarity and filtering; filtering by
  Mahalanobis distance is filtering by a multivariate Gaussian niche. The help explains when
  each fits: Euclidean distances count correlated traits more than once.
  `growth_rates()` and `simulate_metacommunity()` take a niche with several environmental
  axes: a list of traits, one per axis, and, in `simulate_metacommunity(env = )`, a matrix
  with one column per axis; the growth rate is the product of the responses to each axis,
  as in `simulate_static()`.

- `traits_ou()` simulates several correlated traits: `sigma2` can be an evolutionary rate
  matrix, or a named list with one matrix (or one rate) per selective regime, as in the
  multi-rate models of mvMORPH (Clavel et al. 2015). `regimes` can also be a mapped history
  (class `simmap`, such as from `trait_mk()`), so that regimes, optima and rates change within
  branches. With one number as `sigma2` the values are the same as before for the same seed.
  `rate_matrix()` builds rate matrices with the structures that `mvMORPH::mvBM()` fits: any
  matrix (`"default"`), independent traits (`"diagonal"`, `"equaldiagonal"`), equal rates with
  correlations (`"equal"`), and, across regimes, proportional matrices (`"proportional"`),
  shared eigenvectors (`"shared"`), equal variances (`"variance"`) or equal correlations
  (`"correlation"`). `data-raw/validation_mvtraits.R` checks them against mvMORPH 1.2.3: over
  20000 simulations with two regimes changing within branches, the covariances of the tips
  stayed within 2.9 Monte Carlo standard errors of `sum_k kronecker(R_k, C_k)` (mvSIM: 2.8),
  and under OU within 3.1 of the covariance of Hansen (1997); `mvBM()` with the same constraint
  recovered the matrices (100 tips, median of 30 fits within 9% of each true element).

- `sample_individuals()` draws individuals from the populations of a simulated metacommunity,
  with trait values made of the optimum of the species, a fraction of the gap to the
  environment of the site (`local_adaptation`), a deviation of the population that does not
  track the environment (`population_sd`) and variation among individuals
  (`intraspecific_sd`). `tstats()` computes the T-statistics of Violle et al. (2012) per site,
  in base R for the observed values or through `cati::Tstats()` (Taudiere & Violle 2016, now
  in Suggests) for its null models; both give the same values. Sites are passed to cati under
  labels of equal width, because `cati::Tstats()` matches site names as regular expressions
  and would pool `site1` with `site10`. `data-raw/demo_tstats.R` shows that local adaptation
  raises T_IP.IC (from 0.02 to 1.0 as it goes from none to complete) and lowers T_IC.IR and
  T_PC.PR (from 0.54 and 0.89 to below 0.01), means of 20 metacommunities.

- `simulate_damocles()`, `damocles_loglik()`, `fit_damocles()` and `damocles_test()` implement
  the dynamic null model of community assembly of Pigot & Etienne (2015, DAMOCLES): each
  lineage colonises the local community at rate gamma and goes locally extinct at rate mu along
  the tree, and at each allopatric speciation a present ancestor leaves one present daughter.
  The likelihood follows the pruning algorithm of the article, the rates are fitted by maximum
  likelihood with `optim()`, and `damocles_test()` compares MPD and MNTD of each community with
  communities simulated under the fitted rates and with random draws of the same richness. All
  is base R, from the published equations. `data-raw/validation_damocles.R` checks it against
  the DAMOCLES package (2.3): log-likelihoods matched `DAMOCLES_loglik()` to 3e-14 in 180
  comparisons (20 trees, three pairs of rates, three priors at the root, some trees with a root
  edge); maximum likelihoods matched `DAMOCLES_ML()` to 5e-4 in ten communities, ours equal or
  higher in all; and over 4000 simulated communities per setting the frequency of presence of
  each of 40 tips, the mean richness, MPD and MNTD agreed with `DAMOCLES_sim()` (largest
  difference in a tip frequency 0.032, |z| < 3). The rates themselves are poorly identified
  from one community of 60 species (true gamma = mu = 1; estimated mu from zero to thousands).

- `partner_turnover()` counts, for each species of a simulated metacommunity, the entries and
  exits of other species (all of them, or only competitors) in the sites it occupies, per site
  and time step between recorded steps. `turnover_rates()` turns this turnover into rates per
  branch for `trait_rates(edge_rates = )`, `sigma2 * (1 + k * turnover)` over the last
  `window` units of time of each tip branch, linking the assembly of the metacommunity to
  rates of trait evolution (evolution speeds up as the community around a lineage changes).
  The simulation runs among the tips in ecological time, so only tip branches are informed and
  the match between ecological and branch time is set by `window`.

- `plot.food_web()` draws the abundances, mean traits and rates of trait change from
  `food_web_dynamics()` in stacked panels, with introductions marked; with `time` it draws the
  state up to that time on the axes of the whole run. `animate_food_web()` turns these frames
  into a GIF (gifski in `Suggests`); the reference page of `food_web_dynamics()` shows one,
  made by `data-raw/food_web_gif.R`.

- `food_web_dynamics()` simulates the joint dynamics of abundances and mean traits of a
  resource, a consumer and a predator, the three-trophic-level model of chapter 3 of McPeek
  (2017, *Evolutionary Community Ecology*): quadratic intrinsic rates, Holling type II feeding,
  optional omnivory, and attack coefficients that are Gaussian, logistic or Michaelis-Menten
  functions of the traits. Mean traits follow the selection gradients, computed numerically
  from the fitness of a phenotype among individuals at the mean; everything is integrated by
  a fixed-step Runge-Kutta method in base R. Levels can be absent or introduced at given
  times, and the output keeps the selection gradients and the rates of trait change, for
  later use as rates per lineage. `food_web_params()` holds the parameters with the names of
  the book. `data-raw/validation_food_web.R` reproduces figure 3.2 (consumer peak of 10.2 at
  iteration 198, final traits 17.1 and 15.6, attack coefficient 0.091, rates of trait change
  280 and 47 times lower at iteration 3000 than at their maxima), figure 3.3 (the consumer at
  a maximum and the resource at a minimum of their fitness surfaces at the end) and the
  effect of the predator in figure 3.11: with the predator the resource evolves a trait
  closer to its birth optimum of 12 under the three attack forms, from 16.8 to 11.8 with
  Gaussian attack (book: 15.6 to 12.4), from 13.6 to 12.6 with logistic attack (15.4 to 12.9)
  and from 10.1 to 11.6 with Michaelis-Menten attack (10.1 to 11.5). The three-level
  Michaelis-Menten panel is reproduced closely; the consumer and predator traits of the
  logistic panel cannot be reproduced from the parameters of its caption, at which the
  consumer cannot persist. Against the MATLAB code of the resource-consumer models published
  with McPeek (2017, *The American Naturalist*), the fitness functions and gradients are the
  same and the trajectories of figure 3.2 agree within 1e-4 of the author's `ode45`
  integration; that code has no predator, and the chapter's own code was offline and not
  archived. `trait_floor` keeps mean traits at or above a value, as the author's code does
  at zero (the default with Michaelis-Menten attack).

- `history_rates()` turns the history of a discrete character along the tree, such as the
  region or habitat of each lineage from `trait_mk()` or a stochastic map, into the rate of
  each branch for `trait_rates(edge_rates = )`: a rate per state (`"state"`), a burst after
  each change that decays exponentially (`"burst"`), or a rate that depends on the number of
  lineages in the same state (`"sympatry"`, as in the diversity-dependent model with
  biogeography of Drury et al. 2016). The rate of a branch is its mean over time.
  `data-raw/validation_sympatry.R` checks `"sympatry"` against RPANDA: the covariance matrices
  of the tips matched those of the diversity-dependent model with biogeography to 3e-15 in 30
  stochastic maps, and `RPANDA::fit_t_comp()` recovered r from traits simulated with these
  rates (60 tips, three regions; mean estimate -0.30 for a true -0.3, 80% of estimates between
  -0.32 and -0.27, and 0.005 for a true 0; mean sigma2 0.98 and 0.91 for a true 1).

- `trait_rates()`: `shift_rate` draws shifts as a Poisson process along the branches, as under
  the prior of BAMM, so they can fall within branches and the rate of a branch is its mean over
  time; `relative = FALSE` gives each new regime a rate of its own, as in BAMM and bifrost,
  instead of multiplying the rate before the shift; `shift_size` can be a function that draws
  the multipliers. `data-raw/validation_bamm.R` runs BAMM 2.5.0 on 20 trees of 100 tips with
  about seven such shifts each (10 million generations, lowest effective sample size 62): the
  mean posterior rate per branch correlated with the true rate (Spearman 0.51, Pearson 0.67 on
  the log scale), and the marginal probability of a shift was 0.18 on branches with a true
  shift against 0.02 on the others.

- `trait_rates()` simulates a continuous trait whose rate varies among branches and keeps the
  rates, as the known truth for methods that estimate rates per branch. Three models: shifts
  at a few branches, inherited by their descendants or not (`"shift"`); a logarithm of the
  rate that evolves by Brownian motion (`"gbm"`); and independent lognormal rates
  (`"iid"`). It gives the trait at tips and nodes, the rate of each branch and the rescaled
  tree.
- `trait_mk()`: `alpha` multiplies the rates of each branch by a gamma factor with mean 1, the
  model of `phytools::fitgammaMk()`, or `edge_rates` gives the factors; `rate_classes` adds
  hidden rate classes, as in corHMM, and keeps the hidden history in `$hidden`. Without them
  the results are the same as before with the same seed.
- `data-raw/validation_rates.R` compares the truth with RRphylo, `phytools::multirateBM()`,
  traitMedusa (motmot), `fitgammaMk()` and corHMM. With one change per branch, rates per
  branch are weakly identifiable: the Spearman correlation between the true change and the
  true rate of a branch was 0.29 to 0.40, and RRphylo rates followed the change (0.32 to 0.38)
  but hardly the rate (0.07 to 0.13). multirateBM recovered the log rates at nodes with
  r = 0.62 and slope 0.46. traitMedusa found a true eightfold shift in 20 of 30 trees, with a
  median estimated multiplier of 8.3. fitgammaMk gave a median alpha of 0.33 for a true 0.3
  (500 tips), and corHMM preferred two hidden classes in 10 of 12 trees with classes and in
  none of 12 without.

- `simulate_metacommunity()` and `mc_dynamics()`: `immigration = "multinomial"` distributes
  the emigrants of each site among the others by the dispersal kernel and conserves them, as
  mcomsimr does, on spatially explicit landscapes. The default stays `"poisson"`, as in the
  code of Thompson et al. (2020). With extirpation 0.02 and dispersal 0.001 on a grid, mean
  local richness was 1.51 with Poisson, 1.61 with multinomial and 1.61 in mcomsimr.
- purrr moved from `Imports` to `Suggests`: it only gives the progress bar of long loops in
  interactive sessions. The package now imports only ape and base packages.

# comdynamite 0.6.3

- `species_traits()`: `optima`, `sigma` and `r_max` take values from elsewhere, such as mvMORPH,
  by species name or in the order of the species. The traits not given are drawn as before, and
  with the same seed they are the same as without the external ones. `rescale_optima = FALSE`
  keeps external optima on their own scale.
- Herodotools: `as_herodotools()` gives its inputs from a simulation or an `mc_data` object,
  with the true ancestral areas from a history of `trait_mk()`; `region_pool()` limits each
  species to the sites of its region; `arrival_ages()` gives the true ages at which lineages
  arrived in each region. The article "Known biogeographic histories for Herodotools" shows
  the steps and the results of `data-raw/validation_herodotools.R`.
- Trait histories along the tree, to validate reconstructions: `trait_bm(nodes = TRUE)` and
  `species_traits(nodes = TRUE)` give the values at the internal nodes, drawn after everything
  else so that the other values stay the same with the same seed; `traits_ou()` keeps them as
  the attribute `nodes`. `trait_mk()` simulates a discrete trait under the Mk model and keeps
  every change, in the `simmap` format of phytools. `tip_history_metrics()` gives the
  tip-based metrics of Luza et al. (2021), as computed by Herodotools, from a known history,
  plus the number of changes and the time since the last one in continuous time.
- `simulate_metacommunity()` and `mc_dynamics()`: `pool` limits each species to the sites
  where it can occur, such as its range or a biogeographic reconstruction, and `initial` gives
  the abundances or presences at the start. Without them, results with the same seed are the
  same as before.
- `summary_stats()`: new statistic `trait_k`, Blomberg's K of the trait over the species shared
  by `trait` and `tree`, as `picante::Kcalc()` gives it. It is `NA` without a tree or a trait.

# comdynamite 0.6.2

- `mc_data()` and `as_mc_data()`: community, tree, traits, environment, coordinates, areas and
  time series in one object, with names checked once. The analyses and `summary_stats()` take an
  `mc_data` object as first argument and fill what the call leaves out. Simulations keep the tree
  their traits came from, as `tree`.
- `run_scenarios()`: grids of scenarios with replicates paired by what `setup` draws, a seed per
  cell and optional parallel runs.
- `equilibrium_check()`: warns when a simulation is still changing at the end of its record.
- `as.data.frame()` for `mc_sim`: long table from `history`.
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
