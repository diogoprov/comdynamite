# comdynamite

Process-based metacommunity simulations in R, with phylogenetic structure.

comdynamite implements the model of Thompson et al. (2020): density-independent
growth set by the match between each species' environmental optimum and the
local environment, Beverton–Holt competition, Poisson demographic
stochasticity, and dispersal among sites. On top of the original model it
offers:

- **phylogenetic structure**: optima, niche breadths and maximum growth rates
  with phylogenetic signal (Pagel's lambda), and competition that is stronger
  between close relatives;
- **spatially implicit or explicit landscapes**: uniform dispersal among
  exchangeable sites, or a grid (optionally a torus) with an exponential
  dispersal kernel and a spatially autocorrelated environment;
- **restricted null models** for grids: site permutations of a torus that keep
  its spatial structure;
- an environment that can change over time, and abundances recorded at
  chosen steps.

## Installation

```r
remotes::install_github("diogoprov/comdynamite")
```

## Example

```r
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

See `vignette("comdynamite")` for species sorting vs. mass effects, grids and
time-varying environments.

## Reference

Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik, R.,
Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A process-based
metacommunity framework linking local and regional scale community ecology.
*Ecology Letters*, 23, 1314–1329. https://doi.org/10.1111/ele.13568
