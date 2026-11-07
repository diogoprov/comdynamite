# Analysing clades and traits in a metacommunity

This article goes through the two analyses of metacommunity
phylogenetics in the package: the node regression of Leibold et
al. (2010) with their daphniids, and the phylogenetic gradient analysis
of Peres-Neto et al. (2012) with simulated data, which are the methods
checked in the article on validation.

``` r

library(comdynamite)
```

## Which clades respond to the environment?

`daphniids` rebuilds the data of Leibold et al. (2010) from the public
survey: presences of 13 species in lakes of the northeastern United
States, nine predictors and six biogeographic contrasts.
[`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
codes each node of the tree as a contrast between its two clades. This
tree has no branch lengths; with them, `branch_lengths = TRUE` gives the
coding of the supplement of the article.

``` r

nodes <- node_contrasts(daphniids$tree)
dim(nodes)
#> [1] 13 12
```

[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
regresses each contrast on the predictors. The default of 9999
permutations suits the false discovery rate. We use 999 here to keep the
page quick to build, with nodes in parallel on two cores.

``` r

set.seed(1)
out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography, n_perm = 999,
                 cores = 2)
out
#> <mc_node_test> 12 nodes, 11 tested
#> permutations: 999 free of sites, 999 of species; smallest possible p: 0.001 
#> nodes with p (max) <= 0.05: 1 before correction, 0 after
#> nodes not tested: 1 (a clade without individuals, too few sites or a predictor without variance in the sites of the node)
```

The p-values answer different questions, so this choice comes first.
`p_max` tests the association between predictors and phylogeny, and it
needs both permutations to reject. With 13 species, permutations of
species have little power, so nothing passes it here. Permuting sites
alone tests whether a predictor matters for that contrast, whether or
not that comes from phylogeny. The figures below use sites, as this
recovers most of the results of the article, and should be read with
that limit in mind.

[`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)
places the results on the tree, as in Figure 3 of the article.

``` r

plot_node_test(out, daphniids$tree, test = "site")
```

![Phylogeny of daphniids with blue circles for predictors and red
triangles for biogeographic events at the nodes, labelled with the names
of the
variables.](analysing-clades_files/figure-html/unnamed-chunk-5-1.png)

A tree of this size reads best as a phylogram. For very large trees,
`type = "fan"` draws the same plot in less space.

[`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)
gives the direction and size of each association, as in Table S1 of the
article.

``` r

plot_node_heatmap(out, test = "site")
```

![Table of nodes by variables coloured by the semipartial correlation,
with significant values in
bold.](analysing-clades_files/figure-html/unnamed-chunk-6-1.png)

A positive value means that the second clade of the node, in the order
of the tree, is favoured where the predictor is high.

## Coding biogeographic events

[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
builds the events from the area of each site. Events should come from a
hypothesis stated in advance, as the areas on the two sides of a barrier
(Leibold et al. 2010; Leibold & Chase 2018). Here, three regions are
separated by two barriers, and the later one splits the last two
regions.

``` r

region <- rep(c("north", "centre", "south", NA), c(8, 6, 6, 10))
events <- list(barrier = list("north", c("centre", "south")),
               river = list("centre", "south"))
head(biogeography_matrix(region, events, min_sites = 3))
#>      barrier river
#> [1,]       1     0
#> [2,]       1     0
#> [3,]       1     0
#> [4,]       1     0
#> [5,]       1     0
#> [6,]       1     0
```

`contrast = "others"` compares each area with the rest, for example each
refugium with the other refugia. It needs no hypothesis about how areas
relate, so it asks a weaker question.

Sites of an area are often neighbours, and free permutations then reject
too often. We simulated 1000 metacommunities on a torus of 15 by 15
sites, with communities that followed a hidden autocorrelated
environment and two blocks of 25 sites placed at random, so that events
had no effect (`data-raw/validation_events.R`). Free permutations
rejected 45.1% of the tests of a node at the 5% level, and permutations
from
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
given in `site_permutations`, rejected 4.8%. The events in the figures
above used free permutations, so they should be read with this in mind.
Landscapes with edges need a restricted null model chosen for them, and
a matrix of permutations of any origin can be given.

## How much of the phylogenetic structure is due to environment and space?

For this question we need a tree with branch lengths and sites with
coordinates, so we simulate them with the generator of Appendix C of
Peres-Neto et al. (2012). In the first scenario, two independent
gradients act on species.

``` r

set.seed(1)
tree <- treeshapesim::sim_tree_shape(32, beta = 10)
gradients <- logistic_gradients(1)
sim <- simulate_logistic(tree, gradients[c("x", "y")])
keep <- rowSums(sim$abundance) > 0
comm <- sim$abundance[keep, ]
comm <- comm[, colSums(comm) > 0]
x <- gradients$x[keep]
y <- gradients$y[keep]
```

A high `beta` gives a balanced tree, as in that appendix.
[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
tests whether a gradient is related to the phylogeny and, if so, gives a
vector of species values for the next step. Here we ask this for each
gradient.

``` r

first <- phylo_vector(comm, tree, x, n_perm = 199)
second <- phylo_vector(comm, tree, y, n_perm = 199)
c(first$p_max, second$p_max)
#> [1] 0.005 0.010
```

[`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
then splits the variation of these vectors among sites. The vectors are
chosen to fit the gradients, so this step only makes sense after a
significant global test.

``` r

vectors <- cbind(first$pev_brown, second$pev_brown)
partition <- trait_partition(comm, vectors, x, y, n_perm = 199)
partition
#> <mc_trait_partition> mean component, 0.476 of the total variation
#> permutations: 199 free of sites, 199 of species; smallest possible p: 0.005 
#>   fraction     r2 r2_adjusted r2_component p_site p_taxa p_max
#> 1      abc 0.3754       0.331        0.789  0.005  0.005 0.005
#> 2       ab 0.2076       0.180        0.436  0.005  0.005 0.005
#> 3       bc 0.2348       0.208        0.493  0.005  0.005 0.005
#> 4        a 0.1405       0.122        0.295  0.005  0.005 0.005
#> 5        b 0.0671       0.058        0.141     NA     NA    NA
#> 6        c 0.1678       0.150        0.352  0.005  0.005 0.005
#> 7        d 0.6246       0.669        0.211     NA     NA    NA
round(partition$semipartial, 2)
#>        trait1 trait2
#> env1     0.50  -0.17
#> space1  -0.27   0.51
```

``` r

plot_trait_partition(partition, labels = c("Gradient x", "Gradient y"))
```

![Two overlapping bars along the total variation, for environment and
space, sharing the fraction b, followed by the unexplained
fraction.](analysing-clades_files/figure-html/unnamed-chunk-11-1.png)

The test of each trait against each predictor is a different question,
for which ter Braak (2019) is the reference.

## References

Leibold, M. A., & Chase, J. M. (2018). *Metacommunity Ecology*.
Monographs in Population Biology 59. Princeton University Press.
<https://doi.org/10.2307/j.ctt1wf4d24>

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
<https://doi.org/10.1111/j.1461-0248.2010.01523.x>

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
<https://doi.org/10.1890/11-0494.1>

ter Braak, C. J. F. (2019). New robust weighted averaging- and
model-based methods for assessing trait-environment relationships.
*Methods in Ecology and Evolution*, 10, 1962–1971.
<https://doi.org/10.1111/2041-210X.13278>
