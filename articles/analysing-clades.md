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
codes each node of the tree as a contrast between its two clades.
`tree_leibold` is the topology of Figure 3a, without branch lengths;
with them, `branch_lengths = TRUE` gives the coding of the supplement of
the article. (`daphniids$tree` is a dated tree from TimeTree with a
different topology.)

``` r

nodes <- node_contrasts(daphniids$tree_leibold)
dim(nodes)
#> [1] 13 12
```

[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
regresses each contrast on the predictors. The 9999 permutations of the
article suit the false discovery rate. We use the default of 999 here to
keep the page quick to build, with nodes in parallel on two cores.

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

plot_node_test(out, daphniids$tree_leibold, test = "site")
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

A torus does not fit real landscapes, which have edges.
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
splits each area into patches of neighbouring sites and places each
patch around a site drawn at random, which keeps the number and size of
patches but not their shape or place.

``` r

area_moves <- patch_permutations(coords, events, n_perm = 999)
node_test(comm, nodes, env, events, area_permutations = area_moves)
```

We evaluated it in a square of 20 by 20 sites, a strip of 50 by 8 and an
irregular band of 436 sites, with 200 replicates of each scenario
(`data-raw/validation_edges.R`). Areas were two compact blocks,
irregular areas of several patches, or one such area against all other
sites. Without effect of areas, free permutations rejected 44% to 55% of
the tests of a node at the 5% level, and
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
3.0% to 6.3%. Permutations of a torus, used in the square although it
has edges, rejected 5.8% to 7.1%. When areas moved the environment of
their sites, the rejection rate of
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
rose to 11% to 32%. It was lowest for compact blocks and highest for one
area against the others, and the effect of an area does not reach all
nodes, so these values are not the power for a node with an effect.

The hidden environment above had a kernel bandwidth of 3 sites. With
bandwidths of 6, 12 and 25 sites, up to the scale of the landscape, and
no effect of areas,
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
rejected 2.2% to 6.0% and free permutations 49% to 75%. It was
conservative for one area against the others (2.2% to 3.9%). The torus
rejected 5.6% to 9.5% in the square.

These rates average over places of areas. Patches are placed
independently, so the distance between the two sides of an event is not
kept. With areas fixed in the square, the rate was 0.4% to 2% for
adjacent areas and 6% to 10% for areas in opposite corners or halves of
the landscape. The test is thus liberal for sides far apart and
conservative for neighbouring sides, and p-values near the threshold
need caution. With sparse sites the result depends on `gap`, which sets
what counts as neighbouring sites: with 25% to 40% of the sites of the
square and areas in opposite corners, 1.5 times the distance to the
nearest site, the former default, gave 15% to 17%, and the current
default, which follows the fourth nearest site, 5.6% to 5.9%
(`data-raw/validation_gap.R`, 150 replicates and 99 permutations). In
sparse landscapes held out from that choice the current default gave
2.1% to 5.9%.

[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
keeps those distances. It moves the whole event at once, by a rotation
or a reflection and a translation, as the torus-translation test of
Harms et al. (2001) without the wrap-around, so each side keeps its
shape and the sides keep their distance.

``` r

area_moves <- rigid_permutations(coords, events, n_perm = 999)
node_test(comm, nodes, env, events, area_permutations = area_moves)
```

With 150 replicates (`data-raw/validation_rigid.R`) and no effect of
areas, it rejected 3.7% to 4.7% for areas placed at random in the square
and the strip, 2.2% to 3.5% for adjacent and several small areas fixed
in the square, and 4.2% to 4.8% for adjacent areas, an L shape and two
areas at a short distance when the hidden environment had an exponential
covariance. With sparse sites and areas with room to move, it rejected
5.5% to 6.1%, and with areas in opposite corners it had no placement or
fewer than five. Power for adjacent areas was about twice that of
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
(24% against 13%).

The price is that an event that spans the landscape has little room to
move. Two halves, two strips along opposite edges or areas in opposite
corners had one placement besides the observed one, and most random
areas in the band had none or fewer than 20. The function then returns
fewer permutations, with a warning, or stops, and
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
prints how many each event got. With fewer than 19 no node can reach the
5% level. Such an event cannot be tested by moving it whole, by this or
by any other permutation that keeps its shape.

`independent` lets parts of an event farther apart than that distance
move on their own, which gives placements to sides far apart, but
assumes that the community is not spatially autocorrelated beyond it.

``` r

area_moves <- rigid_permutations(coords, events, n_perm = 999, independent = 10)
```

With `independent = 10`, random areas in the three landscapes gave 3.9%
to 5.5%, but areas in opposite corners gave 6% to 7% with
autocorrelation over 3 or 4 sites and 12% over 12 sites, and strips
along opposite edges 8% to 9%, against 2% to 3% with
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md).
We found no null model that is valid for sides far apart at any range of
autocorrelation, so that range has to be argued from outside the data,
and results should be compared across values.

As a guide, on a grid use
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
when the sides have room to move, and compare it with
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
and values of `independent` when they do not. Off a grid,
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
is the safer choice, for the reason below. Coordinates off a regular
lattice get an enumerated grid of rotations, reflections and
translations, with the moved sites taken to the nearest free sites.
There (`data-raw/validation_offgrid.R`, 150 replicates and 99
permutations) it rejected 4.7% and 5.9% for adjacent areas in a square
with jitter and in sites placed at random, but 6.8% to 8.0% for two
blocks placed at random in four landscapes, where it is liberal, against
1.2% to 6.3% with
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md).
With a second generator in an L-shaped region of random sites
(`data-raw/validation_nulls_check.R`) it rejected 7.0% and 8.7%, against
5.5% and 6.6%. Power for adjacent areas was 23% against 12%. A large
part of an event in an irregular region may have few moves that fit; it
then stays nearly in place while small parts move, and the function
warns.

Both are methods of this package without published evaluation, tested
only in these simulations. Coordinates should be projected, and with
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
landscapes with parts far apart may give patches that cross the gaps.

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
#>   fraction     r2 r2_adjusted r2_component statistic p_site p_species p_max
#> 1      abc 0.3754       0.331        0.789    0.3754  0.005     0.005 0.005
#> 2       ab 0.2076       0.180        0.436    0.2076  0.005     0.005 0.005
#> 3       bc 0.2348       0.208        0.493    0.2348  0.005     0.005 0.005
#> 4        a 0.1405       0.122        0.295    0.5828  0.005     0.005 0.005
#> 5        b 0.0671       0.058        0.141    0.0671     NA        NA    NA
#> 6        c 0.1678       0.150        0.352    0.6251  0.005     0.005 0.005
#> 7        d 0.6246       0.669        0.211        NA     NA        NA    NA
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
fraction.](analysing-clades_files/figure-html/unnamed-chunk-14-1.png)

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
