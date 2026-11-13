# Replicating Urban & De Meester (2009)

Urban & De Meester (2009) proposed that local adaptation can strengthen
priority effects. In their individual-based model, two initially
equivalent species start in two source patches and colonise a third,
vacant patch whose environment differs from that of the sources. The
first species to arrive adapts to the new patch and then keeps the later
species out, a process they called community monopolization. When
dispersal is high, gene flow from the sources prevents adaptation, and
the outcome turns on dispersal instead. We asked whether the
quantitative-genetic engine of comdynamite reproduces these results, and
where it departs from them. Like the [replication of Leibold et
al. (2019)](https://provetelab.org/comdynamite/articles/replication-leibold.md),
this article is part of the nesting map of the package.

## Translation of the model

The script `data-raw/replication_urban2009.R` follows the parameters of
Table S1 of their electronic supplementary material, with the
environment and the traits multiplied by 100:

- three sites with environment 40, 40 and 60, connected to each other,
  with emigrants split equally between the other two sites; species 1
  fills site 1 and species 2 fills site 2 at the start, both with
  optimum 40;
- a Gaussian niche with \sigma = 35.4, since their survival is
  \exp(-((z - e)/\omega)^2) with \omega = 0.5, and a non-genetic
  variance of 12.5 for the phenotypic and microsite deviations (standard
  deviation 2.5 each);
- Beverton-Holt growth with r\_{max} = 5 and equal competition within
  and between species, so that a site holds about 1000 individuals, as
  their 1000 microsites;
- one step per generation, taken as 10 of their time steps (the lifespan
  1/b without selection, with baseline mortality b = 0.1), so 500
  generations instead of 5000 steps, delays between the two species of
  -20 to 20 generations instead of -200 to 200 steps, and dispersal 0.05
  per generation;
- the new site is opened to each species at its introduction time
  through the argument `pool_from` of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
  added for this replication, and is then colonised by dispersal from
  the sources; the first species can enter after 20 generations; until
  then the site is not a destination of that species, so its emigrants
  go only to the other source rather than dying in a closed site;
- a genetic variance per population with mutation and drift
  (`mutation_var` in
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
  Lynch & Hill 1986): mutation adds V_m = L \mu a^2 = 0.05 per
  generation for L = 20 loci, \mu = 10^{-4} and an allelic effect a = 5,
  and drift removes 1/N of the variance per generation, as in a haploid
  population of N, so that the equilibrium in a full site is 50;
- evolution by mutation starts without variance, and evolution from
  standing variation starts with L p (1 - p) a^2 = 120 for an allele
  frequency p = 0.4; a species that does not evolve has no genetic
  variance;
- 100 replicates of each setting, as in their figures.

Their model has 20 haploid loci with alleles of 0 or 1, lottery
recruitment into microsites vacated by mortality, overlapping
generations with selection on survival, and microsite variation within
patches. None of these is translated. One departure matters for what
follows: selection does not erode the genetic variance of
`mutation_var`, so standing variation persists and lowers mean fitness
(a variance load), whereas selection on 20 loci in their model erodes
it. We therefore added two scenarios with a constant variance G = 120,
as in the rest of the package.

A first version of this replication let the emigrants that headed for
the new site before it opened die there. That loss of up to half of the
emigrants for 20 generations penalised the better disperser and, since
competition between the two species is neutral in the sources, its
effect persisted; it produced a dispersal result that we first
attributed to the Beverton-Holt engine. `pool_from` now keeps those
emigrants home, and all results below come from this corrected version.

## Colonisation delay

Without evolution, the first species to enter the new site usually
dominated it, but not always (Fig. 1), as in their Fig. 2a. The first
species to arrive ended with more individuals in 65% of the runs with a
lead of one generation and in 84% with a lead of 10 generations. When
species 1 alone evolved by mutation, it dominated the new site at every
delay, with 88% of the final abundance when it arrived 10 generations
after species 2 and 93% when it arrived 10 generations before. This
matches their result that the evolving species wins even when it arrives
last.

When both species evolved by mutation, the outcome became more extreme:
one species held more than 90% of the new site in 34% of the runs,
against 10% without evolution, which gives the U-shaped distribution of
final abundances of their model. However, the first arrival was not more
likely to win than without evolution: it won 66% of the runs with a lead
of one generation and 72% with a lead of 10 generations (65% and 84%
without evolution), and held 67% of the site with a lead of 10
generations (68% without evolution). In their model, a lead of 5 to 10
of their time steps (half a generation to one generation here) was
enough for the first arrival to dominate. Evolution by mutation
therefore amplified the dominance of whichever species took the lead,
but did not tie that lead to the order of arrival as strongly as in
their model.

![Mean final abundance of each species in the new site against the delay
between the introductions of the two species, over 100 runs. Positive
delays: species 1 enters first. Notice that with mutation in both
species the curves cross more steeply than without evolution, and that
with a constant genetic variance they cross most
steeply.](urban_delay.png)

Mean final abundance of each species in the new site against the delay
between the introductions of the two species, over 100 runs. Positive
delays: species 1 enters first. Notice that with mutation in both
species the curves cross more steeply than without evolution, and that
with a constant genetic variance they cross most steeply.

Standing variation did not reproduce their results. When species 1 alone
evolved from standing variation, it held about 70% of the new site at
every delay, whereas in their model it dominated less when it arrived
after its competitor. When both evolved from standing variation, the
priority effect was reversed: the second species to arrive held 60% of
the site with a lag of 10 generations, while in their model
monopolization was stronger than with mutation alone. This seems to come
from the variance load of `mutation_var`, since the reversal disappears
when the variance is held constant, but we have not separated how the
load acts on the first and on the second species to arrive. With a
constant G = 120 for both species, the first arrival held 81% of the
site with a lead of 10 generations and won 90% of the runs, and one
species held more than 90% of the site in 42% of the runs. This is
closer to their result for standing variation. With a constant G for
species 1 alone, species 1 held about 62% of the site at every delay.

## Dispersal

With both species evolving by mutation and entering together, dispersal
had a weaker effect than in their Fig. 4 (Fig. 2). Off the diagonal, the
poorer disperser tended to win: the relative dominance (N_1 - N_2)/1000
reached 0.2 to 0.33 when species 1 dispersed at 0.008 to 0.02 and
species 2 at 0.02 to 0.08, and the symmetric values were negative. When
both species dispersed at 0.08 or more, the dominance stayed between
-0.12 and 0.12. Their Fig. 4 has two regions: a monopolization region,
where the better disperser arrives first and wins, and a region of
gene-flow-induced maladaptation, where high dispersal prevents
adaptation and the poorer disperser wins. We found only the second, and
with a weaker effect.

![Relative proportional dominance of species 1 in the new site after 500
generations, mean of 100 runs, for each combination of dispersal of the
two species (both evolving by mutation). Notice that the poorer
disperser tends to win at intermediate dispersal, and that no region
favours the better disperser.](urban_dispersal.png)

Relative proportional dominance of species 1 in the new site after 500
generations, mean of 100 runs, for each combination of dispersal of the
two species (both evolving by mutation). Notice that the poorer
disperser tends to win at intermediate dispersal, and that no region
favours the better disperser.

Evolution is what turns the advantage around. With dispersal 0.0315 for
species 1 and 0.005 for species 2 (40 runs), species 1 held 694
individuals of the new site and species 2 held 115 without evolution, so
the better disperser won, as in their monopolization region. With both
species evolving by mutation, species 1 held 236 and species 2 held 731.
The better disperser seems to colonise first and then to be held back by
the gene flow it keeps receiving from its source, while the few founders
of the poorer disperser adapt with little gene flow. In their model, the
first colonist also adapts before the gene flow catches up, which our
engine does not reproduce, perhaps because the genetic variance of a new
population starts low with `mutation_var` and builds up only with
population size.

## Gene flow back into the sources

Their Fig. 3 shows that at high dispersal (m = 0.10) species 2 replaced
species 1 in the source of species 1. We ran 40 replicates of 2000
generations with both species evolving by mutation and dispersal 0.01,
0.05 or 0.10, and asked whether the species that won the new site kept
its own source (Fig. 3). With dispersal 0.01, the winner of the new site
held 11% of its own source after 2000 generations, and the mean genetic
optimum of its source population had moved from 40 to about 56, close to
the environment of the new site. Gene flow from the large, adapted
population of the new site seems to have maladapted the shrinking source
population, which was then taken over by the other species. With
dispersal 0.05 and 0.10, the winner kept 68% and 75% of its source, and
the loser of the new site kept only 31% and 24% of its own. The
replacement of a source by gene flow therefore happened here, but at low
rather than high dispersal, and we have not identified why the winner’s
source shrinks at the start when dispersal is low.

![Share of the winner of the new site in its own source site, and mean
genetic optimum of the winner in its source, over 2000 generations;
means of 40 runs for each dispersal.](urban_source.png)

Share of the winner of the new site in its own source site, and mean
genetic optimum of the winner in its source, over 2000 generations;
means of 40 runs for each dispersal.

## Summary

The quantitative-genetic engine reproduces two results of Urban & De
Meester (2009): a species that evolves by mutation wins a new site
against one that does not, whatever the order of arrival, and evolution
in both species makes the new site more often dominated by a single
species. It departs from their model in three ways. With mutation, the
first arrival is not more likely to win than without evolution; with
standing variation, the priority effect is reproduced only with a
constant G, since the size-dependent variance of `mutation_var` carries
a variance load that reverses it; and over the dispersal grid, evolution
makes the poorer disperser win, as in their region of gene-flow-induced
maladaptation, but no region favours the better disperser. Gene flow
from the new site replaced a source population, as in their Fig. 3, but
at low dispersal. Lottery recruitment into vacated microsites,
overlapping generations and a genetic variance built from new mutations
and eroded by selection are the missing mechanisms that seem most likely
to close these gaps. One step was taken as 10 of their time steps, and
the speed of adaptation was set by V_m and the variance rather than
matched to their loci.

## References

Lynch, M., & Hill, W. G. (1986). Phenotypic evolution by neutral
mutation. *Evolution*, 40, 915–935.
<https://doi.org/10.1111/j.1558-5646.1986.tb00561.x>

Urban, M. C., & De Meester, L. (2009). Community monopolization: local
adaptation enhances priority effects in an evolving metacommunity.
*Proceedings of the Royal Society B*, 276, 4129–4138.
<https://doi.org/10.1098/rspb.2009.1382>
