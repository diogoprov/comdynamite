# Block B: more replicates and parameter sweeps of the three replications

The scripts in this folder rerun the three replications of the package (Leibold et al. 2019,
Urban & De Meester 2009, Hajian-Forooshani & Chase 2026) with 500 replicates of the published
designs and sweep the central parameters of each, to check that the results and departures
reported in the replication articles do not depend on the number of replicates or on the
point of parameter space chosen (fifth limitation of the manuscript; `dev/ROADMAP.md`,
block B).

* `run_chunks.R`: chunked, resumable runner (environment variables CORES, CHUNK, OUT,
  MAX_ROWS).
* `blockB_leibold.R`, `blockB_urban.R`, `blockB_hajian.R`: designs and runs.
* `summarise_blockB.R`: the numbers below, and the aggregated tables in
  `blockB_summary.rds`. The raw chunks (about 620,000 runs, 27 MB) are not in git.

The runs took 7 h 20 min on 16 cores of a Mac (2026-11-19), from a frozen copy of the package
at commit 361bbf5. The seeds of the published designs are those of the replication scripts, so
their replicates are a subset of these, and the first 10 (Leibold) or 100 (Urban) replicates
reproduce the numbers of the articles exactly. Intervals are 95% normal intervals over runs.

## Leibold et al. (2019)

**Replicates (500 per setting).** The qualitative results held. Without evolution, 94% of
the patches kept their preadapted species at generation 5000 (`Nm = 10`). With `G >= 50`,
fitness differences within patch types were 0.01 to 0.07 and the adjusted R2 of patch type for
composition was below 0.02, while it stayed at 0.70 to 0.95 for the community mean optimum. At
`Nm = 1`, 42% to 45% of the patches were switched with `G >= 20`; at `Nm = 10`, 33% to 42% of
the patches were shared with `G >= 50`, as in the article (29% to 38%). The loss of a species
was slightly faster than estimated from 10 replicates: with `Nm = 10`, both species persisted
to generation 10,000 in 7% and 12% of the runs with `G` of 200 and 400 (2 of 10 in the
article), and in 6% of the control with a single patch type and `G = 200`.

**Sweep of patch extinction** (0 to 0.05 per generation, `Nm` of 1 to 100, 100 runs each).
Regional neutrality appeared with `G >= 50` at every extinction rate, including none, and at
every dispersal level at which both species persisted (adjusted R2 of composition at most
0.07). The fast loss of a species did not depend on the extinction rate chosen: at
`p_ext = 0`, both species persisted to generation 10,000 in 11% to 33% of the runs with
`Nm = 10` and in none with `Nm = 30`, and the loss became faster as extinction increased (no
run kept both species to generation 10,000 at `p_ext = 0.05` with `G >= 50`).
Local dominance at `Nm = 10` was also stable, with 29% to 40% of shared patches across
extinction rates. Interestingly, `G = 20` with `Nm = 1` gave regional neutrality from
`p_ext = 0.001` upwards, but not without extinction (R2 0.92), so at low dispersal patch
extinction is what lets slow adaptation erase species sorting.

**Sweep of the mutational variance** (`Vm` 0.01 to 1 from no variance, 100 runs each), the
question left open in the article. With `Vm <= 0.1`, the variance built too slowly and the
runs stayed close to species sorting (R2 of composition 0.52 to 0.96). With `Vm` of 0.5 and 1,
regional neutrality emerged (fitness differences 0.01 to 0.03, R2 of composition below 0.01),
and no species was lost in any run. Local dominance did not: 46% to 56% of the patches were
shared with `Nm = 1`, and 80% to 100% with `Nm >= 3`. A higher rate of mutation therefore
gives the same outcome as the size-dependent variance of the article, regional neutrality
without local dominance.

## Urban & De Meester (2009)

**Replicates (500 per setting).** The numbers of the article fell within the intervals of
the 500 replicates except three, all on the advantage of the first arrival. Without evolution,
it won 76% of the runs with a lead of 10 generations (interval 74% to 79%, against 84% in the
first 100) and held 64% of the new site (68% in the article). With both species evolving by
mutation, it won 61% of the runs with a lead of one generation (66% in the article) and 74%
with a lead of 10, against 64% and 76% without evolution. Evolution by mutation therefore did
not make the first arrival more likely to win, as stated in the article, and with 500
replicates the two settings are closer than the article suggests. One species held more than
90% of the new site in 34% of the runs with mutation in both species and in 9% without
evolution. With standing variation in both species, the first arrival won 29% to 39% of the
runs (the reversal), and with a constant `G` 68% to 93%.

Over the dispersal grid, the poorer disperser won in 65 of the 110 cells off the diagonal
(interval of the dominance excluding zero) and the better disperser in 3. Two of these are
adjacent (species 1 at 0.126 and 0.199 against 0.005; dominance 0.09 and 0.14), and their
mirror cells have the same sign with intervals that include zero. With 110 cells about three
would exceed the interval by chance, so this is at most weak evidence for a small
monopolization region at the edge of the grid. Gene flow back into the sources gave the same
result as with 40 replicates: the winner of the new site held 11% of its own source at
dispersal 0.01 and 69% and 75% at 0.05 and 0.10.

**Sweep** (`Vm` 0.01 to 0.25, environment of the new site 50 to 70 with sources at 40, four
genetic regimes, delays -20 to 20 by 5, 200 runs each). The advantage of a species that
evolves by mutation against one that does not held when the new site differed from the
sources by 15 units or more and `Vm <= 0.1`: arriving 10 generations late, it held 63% to 94%
of the site. It disappeared with a contrast of 10 units (28% to 49%), where adaptation adds
little and the first arrival keeps the site, and weakened with `Vm = 0.25` (8% to 82%). The reversal with
standing variation depended on the environmental contrast. At contrasts of 10 to 20 units
(the published design is 20), the first arrival won 10% to 49% of the runs; at 25 and 30 units
it won 84% to nearly 100% with `Vm <= 0.1` (but 18% at 25 units with `Vm = 0.25`), so the
priority effect of standing variation in their model reappears without erosion of the
variance when the new environment is far enough from the sources. Monopolization with mutation in both species (one species above 90%) increased
with the contrast, from 6% to 18% at 10 units to 27% to 54% at 25 and 30 (`Vm <= 0.1`).

## Hajian-Forooshani & Chase (2026)

**Replicates (500 per setting, new seeds).** The pattern of their Fig. 2 held (niche 0.25,
step 700). With dispersal limitation (`1e-5`), beta diversity increased at every extent, up to
0.60 at 85%. With a dispersal surplus (`1e-3`), it decreased at every extent, from -0.08 to
-0.93. At intermediate dispersal (`1e-4`), it decreased with 0% to 10% of the landscape
disturbed, did not change at 15%, increased from 20% to 95% (up to 0.20 at 80%) and decreased
at 100%; the article reported no change at 20%. The undisturbed control confirmed the decline
of beta diversity without disturbance at intermediate and high dispersal (-0.05 to -0.06 and
-0.07 to -0.10 across niche breadths). Relative to the control, disturbance with a dispersal
surplus lowered beta diversity by less than 0.03 up to half of the landscape disturbed, so
most of the homogenisation at small extents is background change, and at intermediate
dispersal disturbance raised beta diversity at every extent from 5% to 95%.

**Sweep of the pulse** (duration 10 to 400 steps, disturbed environment 0.1 to 0.5, dispersal
`1e-5` to `1e-1`, extent 1% to 100%, 100 runs each, beta diversity 100 steps after the end of
the pulse). Dispersal set the direction of the change whatever the pulse: beta diversity never
decreased with dispersal limitation (139 of 150 settings positive) and always decreased at
`1e-3` and `1e-2` (150 of 150). At intermediate dispersal the direction depended on the pulse:
short and severe pulses increased beta diversity, long and mild ones decreased it (at 40% of
the landscape, from 0.09 for 50 steps at 0.1 to -0.10 for 400 steps at 0.4). At `1e-1`, the
metacommunity was well mixed and the change was zero in 97 of 150 settings and small where
negative.

## What this changes in the articles and the manuscript

Nothing reverses a conclusion. The changes to consider are: (i) Leibold: the persistence
numbers with 500 replicates, the robustness to patch extinction, and the answer to the open
question on `Vm`; (ii) Urban: the advantage of the first arrival with and without evolution
(76% and 74% with a lead of 10 generations, not 84% and 72%), and the dependence of the standing-variation reversal and of the mutation-evolver
advantage on the environmental contrast, which qualifies the sentence of the manuscript on
standing variation; (iii) Hajian: no change in beta diversity at 15% rather than 15% to 20%,
the comparison with the control, and the robustness of the direction set by dispersal to the
duration and severity of the pulse.

## References

Hajian-Forooshani, Z., & Chase, J. M. (2026), Leibold, M. A., et al. (2019) and Urban, M. C.,
& De Meester, L. (2009): see the replication articles in `vignettes/articles/`.
