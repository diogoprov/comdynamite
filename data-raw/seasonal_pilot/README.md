# Identifiability pilot for a seasonal simulator (issue #97)

Prototype only; nothing here is exported. `seasonal_model.R` simulates pond-breeding anurans
with a monthly step (breeding windows, larval period, metamorphosis, drying of temporary ponds,
larval Beverton-Holt competition, dispersal of metamorphs through an exponential kernel,
imperfect detection of tadpoles and calling males). The four community models (SS species
sorting, NT neutral, PE strict priority effects, PM mixed priority effects) differ only in the
competition matrix.

`run_pilot.R` simulates 3000 runs per model on one design of 13 ponds and 15 species (about
3.5 min on 16 cores; `pilot.rds`, 5.8 MB, is not tracked), each summarised for 1, 2 and 3
field years. `analyse_pilot.R` runs ABC-RF model choice (`abc_choose()`, 500 trees).

Result (2026-11-17), out-of-bag prior error rate:

| Field years | Four models, tadpoles | Four models, + calling males | SS vs PM, tadpoles | SS vs PM, + calling males |
|---|---|---|---|---|
| 1 | 0.561 | 0.547 | 0.460 | 0.454 |
| 2 | 0.537 | 0.529 | 0.438 | 0.436 |
| 3 | 0.518 | 0.507 | 0.432 | 0.447 |

Chance is 0.75 for four models and 0.50 for two. PM was misclassified in 73% of runs (3 years,
both data sets), mostly as SS or PE; NT was the best identified (33% error). Error rose with
the strength of larval competition (0.43, 0.50 and 0.59 in terciles of `a0`). Prior predictive
ranges of monthly colonisation, extinction and turnover overlapped those of a field data set
of this size, but persistence of the dominant species was higher.

Caveats: one design, one prototype structure (niche acting on larval survival, competition
only among larvae) and one set of 15 statistics per data set; statistics targeted at priority
effects (for example colonisation conditional on the tenure of residents) were not tried.

Outcome: the issue was closed as not planned. Monthly data from a single landscape of this
size barely separate the community models, and further field years or calling-male counts add
little, so the seasonal simulator was not moved into `R/`.
