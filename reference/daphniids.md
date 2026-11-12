# Daphniids of the northeastern United States

Incidence of 13 daphniids in lakes of the EMAP Northeastern Lakes Survey
of 1991 to 1994 (first visit to each lake), rebuilt to follow Leibold et
al. (2010). The tables are public data of the US Environmental
Protection Agency, and the script is `data-raw/daphniids.R`.

## Usage

``` r
daphniids
```

## Format

List with `comm`, a matrix of lakes by species; `env`, a matrix with
colour, conductivity, dissolved organic carbon, total nitrogen, pH,
total phosphorus, mean depth, elevation and the ratio of nitrogen to
phosphorus; `biogeography`, a matrix of contrasts with values of 1, -1
and 0; `tree`, the dated `phylo` object; `tree_leibold`, the topology of
Figure 3a; and `coords`, a matrix with the longitude and latitude of
each lake in decimal degrees.

## Source

<https://archive.epa.gov/emap/archive-emap/web/html/nelakes.html>

## Details

The article gives 327 lakes and this rebuilt set has 314, of which 239
have some species. The survey lists *D. longiremis* and not the *D.
longispina* of Figure 3a of the article, which we took to be the same
taxon.

`tree_leibold` is the topology of Figure 3a, without branch lengths, for
analyses that reproduce the article. `tree` has branch lengths in
millions of years from TimeTree 5 (Kumar et al. 2022), which dates nine
of the 13 species, mostly from mitochondrial genomes (Xu et al. 2021).
Its topology differs from that of Figure 3a: *D. laevis* and *D. dubia*
are the sister of *D. galeata mendotae* and *D. rosea* rather than of
all other *Daphnia*. The four species without dates were placed as in
Figure 3a, at arbitrary ages: *D. catawba* and *D. parvula* split the
branch of *D. retrocurva* into thirds, *D. rosea* joins *D. galeata
mendotae* halfway along its branch, and *S. serrulatus* joins *S.
vetulus* at the crown age of *Simocephalus* in TimeTree. The absolute
ages are uncertain: for the ancestor of the *D. longispina* and *D.
pulex* groups, Cornetti et al. (2019) found 105 Ma with mitochondrial
genes and 31.5 Ma with nuclear genes, so relative ages are more reliable
than absolute ones. The article does not give the area of each lake, so
rivers in `biogeography` are approximated by longitude and these
contrasts differ from the original ones: S1, east and west of the
Hudson; S2, Adirondacks and other lakes to the west of the Hudson; S3,
east and west of the Connecticut River, east of the Hudson; S4, Maine
and other lakes east of the Hudson; S5, New Hampshire and lakes to its
west but east of the Hudson; S6, New Hampshire and lakes in Rhode Island
and Massachusetts.

## References

Cornetti, L., Fields, P. D., Van Damme, K., & Ebert, D. (2019). A
fossil-calibrated phylogenomic analysis of *Daphnia* and the Daphniidae.
*Molecular Phylogenetics and Evolution*, 137, 250–262.
[doi:10.1016/j.ympev.2019.05.018](https://doi.org/10.1016/j.ympev.2019.05.018)

Kumar, S., Suleski, M., Craig, J. M., Kasprowicz, A. E., Sanderford, M.,
Li, M., Stecher, G., & Hedges, S. B. (2022). TimeTree 5: an expanded
resource for species divergence times. *Molecular Biology and
Evolution*, 39, msac174.
[doi:10.1093/molbev/msac174](https://doi.org/10.1093/molbev/msac174)

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

Taylor, D. J., Finston, T. L., & Hebert, P. D. N. (1998). Biogeography
of a widespread freshwater crustacean: pseudocongruence and cryptic
endemism in the North American *Daphnia laevis* complex. *Evolution*,
52, 1648–1670.
[doi:10.1111/j.1558-5646.1998.tb02245.x](https://doi.org/10.1111/j.1558-5646.1998.tb02245.x)

Xu, S.-L., Han, B.-P., Martínez, A., Schwentner, M., Fontaneto, D.,
Dumont, H. J., & Kotov, A. A. (2021). Mitogenomics of Cladocera
(Branchiopoda): marked gene order rearrangements and independent
predation roots. *Molecular Phylogenetics and Evolution*, 164, 107275.
[doi:10.1016/j.ympev.2021.107275](https://doi.org/10.1016/j.ympev.2021.107275)
