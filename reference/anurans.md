# Anurans of urban ponds in two Brazilian cities

Counts of 44 species of anurans in 40 ponds, 20 in Campo Grande (Mato
Grosso do Sul) and 20 in São José dos Campos (São Paulo), from the
fieldwork of L. B. C. Menezes for chapter 1 of an MSc dissertation (L.
B. C. Menezes, M. R. Severgnini & D. B. Provete, unpublished), with the
regional pool of each city. The script is `data-raw/anurans.R`. The
counts come from surveys of calling and visual encounter, so their
totals depend on effort and detection, and
[`vignette("abc")`](https://provetelab.org/comdynamite/articles/abc.md)
uses incidence only.

## Usage

``` r
anurans
```

## Format

List with `comm`, a matrix of ponds by species with counts; `env`, a
matrix with the hydrology score, area (m²), depth (cm), floating
vegetation (%), water temperature (°C) and the urban cover within 1 km
in 2022 (%); `city`, a factor; `coords`, a matrix with the longitude and
latitude of each pond in decimal degrees; `pool`, a 0/1 matrix of ponds
by the 155 species of the regional pools; `tree`, a `phylo` object with
those species; and `grafted`, the species added to the tree.

## Source

GBIF.org (11 November 2026), occurrence search API; see
`data-raw/anurans_pool`.

## Details

The hydrology score is the first principal component of log(x + 1) area,
depth, floating vegetation and water temperature, standardised across
the 40 ponds and rescaled to 0–100, with larger and deeper ponds at
higher scores.

The regional pool of each city has the anurans that breed in lentic
water and were recorded both in its ecoregion (Dinerstein et al. 2017:
Cerrado for Campo Grande; Alto Paraná Atlantic forests and Serra do Mar
coastal forests for São José dos Campos) and within 200 km of its ponds,
in occurrences of GBIF (gbif.org, accessed on 11 November 2026, through
the occurrence search API; queries in `data-raw/anurans_pool`), plus the
species recorded in its ponds. A species breeds in lentic water if one
of its reproductive modes in Nunes-de-Almeida et al. (2021), or in da
Silva et al. (2012) for species missing from it, has exotrophic tadpoles
in lentic water, or, for species in neither, if at least three congeners
are classified and all agree. GBIF records of *Leptodactylus latrans*
were taken as *L. luctator*, split from it by Magalhães et al. (2020),
and the introduced *Lithobates catesbeianus* was left out.

The tree is the time-calibrated tree of Portik et al. (2023) pruned to
the species of the pools, with branch lengths in millions of years and
the names of GBIF (so *Pithecopus* and *Ololygon* instead of
*Phyllomedusa* and *Scinax*; Duellman et al. 2016). The tip of *L.
latrans* stands for *L. luctator*. The six species missing from it
(`grafted`) were joined at the crown of their genus, or of the species
of *Ololygon* in the tree, and the polytomies were resolved at random
with branches of length zero.

## References

da Silva, F. R., Almeida-Neto, M., do Prado, V. H. M., Haddad, C. F. B.,
& Rossa-Feres, D. C. (2012). Humidity levels drive reproductive modes
and phylogenetic diversity of amphibians in the Brazilian Atlantic
Forest. *Journal of Biogeography*, 39, 1720–1732.
[doi:10.1111/j.1365-2699.2012.02726.x](https://doi.org/10.1111/j.1365-2699.2012.02726.x)

Dinerstein, E., Olson, D., Joshi, A., Vynne, C., Burgess, N. D.,
Wikramanayake, E., et al. (2017). An ecoregion-based approach to
protecting half the terrestrial realm. *BioScience*, 67, 534–545.
[doi:10.1093/biosci/bix014](https://doi.org/10.1093/biosci/bix014)

Duellman, W. E., Marion, A. B., & Hedges, S. B. (2016). Phylogenetics,
classification, and biogeography of the treefrogs (Amphibia: Anura:
Arboranae). *Zootaxa*, 4104.
[doi:10.11646/zootaxa.4104.1.1](https://doi.org/10.11646/zootaxa.4104.1.1)

Magalhães, F. de M., Lyra, M. L., de Carvalho, T. R., Baldo, D.,
Brusquetti, F., Burella, P., Colli, G. R., Gehara, M. C., Giaretta, A.
A., Haddad, C. F. B., Langone, J. A., López, J. A., Napoli, M. F.,
Santana, D. J., de Sá, R. O., & Garda, A. A. (2020). Taxonomic review of
South American butter frogs: phylogeny, geographic patterns, and species
delimitation in the *Leptodactylus latrans* species group (Anura:
Leptodactylidae). *Herpetological Monographs*, 34, 131–177.
[doi:10.1655/0733-1347-31.4.131](https://doi.org/10.1655/0733-1347-31.4.131)

Nunes-de-Almeida, C. H. L., Haddad, C. F. B., & Toledo, L. F. (2021). A
revised classification of the amphibian reproductive modes.
*Salamandra*, 57, 413–427.
[doi:10.5281/zenodo.5205572](https://doi.org/10.5281/zenodo.5205572)

Portik, D. M., Streicher, J. W., & Wiens, J. J. (2023). Frog phylogeny:
a time-calibrated, species-level tree based on hundreds of loci and
5,242 species. *Molecular Phylogenetics and Evolution*, 188, 107907.
[doi:10.1016/j.ympev.2023.107907](https://doi.org/10.1016/j.ympev.2023.107907)
