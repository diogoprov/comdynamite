# Regional pools of `anurans`

`pools.csv` lists, for each city, the anurans recorded by GBIF in its ecoregion (RESOLVE 2017,
Dinerstein et al. 2017) with the number of records, the number within 200 km of its ponds, the
lentic classification of the mode of reproduction and its source, and whether the species is in
the tree of Portik et al. (2023). `data-raw/anurans.R` selects the pools from it.

The scripts `01_ecoregions.R` to `05_build.R` built it on 11 November 2026 in a separate
working folder (paths inside are those of that folder): ecoregion polygons from the UNEP-WCMC
feature server of RESOLVE 2017; species facets from the GBIF occurrence search API (no download
DOI; the bodies of the POST queries are `query_eco*.json`, the URLs of the 200-km queries are in
`gbif_query_info.txt`); reproductive modes from the supplementary table S1 of
Nunes-de-Almeida et al. (2021, CC-BY 4.0; definitions in `nunes2021_mode_definitions.csv`) and,
for species missing from it, from da Silva et al. (2012, Dryad, CC0), whose modes follow Haddad &
Prado (2005). The correspondence of the modes of Haddad & Prado (2005) to lentic breeding was
not checked against the original table.

## GBIF download with a DOI

The facets above came from the occurrence search API, which gives no DOI. On 25 November 2026,
`06_download.R` requested one GBIF download with the five filters of `query_eco*.json` and the
three ecoregion polygons joined by OR: GBIF.org (25 November 2026) GBIF Occurrence Download
https://doi.org/10.15468/dl.g4gz8u (115,216 records from 63 datasets, against 115,023 in the
queries of 11 November; key and request in `gbif_download.txt`). The file itself is not kept
here. Assigned to the ecoregion polygons, its counts by species and ecoregion matched those of
`pools.csv` in 1,296 of 1,427 pairs (Spearman correlation of the totals by species and city
0.9998); 23 pairs had 1 to 10 records fewer, 105 had 1 to 10 more and 3 had 11 to 17 more, and
one species with a single record in the Cerrado (*Pristimantis moa*, not in any pool) was no
longer there. With the threshold of three records, three species would now pass it: of these,
only *Rhinella arenarum* in Campo Grande breeds in lentic water and would enter its pool. The
pools of `anurans` were not rebuilt, so they rest on the queries of 11 November; the download is
the citable record of the same query.
