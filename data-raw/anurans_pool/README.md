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
