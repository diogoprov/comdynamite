#' Daphniids of the northeastern United States
#'
#' Incidence of 13 daphniids in lakes of the EMAP Northeastern Lakes Survey of 1991 to
#' 1994 (first visit to each lake), rebuilt to follow Leibold et al. (2010). The tables
#' are public data of the US Environmental Protection Agency, and the script is
#' `data-raw/daphniids.R`.
#'
#' The article gives 327 lakes and this rebuilt set has 314, of which 239 have some
#' species. The survey lists *D. longiremis* and not the *D. longispina* of Figure 3a of the
#' article, which we took to be the same taxon.
#'
#' `tree_leibold` is the topology of Figure 3a, without branch lengths, for analyses that
#' reproduce the article. `tree` has branch lengths in millions of years from TimeTree 5 (Kumar et al. 2022), which
#' dates nine of the 13 species, mostly from mitochondrial genomes (Xu et al. 2021). Its
#' topology differs from that of Figure 3a: *D. laevis* and *D. dubia* are the sister of
#' *D. galeata mendotae* and *D. rosea* rather than of all other *Daphnia*. The four species
#' without dates were placed as in Figure 3a, at arbitrary ages: *D. catawba* and
#' *D. parvula* split the branch of *D. retrocurva* into thirds, *D. rosea* joins
#' *D. galeata mendotae* halfway along its branch, and *S. serrulatus* joins *S. vetulus* at
#' the crown age of *Simocephalus* in TimeTree. The absolute ages are uncertain: for the
#' ancestor of the *D. longispina* and *D. pulex* groups, Cornetti et al. (2019) found
#' 105 Ma with mitochondrial genes and 31.5 Ma with nuclear genes, so relative ages are more
#' reliable than absolute ones. The article does not give the area of each lake, so rivers
#' in `biogeography` are approximated by longitude and these contrasts differ from the
#' original ones: S1, east and west of the Hudson; S2, Adirondacks and other lakes to the
#' west of the Hudson; S3, east and west of the Connecticut River, east of the Hudson; S4,
#' Maine and other lakes east of the Hudson; S5, New Hampshire and lakes to its west but
#' east of the Hudson; S6, New Hampshire and lakes in Rhode Island and Massachusetts.
#'
#' @format List with `comm`, a matrix of lakes by species; `env`, a matrix with colour,
#'   conductivity, dissolved organic carbon, total nitrogen, pH, total phosphorus, mean
#'   depth, elevation and the ratio of nitrogen to phosphorus; `biogeography`, a matrix of
#'   contrasts with values of 1, -1 and 0; `tree`, the dated `phylo` object;
#'   `tree_leibold`, the topology of Figure 3a; and `coords`, a matrix with
#'   the longitude and latitude of each lake in decimal degrees.
#' @source \url{https://archive.epa.gov/emap/archive-emap/web/html/nelakes.html}
#' @references
#' Cornetti, L., Fields, P. D., Van Damme, K., & Ebert, D. (2019). A fossil-calibrated
#' phylogenomic analysis of *Daphnia* and the Daphniidae. *Molecular Phylogenetics and
#' Evolution*, 137, 250--262. \doi{10.1016/j.ympev.2019.05.018}
#'
#' Kumar, S., Suleski, M., Craig, J. M., Kasprowicz, A. E., Sanderford, M., Li, M., Stecher, G.,
#' & Hedges, S. B. (2022). TimeTree 5: an expanded resource for species divergence times.
#' *Molecular Biology and Evolution*, 39, msac174. \doi{10.1093/molbev/msac174}
#'
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity phylogenetics:
#' separating the roles of environmental filters and historical biogeography. *Ecology
#' Letters*, 13, 1290--1299. \doi{10.1111/j.1461-0248.2010.01523.x}
#'
#' Taylor, D. J., Finston, T. L., & Hebert, P. D. N. (1998). Biogeography of a widespread
#' freshwater crustacean: pseudocongruence and cryptic endemism in the North American
#' *Daphnia laevis* complex. *Evolution*, 52, 1648--1670.
#' \doi{10.1111/j.1558-5646.1998.tb02245.x}
#'
#' Xu, S.-L., Han, B.-P., Martínez, A., Schwentner, M., Fontaneto, D., Dumont, H. J., & Kotov,
#' A. A. (2021). Mitogenomics of Cladocera (Branchiopoda): marked gene order rearrangements and
#' independent predation roots. *Molecular Phylogenetics and Evolution*, 164, 107275.
#' \doi{10.1016/j.ympev.2021.107275}
"daphniids"

#' Anurans of urban ponds in two Brazilian cities
#'
#' Counts of 44 species of anurans in 40 ponds, 20 in Campo Grande (Mato Grosso do Sul) and 20
#' in São José dos Campos (São Paulo), from the fieldwork of L. B. C. Menezes for chapter 1 of
#' an MSc dissertation (L. B. C. Menezes, M. R. Severgnini & D. B. Provete, unpublished), with
#' the regional pool of each city. The script is `data-raw/anurans.R`. The counts come from
#' surveys of calling and visual encounter, so their totals depend on effort and detection, and
#' `vignette("abc")` uses incidence only.
#'
#' The hydrology score is the first principal component of log(x + 1) area, depth, floating
#' vegetation and water temperature, standardised across the 40 ponds and rescaled to 0–100,
#' with larger and deeper ponds at higher scores.
#'
#' The regional pool of each city has the anurans that breed in lentic water and were recorded
#' both in its ecoregion (Dinerstein et al. 2017: Cerrado for Campo Grande; Alto Paraná Atlantic
#' forests and Serra do Mar coastal forests for São José dos Campos) and within 200 km of its
#' ponds, in occurrences of GBIF (gbif.org, accessed on 11 November 2026, through the occurrence
#' search API; queries in `data-raw/anurans_pool`), plus the species recorded in its ponds. A
#' species breeds in lentic water if one of its reproductive modes in Nunes-de-Almeida et al.
#' (2021), or in da Silva et al. (2012) for species missing from it, has exotrophic tadpoles in
#' lentic water, or, for species in neither, if at least three congeners are classified and
#' all agree. GBIF records of *Leptodactylus latrans* were taken as *L. luctator*, split from it
#' by Magalhães et al. (2020), and the introduced *Lithobates catesbeianus* was left out.
#'
#' The tree is the time-calibrated tree of Portik et al. (2023) pruned to the species of the
#' pools, with branch lengths in millions of years and the names of GBIF (so *Pithecopus* and
#' *Ololygon* instead of *Phyllomedusa* and *Scinax*; Duellman et al. 2016). The tip of
#' *L. latrans* stands for *L. luctator*. The six species missing from it (`grafted`) were
#' joined at the crown of their genus, or of the species of *Ololygon* in the tree, and the
#' polytomies were resolved at random with branches of length zero.
#'
#' @format List with `comm`, a matrix of ponds by species with counts; `env`, a matrix with
#'   the hydrology score, area (m²), depth (cm), floating vegetation (%), water temperature (°C)
#'   and the urban cover within 1 km in 2022 (%); `city`, a factor; `coords`, a matrix with the
#'   longitude and latitude of each pond in decimal degrees; `pool`, a 0/1 matrix of ponds by
#'   the 155 species of the regional pools; `tree`, a `phylo` object with those species; and
#'   `grafted`, the species added to the tree.
#' @source GBIF.org (11 November 2026), occurrence search API; see `data-raw/anurans_pool`.
#' @references
#' da Silva, F. R., Almeida-Neto, M., do Prado, V. H. M., Haddad, C. F. B., & Rossa-Feres, D. C.
#' (2012). Humidity levels drive reproductive modes and phylogenetic diversity of amphibians in
#' the Brazilian Atlantic Forest. *Journal of Biogeography*, 39, 1720–1732.
#' \doi{10.1111/j.1365-2699.2012.02726.x}
#'
#' Dinerstein, E., Olson, D., Joshi, A., Vynne, C., Burgess, N. D., Wikramanayake, E., et al.
#' (2017). An ecoregion-based approach to protecting half the terrestrial realm. *BioScience*,
#' 67, 534–545. \doi{10.1093/biosci/bix014}
#'
#' Duellman, W. E., Marion, A. B., & Hedges, S. B. (2016). Phylogenetics, classification, and
#' biogeography of the treefrogs (Amphibia: Anura: Arboranae). *Zootaxa*, 4104.
#' \doi{10.11646/zootaxa.4104.1.1}
#'
#' Magalhães, F. de M., Lyra, M. L., de Carvalho, T. R., Baldo, D., Brusquetti, F., Burella,
#' P., Colli, G. R., Gehara, M. C., Giaretta, A. A., Haddad, C. F. B., Langone, J. A.,
#' López, J. A., Napoli, M. F., Santana, D. J., de Sá, R. O., & Garda, A. A. (2020). Taxonomic
#' review of South American butter frogs: phylogeny, geographic patterns, and species
#' delimitation in the *Leptodactylus latrans* species group (Anura: Leptodactylidae).
#' *Herpetological Monographs*, 34, 131–177. \doi{10.1655/0733-1347-31.4.131}
#'
#' Nunes-de-Almeida, C. H. L., Haddad, C. F. B., & Toledo, L. F. (2021). A revised
#' classification of the amphibian reproductive modes. *Salamandra*, 57, 413–427.
#' \doi{10.5281/zenodo.5205572}
#'
#' Portik, D. M., Streicher, J. W., & Wiens, J. J. (2023). Frog phylogeny: a time-calibrated,
#' species-level tree based on hundreds of loci and 5,242 species. *Molecular Phylogenetics and
#' Evolution*, 188, 107907. \doi{10.1016/j.ympev.2023.107907}
"anurans"
