# Daphniids of the EMAP Northeastern Lakes Survey, as in Leibold et al. (2010). The tables
# are public data of the US Environmental Protection Agency. The tree is the topology of
# Figure 3a. The biogeographic contrasts follow the text, with rivers approximated by
# longitude because the article does not give the assignment of each lake.
base <- "https://archive.epa.gov/emap/archive-emap/web/txt/"
read_emap <- function(file) {
  lines <- readLines(paste0(base, file, ".txt"), warn = FALSE)
  out <- read.csv(text = lines[grep('^"', lines)[1]:length(lines)], stringsAsFactors = FALSE)
  out[] <- lapply(out, function(x) if (is.character(x)) trimws(x) else x)
  out
}
counts <- read_emap("zoocnt")
chemistry <- read_emap("chmval-2")
design <- read_emap("ldesign")
watershed <- read_emap("watchar")

# first visit of each lake
visits <- unique(counts[counts$SAMPLED == "Yes", c("LAKE_ID", "YEAR", "VISIT_NO")])
visits <- visits[order(visits$YEAR, visits$VISIT_NO), ]
visits <- visits[!duplicated(visits$LAKE_ID), ]
number <- function(x) suppressWarnings(as.numeric(x))
lakes <- merge(visits, chemistry[, c("LAKE_ID", "YEAR", "VISIT_NO", "COLOR", "COND", "DOC", "NTL",
                                     "PHEQ", "PTL")])
# the design has one row per lake and year, and some leave the ecoregion blank
design <- design[order(design$ECO_ATH == ""), ]
design <- design[!duplicated(design$LAKE_ID), ]
lakes <- merge(lakes, design[, c("LAKE_ID", "STATE", "ECO_ATH", "LAT_DD", "LON_DD")])
lakes <- merge(lakes, watershed[, c("LAKE_ID", "AV_DEP", "ELEV")])
variables <- c("COLOR", "COND", "DOC", "NTL", "PHEQ", "PTL", "AV_DEP", "ELEV")
lakes[c(variables, "LAT_DD", "LON_DD")] <- lapply(lakes[c(variables, "LAT_DD", "LON_DD")], number)
lakes$TN_TP <- lakes$NTL / lakes$PTL
lakes <- lakes[complete.cases(lakes[c(variables, "TN_TP")]) & is.finite(lakes$TN_TP), ]
lakes <- lakes[order(lakes$LAKE_ID), ]
rownames(lakes) <- lakes$LAKE_ID

# DAPLONG is D. longiremis in the survey; the only taxon that fits D. longispina of the figure
taxa <- c(DAPRETR = "Daphnia_retrocurva", DAPPARV = "Daphnia_parvula", DAPCATA = "Daphnia_catawba",
          DAPPULI = "Daphnia_pulicaria", DAPPULX = "Daphnia_pulex", DAPAMBI = "Daphnia_ambigua",
          DAPLONG = "Daphnia_longispina", DAPGALM = "Daphnia_galeata_mendotae",
          DAPROSE = "Daphnia_rosea", DAPLAEV = "Daphnia_laevis", DAPDUBI = "Daphnia_dubia",
          SIMOSER = "Simocephalus_serrulatus", SIMOVET = "Simocephalus_vetulus")
comm <- matrix(0L, nrow(lakes), length(taxa), dimnames = list(lakes$LAKE_ID, unname(taxa)))
found <- merge(counts, lakes[, c("LAKE_ID", "YEAR", "VISIT_NO")])
found <- found[found$TAXACODE %in% names(taxa) & number(found$COUNT) > 0, ]
comm[cbind(found$LAKE_ID, taxa[found$TAXACODE])] <- 1L

tree <- ape::read.tree(text = paste0(
  "(((((((Daphnia_retrocurva,Daphnia_parvula),Daphnia_catawba),",
  "(Daphnia_pulicaria,Daphnia_pulex)),Daphnia_ambigua),",
  "(Daphnia_longispina,(Daphnia_galeata_mendotae,Daphnia_rosea))),",
  "(Daphnia_laevis,Daphnia_dubia)),(Simocephalus_serrulatus,Simocephalus_vetulus));"))

east_hudson <- lakes$STATE %in% c("CT", "MA", "ME", "NH", "RI", "VT") |
  (lakes$STATE == "NY" & lakes$LON_DD > -73.8)
east_connecticut <- lakes$STATE %in% c("ME", "NH", "RI") |
  (lakes$STATE %in% c("CT", "MA") & lakes$LON_DD > -72.6)
adirondacks <- lakes$ECO_ATH == "Adirondacks"
maine <- lakes$STATE == "ME"
new_hampshire <- lakes$STATE == "NH"
contrast <- function(one, other) ifelse(one, 1L, ifelse(other, -1L, 0L))
biogeography <- cbind(
  S1 = contrast(east_hudson, !east_hudson),
  S2 = contrast(!east_hudson & adirondacks, !east_hudson & !adirondacks),
  S3 = contrast(east_hudson & east_connecticut, east_hudson & !east_connecticut),
  S4 = contrast(maine, east_hudson & !maine),
  S5 = contrast(new_hampshire, east_hudson & !east_connecticut),
  S6 = contrast(new_hampshire, lakes$STATE %in% c("RI", "MA"))
)
rownames(biogeography) <- lakes$LAKE_ID

env <- as.matrix(lakes[, c(variables, "TN_TP")])
coords <- as.matrix(lakes[, c("LON_DD", "LAT_DD")])
colnames(coords) <- c("lon", "lat")
daphniids <- list(comm = comm, env = env, biogeography = biogeography, tree = tree,
                  coords = coords)
usethis::use_data(daphniids, overwrite = TRUE, compress = "xz")
