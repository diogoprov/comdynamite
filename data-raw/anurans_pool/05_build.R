suppressPackageStartupMessages({library(readxl); library(ape); library(jsonlite); library(sf)})
sf_use_s2(FALSE)
load("/Users/diogoborgesprovete/GitHub/comdynamite/.claude/worktrees/pcps-sim-package-continuation-5c5aa1/data/anurans.rda")
thr <- 3
g <- read.csv("out/gbif_facets_raw.csv"); g$sp <- gsub(" ", "_", g$gbif_name)
syn <- read.csv("out/gbif_synonyms.csv"); syn$synonym <- gsub(" ", "_", syn$synonym)
obm <- read.csv("out/observed_gbif_match.csv")
tr <- read.tree("/Users/diogoborgesprovete/Library/CloudStorage/Box-Box/Meus Artigos/2026-Leandro Cap 1/data/raw/legacy/TreePL-Rooted_Anura_bestTree.tre")
tips <- tr$tip.label; tipgen <- unique(sub("_.*", "", tips))

## ---- local (200 km) counts via GET (exact URLs recorded)
cen <- aggregate(anurans$coords, list(city = anurans$city), mean)
get_urls <- c(); loc <- list()
for (i in 1:2) {
  p <- st_sfc(st_point(c(cen$lon[i], cen$lat[i])), crs = 4326) |> st_transform(5880) |> st_buffer(2e5, nQuadSegs = 6) |> st_transform(4326)
  m <- round(st_coordinates(p)[, 1:2], 3); n <- nrow(m)
  if (sum(m[-n, 1] * m[-1, 2] - m[-1, 1] * m[-n, 2]) < 0) m <- m[n:1, ]
  w <- paste0("POLYGON((", paste(m[, 1], m[, 2], collapse = ","), "))")
  u <- paste0("https://api.gbif.org/v1/occurrence/search?taxonKey=952&hasCoordinate=true&hasGeospatialIssue=false&occurrenceStatus=PRESENT",
    "&basisOfRecord=PRESERVED_SPECIMEN&basisOfRecord=HUMAN_OBSERVATION&basisOfRecord=MATERIAL_SAMPLE&facet=speciesKey&facetLimit=1000&limit=0&geometry=",
    URLencode(w, reserved = TRUE))
  r <- fromJSON(u); f <- r$facets$counts[[1]]
  cat(as.character(cen$city[i]), "200 km: records", r$count, "species", nrow(f), "\n")
  loc[[as.character(cen$city[i])]] <- setNames(f$count, f$name); get_urls[as.character(cen$city[i])] <- u
}
writeLines(c(readLines("out/gbif_query_info.txt")[1:5], "", paste("Local 200-km buffer queries (GET), access", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
             paste0("## ", names(get_urls), "\n", get_urls)), "out/gbif_query_info.txt")

## ---- traits
x <- as.data.frame(suppressMessages(read_excel("downloads/nunes2021/Nunes-de-Almeida_SupplementaryTable.xlsx", col_names = FALSE, .name_repair = "minimal", skip = 3)))
x <- x[!is.na(x[[1]]) & x[[1]] == "Anura", ]; rm <- as.matrix(x[, 4:77]); storage.mode(rm) <- "numeric"; rm[is.na(rm)] <- 0
defs <- read.csv("out/nunes2021_mode_definitions.csv"); defs <- defs[!duplicated(defs$RM), ]
ll <- grepl("larvae develop in lentic water", defs$definition); cont <- grepl("water-filled container", defs$definition); exo <- grepl(", exotrophic", defs$definition)
pondRM <- defs$RM[ll & !cont & exo]
nu <- data.frame(name = gsub(" ", "_", trimws(x[[3]])), genus = sub(" .*", "", trimws(x[[3]])),
                 pond = rowSums(rm[, pondRM, drop = FALSE]) > 0, modes = apply(rm, 1, function(r) paste(which(r > 0), collapse = ";")))
ds <- as.data.frame(suppressMessages(read_excel("downloads/dasilva2012/Data_amphibian.xls")))[, 1:2]
names(ds) <- c("name", "MR"); ds$name <- gsub(" ", "_", trimws(ds$name)); ds$mode <- as.integer(sub("mode ", "", ds$MR))
hp_T <- c(1, 10, 11, 13, 14, 24, 28, 30); hp_NA <- c(4, 12, 17, 18, 29, 33)
ds$pond <- ifelse(ds$mode %in% hp_T, TRUE, ifelse(ds$mode %in% hp_NA, NA, FALSE))
cc <- merge(ds, nu, by = "name"); cat("da Silva vs Nunes concordance (non-NA):", mean(cc$pond.x == cc$pond.y, na.rm = TRUE), "n =", sum(!is.na(cc$pond.x)), "\n")

## ---- species table
manual_tree <- c(Leptodactylus_luctator = "Leptodactylus_latrans", Pithecopus_azureus = "Phyllomedusa_azurea")
cat("manual tree labels present:", manual_tree %in% tips, "\n")
cat(grep("Pithecopus_az|Phyllomedusa_az|Leptodactylus_latr|Leptodactylus_luct|Brachycephalus_(roten|ephipp)", tips, value = TRUE), "\n")
sp <- unique(g[, c("speciesKey", "sp", "family", "genus")])
alln <- function(k, s) unique(c(s, syn$synonym[syn$speciesKey == k]))
sp$tree_label <- NA; sp$nunes <- NA; sp$nunes_modes <- NA; sp$nunes_name <- NA; sp$ds <- NA; sp$ds_mode <- NA
for (i in seq_len(nrow(sp))) {
  nn <- alln(sp$speciesKey[i], sp$sp[i])
  t <- intersect(nn, tips); if (!length(t) && sp$sp[i] %in% names(manual_tree) && manual_tree[sp$sp[i]] %in% tips) t <- manual_tree[[sp$sp[i]]]
  if (length(t)) sp$tree_label[i] <- t[1]
  j <- match(nn, nu$name); j <- j[!is.na(j)]
  if (length(j)) { sp$nunes[i] <- any(nu$pond[j]); sp$nunes_modes[i] <- paste(nu$modes[j], collapse = ";"); sp$nunes_name[i] <- nu$name[j[1]] }
  j <- match(nn, ds$name); j <- j[!is.na(j)]
  if (length(j)) { sp$ds[i] <- ds$pond[j[1]]; sp$ds_mode[i] <- ds$mode[j[1]] }
}
# add observed species missing from GBIF lists
for (o in setdiff(colnames(anurans$comm), sp$sp)) {
  j <- match(o, nu$name)
  sp <- rbind(sp, data.frame(speciesKey = NA, sp = o, family = NA, genus = sub("_.*", "", o), tree_label = if (o %in% tips) o else NA,
    nunes = if (!is.na(j)) nu$pond[j] else NA, nunes_modes = if (!is.na(j)) nu$modes[j] else NA, nunes_name = NA, ds = NA, ds_mode = NA))
}
sp$lentic <- ifelse(!is.na(sp$nunes), sp$nunes, sp$ds)
sp$lentic_source <- ifelse(!is.na(sp$nunes), paste0("Nunes-de-Almeida2021 RM ", sp$nunes_modes),
                    ifelse(!is.na(sp$ds), paste0("daSilva2012 H&P mode ", sp$ds_mode), NA))
# genus-level inference: all classified congeners (Nunes global list + pool species classified at species level)
gen_tab <- rbind(data.frame(genus = nu$genus, pond = nu$pond, id = nu$name),
                 data.frame(genus = sp$genus, pond = sp$lentic, id = ifelse(is.na(sp$nunes_name), sp$sp, sp$nunes_name))[!is.na(sp$lentic), ])
gen_tab <- gen_tab[!duplicated(gen_tab$id), ]
ga <- aggregate(pond ~ genus, gen_tab, function(v) c(n = length(v), k = sum(v)))
ga <- data.frame(genus = ga$genus, n = ga$pond[, "n"], k = ga$pond[, "k"])
sp$lentic_inferred <- sp$lentic
for (i in which(is.na(sp$lentic))) {
  r <- ga[ga$genus == sp$genus[i], ]
  if (nrow(r) && r$n >= 3 && (r$k == r$n || r$k == 0)) {
    sp$lentic_inferred[i] <- r$k == r$n
    sp$lentic_source[i] <- sprintf("genus-inferred (%d/%d classified congeners lentic)", r$k, r$n)
  } else sp$lentic_source[i] <- if (nrow(r)) sprintf("unknown (genus mixed or <3: %d/%d congeners lentic)", r$k, r$n) else "unknown (no classified congener)"
}
sp$in_tree <- !is.na(sp$tree_label)
sp$congener_in_tree <- sp$genus %in% tipgen | sub("_.*", "", ifelse(is.na(sp$tree_label), "zz", sp$tree_label)) %in% tipgen

## ---- per-city pools
city_eco <- list("Campo Grande" = "Cerrado", "São José dos Campos" = c("Alto Paraná Atlantic forests", "Serra do Mar coastal forests"))
out <- list()
for (ci in names(city_eco)) {
  ge <- g[g$eco_name %in% city_eco[[ci]], ]
  n <- tapply(ge$n_records, ge$sp, sum)
  by_eco <- sapply(city_eco[[ci]], function(e) { v <- setNames(ge$n_records[ge$eco_name == e], ge$sp[ge$eco_name == e]); v })
  obs <- colnames(anurans$comm)[colSums(anurans$comm[anurans$city == ci, , drop = FALSE]) > 0]
  spp <- union(names(n), obs)
  d <- sp[match(spp, sp$sp), ]
  d$city <- ci; d$n_records <- ifelse(is.na(n[spp]), 0, n[spp]); d$observed <- spp %in% obs
  d$n_by_ecoregion <- sapply(spp, function(s) paste(sapply(city_eco[[ci]], function(e) { v <- ge$n_records[ge$eco_name == e & ge$sp == s]; paste0(e, "=", if (length(v)) v else 0) }), collapse = "; "))
  d$n_records_200km <- { v <- loc[[ci]][as.character(d$speciesKey)]; ifelse(is.na(v), 0, v) }
  d$passes_threshold <- d$n_records >= thr
  d$in_pool_strict <- d$observed | (d$passes_threshold & d$lentic %in% TRUE)
  d$in_pool_inferred <- d$observed | (d$passes_threshold & d$lentic_inferred %in% TRUE)
  d$notes <- paste0(d$lentic_source, ifelse(!is.na(d$tree_label) & d$tree_label != d$sp, paste0("; tree tip = ", d$tree_label), ""),
                    ifelse(d$observed & d$n_records == 0, "; observed but no GBIF record in ecoregion query", ""),
                    ifelse(d$observed & d$n_records > 0 & d$n_records < thr, "; observed but below record threshold", ""),
                    ifelse(d$observed & !(d$lentic_inferred %in% TRUE), "; observed in ponds though not classified as lentic breeder", ""))
  out[[ci]] <- d
}
P <- do.call(rbind, out)
P <- data.frame(city = P$city, species = P$sp, tree_label = P$tree_label, family = P$family, gbif_speciesKey = P$speciesKey,
  n_records = P$n_records, n_by_ecoregion = P$n_by_ecoregion, n_records_200km = P$n_records_200km, passes_threshold = P$passes_threshold,
  lentic = P$lentic, lentic_incl_genus_inferred = P$lentic_inferred, in_tree = P$in_tree, congener_in_tree = P$congener_in_tree,
  observed = P$observed, in_pool_strict = P$in_pool_strict, in_pool_inferred = P$in_pool_inferred, source_notes = P$notes)
P <- P[order(P$city, -P$in_pool_inferred, -P$n_records), ]
write.csv(P, "pools.csv", row.names = FALSE)

## ---- summaries
for (ci in names(city_eco)) {
  d <- P[P$city == ci, ]; k <- d[d$passes_threshold | d$observed, ]
  cat("\n=====", ci, "\n all GBIF spp:", sum(d$n_records > 0), "| >=", thr, ":", sum(d$passes_threshold), "| observed:", sum(d$observed), "\n")
  cat(" lentic (species-level) among >=thr:"); print(table(d$lentic[d$passes_threshold], useNA = "always"))
  cat(" lentic (incl. genus-inferred) among >=thr:"); print(table(d$lentic_incl_genus_inferred[d$passes_threshold], useNA = "always"))
  for (v in c("in_pool_strict", "in_pool_inferred")) { q <- d[d[[v]], ]
    cat(sprintf(" %s: %d spp | in tree %d | not in tree but congener %d | no congener %d | with >=1 rec within 200 km: %d (>=3: %d)\n", v, nrow(q), sum(q$in_tree),
        sum(!q$in_tree & q$congener_in_tree), sum(!q$in_tree & !q$congener_in_tree), sum(q$n_records_200km > 0 | q$observed), sum(q$n_records_200km >= 3 | q$observed))) }
  cat(" observed flags:\n"); print(d[d$observed & (!d$passes_threshold | !(d$lentic_incl_genus_inferred %in% TRUE) | !d$in_tree), c("species", "n_records", "lentic", "lentic_incl_genus_inferred", "in_tree", "source_notes")], row.names = FALSE)
  q <- d[d$in_pool_inferred & !d$in_tree, ]; cat(" pool spp not in tree:", paste(q$species, collapse = ", "), "\n")
}
