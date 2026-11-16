# Agreement between the species lists of Benicio et al. (2021) and ATLANTIC AMPHIBIANS (Vancine et
# al. 2018) in the same protected areas, to set the false-absence rate of the simulations of
# imperfect detection in the DAMOCLES paper (dev/plano-artigo-damocles.md, S2 and S6).
# - Benicio et al. sampled the 14 areas themselves with a standard protocol (December 2014 to
#   February 2017, 27 days per area, six breeding sites, calls, tadpoles and visual search).
# - Sites of ATLANTIC AMPHIBIANS within CUTOFF km of each area (default 2 and 15). Areas whose
#   records come from the same project are left out (circular): PEVA (Benicio & Da Silva 2017)
#   and PECU, PESB and PESV (Da Silva et al. 2017a).
# - Names: the canonical names of data-raw/benicio2021_taxonomy.csv (matrix) and
#   data-raw/benicio2021_regional_pool.csv (ATLANTIC AMPHIBIANS in Sao Paulo); Ischnocnema
#   guentheri and the lineages Ischnocnema_sp1 to sp3 count as one taxon (the complex), and
#   records of ATLANTIC AMPHIBIANS without a valid name are dropped.
# - For each area: species in both lists, only in Benicio et al., only in ATLANTIC AMPHIBIANS.
#   The share of the union missing from each list is a disagreement, not an absolute omission:
#   the sites are not identical, the periods differ (real turnover counts as disagreement) and
#   both lists miss species.
# Data: DATA_DIR/benicio2021/ (Dryad) and ATLANTIC_DIR as in data-raw/benicio2021_trees.R. Output:
# data-raw/benicio2021_omission.rds.
data_dir <- Sys.getenv("DATA_DIR", "data-raw")
atlantic_dir <- Sys.getenv("ATLANTIC_DIR", file.path(data_dir, "atlantic_amphibians"))
cutoffs <- as.numeric(strsplit(Sys.getenv("CUTOFF", "2,15"), ",")[[1]])
norm <- function(x) sub("_+$", "", gsub("'", "", x))

comm <- read.csv(file.path(data_dir, "benicio2021", "frog_data.csv"), row.names = 1,
                 check.names = FALSE)
colnames(comm) <- norm(colnames(comm))
coords <- read.csv(file.path(data_dir, "benicio2021", "coordinates.csv"), row.names = 1)
taxonomy <- read.csv("data-raw/benicio2021_taxonomy.csv")
pool <- read.csv("data-raw/benicio2021_regional_pool.csv")
sites <- read.csv(file.path(atlantic_dir, "ATLANTIC_AMPHIBIANS_sites.csv"))
species <- read.csv(file.path(atlantic_dir, "ATLANTIC_AMPHIBIANS_species.csv"))
refs <- read.csv(file.path(atlantic_dir, "ATLANTIC_AMPHIBIANS_references.csv"))
species <- species[species$order %in% "Anura" & !is.na(species$valid_name), ]
species$valid_name <- gsub(" ", "_", species$valid_name)

complex <- "Ischnocnema_guentheri_complex"
lump <- function(x) ifelse(x %in% c("Ischnocnema_guentheri", paste0("Ischnocnema_sp", 1:3)),
                           complex, x)
comm_name <- taxonomy$name[match(colnames(comm), taxonomy$comm_name)]
stopifnot(!anyNA(comm_name))
colnames(comm) <- lump(comm_name)
atlantic_name <- pool$name[match(species$valid_name, pool$atlantic_name)]
species$name <- lump(ifelse(is.na(atlantic_name), species$valid_name, atlantic_name))

same_project <- c("PEVA", "PECU", "PESB", "PESV")
km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  6371 * 2 * asin(sqrt(a))
}
compare <- function(area, cutoff) {
  d <- km(coords[area, "Latitude"], coords[area, "Longitude"], sites$latitude, sites$longitude)
  ids <- sites$id[d <= cutoff]
  b <- unique(colnames(comm)[comm[area, ] > 0])
  a <- unique(species$name[species$id %in% ids])
  # Names of ATLANTIC AMPHIBIANS outside the regional pool table (outside Sao Paulo or dropped).
  data.frame(area = area, cutoff = cutoff, n_sites = length(ids),
             studies = length(unique(sites$reference_number[sites$id %in% ids])),
             years = if (length(ids)) paste(range(c(sites$year_start[sites$id %in% ids],
                                                    sites$year_finish[sites$id %in% ids]),
                                                  na.rm = TRUE), collapse = "-") else NA,
             benicio = length(b), atlantic = length(a), both = length(intersect(a, b)),
             only_benicio = length(setdiff(b, a)), only_atlantic = length(setdiff(a, b)),
             missing_benicio = length(setdiff(a, b)) / length(union(a, b)),
             missing_atlantic = length(setdiff(b, a)) / length(union(a, b)),
             jaccard = length(intersect(a, b)) / length(union(a, b)),
             only_benicio_names = paste(sort(setdiff(b, a)), collapse = ";"),
             only_atlantic_names = paste(sort(setdiff(a, b)), collapse = ";"))
}
areas <- setdiff(rownames(comm), same_project)
res <- do.call(rbind, lapply(cutoffs, function(k) do.call(rbind, lapply(areas, compare, cutoff = k))))
res <- res[res$n_sites > 0, ]
options(width = 200)
print(res[, c("area", "cutoff", "n_sites", "studies", "years", "benicio", "atlantic", "both",
              "only_benicio", "only_atlantic", "missing_benicio", "missing_atlantic",
              "jaccard")], row.names = FALSE, digits = 2)
pooled <- do.call(rbind, lapply(split(res, res$cutoff), function(r) data.frame(
  cutoff = r$cutoff[1], areas = nrow(r),
  missing_benicio = sum(r$only_atlantic) / sum(r$both + r$only_benicio + r$only_atlantic),
  missing_atlantic = sum(r$only_benicio) / sum(r$both + r$only_benicio + r$only_atlantic),
  median_missing_atlantic = stats::median(r$missing_atlantic),
  median_missing_benicio = stats::median(r$missing_benicio))))
cat("\nPooled over areas (share of the union of the two lists):\n")
print(pooled, row.names = FALSE, digits = 2)
# Species most often only in one list (2 km).
r2 <- res[res$cutoff == min(cutoffs), ]
tab <- function(x) sort(table(unlist(strsplit(x[nzchar(x)], ";"))), decreasing = TRUE)
cat("\nOnly in Benicio et al. (", min(cutoffs), " km), most frequent:\n", sep = "")
print(utils::head(tab(r2$only_benicio_names), 12))
cat("\nOnly in ATLANTIC AMPHIBIANS (", min(cutoffs), " km), most frequent:\n", sep = "")
print(utils::head(tab(r2$only_atlantic_names), 12))
saveRDS(list(results = res, pooled = pooled), "data-raw/benicio2021_omission.rds")
