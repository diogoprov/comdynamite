suppressPackageStartupMessages({library(sf); library(jsonlite)})
sf_use_s2(FALSE)
load("/Users/diogoborgesprovete/GitHub/comdynamite/.claude/worktrees/pcps-sim-package-continuation-5c5aa1/data/anurans.rda")
eco <- st_read("out/ecoregions_full.gpkg", quiet = TRUE)
# WKT for the URL: drop parts < 500 km2, simplify 0.1 deg, CCW rings, 3 decimals
mkwkt <- function(g, tol = 0.1) {
  p <- st_cast(st_sfc(g, crs = 4326), "POLYGON")
  a <- as.numeric(st_area(st_transform(p, "ESRI:54034"))) / 1e6
  p <- p[a > 500]
  p <- st_simplify(p, dTolerance = tol, preserveTopology = TRUE)
  p <- st_make_valid(p)
  p <- st_cast(st_cast(st_union(p), "MULTIPOLYGON"), "POLYGON")
  # drop holes (small enclaves of other ecoregions) smaller than 2000 km2? keep holes: they matter
  p <- st_sfc(lapply(p, function(x) st_polygon(lapply(x, function(r) round(r, 3)))), crs = 4326)
  p <- st_make_valid(p)
  p <- st_cast(st_cast(st_union(p), "MULTIPOLYGON"), "POLYGON")
  p <- st_sfc(lapply(p, function(x) {  # enforce CCW exterior, CW holes
    st_polygon(lapply(seq_along(x), function(i) {
      r <- x[[i]]; n <- nrow(r)
      sa <- sum(r[-n, 1] * r[-1, 2] - r[-1, 1] * r[-n, 2]) / 2
      if ((i == 1 && sa < 0) || (i > 1 && sa > 0)) r[n:1, ] else r
    }))
  }), crs = 4326)
  st_union(p) |> st_cast("MULTIPOLYGON")
}
base <- "https://api.gbif.org/v1/occurrence/search/predicate"
dir.create("out/gbif_queries", showWarnings = FALSE)
res <- list(); wkts <- list(); tot <- c()
for (i in seq_len(nrow(eco))) {
  w <- mkwkt(st_geometry(eco)[[i]], tol = 0.05)
  wkts[[i]] <- w
  txt <- gsub(", ", ",", st_as_text(w, digits = 7))
  cat(eco$eco_name[i], "WKT chars:", nchar(txt), " parts:", length(st_cast(w, "POLYGON")), "\n")
  body <- list(predicate = list(type = "and", predicates = list(
      list(type = "equals", key = "TAXON_KEY", value = "952"),
      list(type = "equals", key = "HAS_COORDINATE", value = "true"),
      list(type = "equals", key = "HAS_GEOSPATIAL_ISSUE", value = "false"),
      list(type = "equals", key = "OCCURRENCE_STATUS", value = "PRESENT"),
      list(type = "in", key = "BASIS_OF_RECORD", values = c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE")),
      list(type = "within", geometry = txt))),
    limit = 0, facets = list("SPECIES_KEY"), facetLimit = 2000)
  f <- sprintf("out/gbif_queries/query_eco%s.json", eco$eco_id[i])
  writeLines(toJSON(body, auto_unbox = TRUE), f)
  r <- httr::POST(base, body = httr::upload_file(f, type = "application/json"))
  stopifnot(httr::status_code(r) %in% c(200, 201))
  r <- fromJSON(httr::content(r, "text", encoding = "UTF-8"))
  cat("  total records:", r$count, "\n"); tot[eco$eco_name[i]] <- r$count
  out <- r$facets$counts[[1]]
  stopifnot(nrow(out) < 2000)
  out$eco_name <- eco$eco_name[i]; out$eco_id <- eco$eco_id[i]
  res[[i]] <- out
}
res <- do.call(rbind, res); names(res)[1:2] <- c("speciesKey", "n_records")
qg <- st_sf(eco_name = eco$eco_name, eco_id = eco$eco_id, geometry = do.call(c, wkts))
st_write(qg, "out/ecoregions_gbif_query_polygons.gpkg", delete_dsn = TRUE, quiet = TRUE)
writeLines(c(paste("Access date:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "Endpoint: POST https://api.gbif.org/v1/occurrence/search/predicate (Content-Type: application/json)",
  "Bodies: out/gbif_queries/query_eco<eco_id>.json",
  "Reproduce: curl -X POST -H 'Content-Type: application/json' -d @query_eco567.json https://api.gbif.org/v1/occurrence/search/predicate",
  paste("Total records:", paste(names(tot), tot, collapse = "; "))), "out/gbif_query_info.txt")
# names
keys <- unique(res$speciesKey)
nm <- do.call(rbind, lapply(keys, function(k) {
  s <- fromJSON(paste0("https://api.gbif.org/v1/species/", k))
  data.frame(speciesKey = k, gbif_name = s$canonicalName %||% NA, status = s$taxonomicStatus %||% NA,
             rank = s$rank %||% NA, family = s$family %||% NA, genus = s$genus %||% NA)
}))
res <- merge(res, nm)
write.csv(res, "out/gbif_facets_raw.csv", row.names = FALSE)
print(table(res$eco_name)); print(table(res$status, res$rank))
for (e in unique(res$eco_name)) { x <- res$n_records[res$eco_name == e]
  cat(e, ": n spp", length(x), "| quantiles", quantile(x, c(0, .1, .25, .5, .75, .9, 1)),
      "| >=1/3/5/10/20:", sapply(c(1, 3, 5, 10, 20), function(t) sum(x >= t)), "\n") }
