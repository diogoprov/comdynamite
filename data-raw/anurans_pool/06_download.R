# GBIF download with a DOI for the regional pools of `anurans`. The species facets of
# 02_gbif.R came from the occurrence search API, which gives no DOI; this script requests one
# download with the same filters and the three ecoregion polygons of query_eco*.json joined by
# OR, so that the records behind pools.csv can be cited. It needs a GBIF account, read from
# GBIF_USER, GBIF_PWD and GBIF_EMAIL in the environment (for example in ~/.Renviron); the
# credentials are never printed. Run from the root of the package:
# Rscript data-raw/anurans_pool/06_download.R. The download key and the request are written to
# data-raw/anurans_pool/gbif_download.txt; the DOI is added to it when GBIF has prepared the file.
dir <- "data-raw/anurans_pool"
queries <- file.path(dir, paste0("query_eco", c(567, 500, 439), ".json"))
predicates <- lapply(queries, function(f) jsonlite::fromJSON(f, simplifyVector = FALSE)$predicate)
filters <- predicates[[1]]$predicates[1:5]
stopifnot(all(vapply(predicates, function(p) identical(p$predicates[1:5], filters), TRUE)))
polygons <- lapply(predicates, function(p) p$predicates[[6]])
body <- list(creator = Sys.getenv("GBIF_USER"),
             notificationAddresses = list(Sys.getenv("GBIF_EMAIL")),
             sendNotification = TRUE, format = "SIMPLE_CSV",
             predicate = list(type = "and",
                              predicates = c(filters, list(list(type = "or",
                                                                predicates = polygons)))))
if (!all(nzchar(Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))))) {
  stop("GBIF_USER, GBIF_PWD and GBIF_EMAIL must be set", call. = FALSE)
}
request <- file.path(tempdir(), "gbif_request.json")
writeLines(jsonlite::toJSON(body, auto_unbox = TRUE), request)
response <- httr::POST("https://api.gbif.org/v1/occurrence/download/request",
                       httr::authenticate(Sys.getenv("GBIF_USER"), Sys.getenv("GBIF_PWD")),
                       body = httr::upload_file(request, type = "application/json"))
if (!httr::status_code(response) %in% c(200, 201)) {
  stop("GBIF answered ", httr::status_code(response), call. = FALSE)
}
key <- trimws(httr::content(response, "text", encoding = "UTF-8"))
writeLines(c(paste("Download key:", key),
             paste("Requested:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
             paste0("Status: https://api.gbif.org/v1/occurrence/download/", key),
             "Predicate: the five filters of query_eco*.json and their three polygons joined by OR"),
           file.path(dir, "gbif_download.txt"))
cat("Download key:", key, "\n")
