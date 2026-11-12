suppressPackageStartupMessages(library(jsonlite))
g <- read.csv("out/gbif_facets_raw.csv")
keys <- unique(g$speciesKey)
syn <- do.call(rbind, lapply(keys, function(k) {
  r <- tryCatch(fromJSON(paste0("https://api.gbif.org/v1/species/", k, "/synonyms?limit=300")), error = function(e) NULL)
  s <- if (!is.null(r) && length(r$results)) r$results$canonicalName else character()
  s <- unique(s[!is.na(s) & lengths(strsplit(s, " ")) == 2])
  if (length(s)) data.frame(speciesKey = k, synonym = s) else NULL
}))
write.csv(syn, "out/gbif_synonyms.csv", row.names = FALSE)
# observed species: match to backbone
load("/Users/diogoborgesprovete/GitHub/comdynamite/.claude/worktrees/pcps-sim-package-continuation-5c5aa1/data/anurans.rda")
ob <- do.call(rbind, lapply(colnames(anurans$comm), function(n) {
  m <- fromJSON(paste0("https://api.gbif.org/v1/species/match?kingdom=Animalia&order=Anura&name=", URLencode(gsub("_", " ", n))))
  data.frame(observed = n, usageKey = m$usageKey %||% NA, matched = m$canonicalName %||% NA, status = m$status %||% NA,
             matchType = m$matchType %||% NA, acceptedKey = m$acceptedUsageKey %||% NA, speciesKey = m$speciesKey %||% NA, species = m$species %||% NA)
}))
write.csv(ob, "out/observed_gbif_match.csv", row.names = FALSE)
print(ob[ob$status != "ACCEPTED" | ob$matchType != "EXACT" | gsub("_", " ", ob$observed) != ob$species, ])
cat("n synonyms:", nrow(syn), "\n")
