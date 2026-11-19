# NAAMP species (>= 50 stops, no complexes) with maximum longevity (union Stark 2018, Stark 2020, TetrapodTraits).
D <- "data-raw"; k <- function(x) gsub(" ", "_", trimws(x))
s18 <- read.csv(file.path(D, "stark2018/geb12804-sup-0001-appendixs1.csv"), check.names = FALSE, fileEncoding = "latin1")
s20 <- suppressWarnings(readxl::read_excel(file.path(D, "stark2020/AppendixS1.xlsx"))); s20 <- subset(s20, Class == "Amphibia")
tt <- read.csv(file.path(D, "tetrapodtraits/TetrapodTraits_v3.0.1.csv")); tt <- subset(tt, Class == "Amphibia")
L <- rbind(data.frame(sp = k(s18$Species), lon = suppressWarnings(as.numeric(s18[["Maximum Longevity (in years)"]]))),
           data.frame(sp = k(tt$Scientific.Name), lon = tt$MaxLongevity),
           data.frame(sp = k(s20$Species), lon = suppressWarnings(as.numeric(s20[["Maximum Longevity (in years)"]]))))
L <- aggregate(lon ~ sp, L[!is.na(L$lon), ], max)
runs <- read.csv(file.path(D, "naamp/Runs.csv")); names(runs) <- trimws(names(runs)); runs <- subset(runs, UnifiedProtocol == 1)
st <- merge(read.csv(file.path(D, "naamp/Stops.csv"))[, c("RunID", "StopNumber", "SiteID", "SkippedStop")],
            runs[, c("RunID", "RouteNumber", "State")]); st <- st[st$SkippedStop == 0, ]
ct <- read.csv(file.path(D, "naamp/Counts.csv"))
c2 <- merge(ct[ct$CallingIndex > 0, c("RunID", "StopNumber", "Species")], st)
sp <- aggregate(SiteID ~ Species, c2, function(v) length(unique(v))); names(sp)[2] <- "stops"
sp$routes <- aggregate(RouteNumber ~ Species, c2, function(v) length(unique(v)))$RouteNumber
sp$states <- aggregate(State ~ Species, c2, function(v) length(unique(v)))$State
sp <- subset(sp, stops >= 50 & !grepl("/| sp|x ", Species))
key <- k(sp$Species); alt1 <- sub("^Lithobates_", "Rana_", key)
sp$lon <- L$lon[match(key, L$sp)]; m <- is.na(sp$lon); sp$lon[m] <- L$lon[match(alt1[m], L$sp)]
sp <- sp[order(sp$lon), ]
saveRDS(sp, "data-raw/persistence_proxy/species_table.rds")
print(sp, row.names = FALSE)
