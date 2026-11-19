D <- "data-raw"
s18 <- read.csv(file.path(D, "stark2018/geb12804-sup-0001-appendixs1.csv"), check.names = FALSE, fileEncoding = "latin1")
s20 <- suppressWarnings(readxl::read_excel(file.path(D, "stark2020/AppendixS1.xlsx")))
s20 <- subset(s20, Class == "Amphibia")
tt <- read.csv(file.path(D, "tetrapodtraits/TetrapodTraits_v3.0.1.csv")); tt <- subset(tt, Class == "Amphibia")
k <- function(x) gsub(" ", "_", trimws(x))
lon <- unique(c(k(s18$Species), k(s20$Species), k(tt$Scientific.Name[!is.na(tt$MaxLongevity)])))
cat("amphibians with longevity: Stark2018", nrow(s18), " Stark2020", nrow(s20), " TT", sum(!is.na(tt$MaxLongevity)), " union", length(lon), "\n")
cat("Stark2018 data quality:\n"); print(table(s18[["Data quality"]])); print(table(s18[["Origin of data"]]))
pool <- read.csv("data-raw/benicio2021_regional_pool.csv")
syn <- c(Boana = "Hypsiboas"); alt <- sub("^Boana_", "Hypsiboas_", pool$name)
alt <- sub("a$", "us", alt[grepl("^Hypsiboas", alt)])  # rough gender fix
cands <- function(i) unique(c(pool$name[i], pool$portik_name[i], pool$atlantic_name[i],
  sub("^Boana_", "Hypsiboas_", pool$name[i]), sub("a$", "us", sub("^Boana_", "Hypsiboas_", pool$name[i])),
  sub("^Scinax_", "Ololygon_", pool$name[i]), sub("^Phyllomedusa_", "Pithecopus_", pool$name[i]),
  sub("^Physalaemus_nattereri", "Eupemphix_nattereri", pool$name[i])))
hit <- vapply(seq_len(nrow(pool)), function(i) any(cands(i) %in% lon), logical(1))
cat("pool with longevity:", sum(hit), "/", nrow(pool), "; community:", sum(hit & pool$in_community), "/", sum(pool$in_community), "\n")
print(pool$name[hit])
in_s18 <- vapply(seq_len(nrow(pool)), function(i) any(cands(i) %in% k(s18$Species)), logical(1))
cat("of which in Stark2018:", sum(in_s18), "\n")
# NAAMP
runs <- read.csv(file.path(D, "naamp/Runs.csv")); names(runs) <- trimws(names(runs))
stops <- read.csv(file.path(D, "naamp/Stops.csv"))
counts <- read.csv(file.path(D, "naamp/Counts.csv"))
runs_u <- subset(runs, UnifiedProtocol == 1)
cat("\nNAAMP runs:", nrow(runs), " unified:", nrow(runs_u), " years:", range(runs_u$SurveyYear), "\n")
cat("states:", length(unique(runs_u$State)), " routes:", length(unique(runs_u$RouteNumber)), "\n")
st <- merge(stops[stops$SkippedStop == 0, c("RunID", "StopNumber", "SiteID")], runs_u[, c("RunID", "RouteNumber", "SurveyYear")])
cat("stops (SiteID):", length(unique(st$SiteID)), "\n")
yrs <- tapply(st$SurveyYear, st$SiteID, function(v) length(unique(v)))
cat("years surveyed per stop: \n"); print(quantile(yrs, c(.1, .25, .5, .75, .9)))
cat("stops with >= 5 years:", sum(yrs >= 5), " >= 10:", sum(yrs >= 10), "\n")
vis <- aggregate(RunID ~ SiteID + SurveyYear, st, length)
cat("visits per stop-year:\n"); print(table(pmin(vis$RunID, 5)))
c2 <- merge(counts[counts$CallingIndex > 0, c("RunID", "StopNumber", "Species")], st)
sp <- aggregate(cbind(stops = SiteID) ~ Species, c2, function(v) length(unique(v)))
sp$routes <- aggregate(RouteNumber ~ Species, c2, function(v) length(unique(v)))$RouteNumber
sp <- sp[order(-sp$stops), ]
cat("\nspecies detected:", nrow(sp), " (with /, sp., hybrids:", sum(grepl("/| sp|x ", sp$Species)), ")\n")
cat("species detected at >= 50 stops:", sum(sp$stops >= 50), " >= 200:", sum(sp$stops >= 200), "\n")
print(head(sp, 10)); print(tail(sp, 5))
lon_naamp <- k(sp$Species) %in% lon | sub("^Lithobates_", "Rana_", k(sp$Species)) %in% lon | sub("^Rana_", "Lithobates_", k(sp$Species)) %in% lon
cat("NAAMP species (>=50 stops) with longevity:", sum(lon_naamp & sp$stops >= 50), "/", sum(sp$stops >= 50), "\n")
