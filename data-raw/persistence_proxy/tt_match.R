tt <- read.csv("data-raw/tetrapodtraits/TetrapodTraits_v3.0.1.csv")
tt <- subset(tt, Class == "Amphibia" & Order == "Anura")
tt$key <- gsub(" ", "_", tt$Scientific.Name)
pool <- read.csv("data-raw/benicio2021_regional_pool.csv")
pool$key <- NA
for (col in c("name", "portik_name", "atlantic_name")) {
  m <- is.na(pool$key) & pool[[col]] %in% tt$key
  pool$key[m] <- pool[[col]][m]
}
syn <- c(Boana_albomarginata="Hypsiboas_albomarginatus", Boana_albopunctata="Hypsiboas_albopunctatus", Boana_bandeirantes="Hypsiboas_bandeirantes", Boana_bischoffi="Hypsiboas_bischoffi", Boana_caingua="Hypsiboas_caingua", Boana_caipora="Hypsiboas_caipora", Boana_faber="Hypsiboas_faber", Boana_geographica="Hypsiboas_geographicus", Boana_latistriata="Hypsiboas_latistriatus", Boana_lundii="Hypsiboas_lundii", Boana_pardalis="Hypsiboas_pardalis", Boana_polytaenia="Hypsiboas_polytaenius", Boana_prasina="Hypsiboas_prasinus", Boana_pulchella="Hypsiboas_pulchellus", Boana_punctata="Hypsiboas_punctatus", Boana_raniceps="Hypsiboas_raniceps", Boana_semilineata="Hypsiboas_semilineatus", Trachycephalus_typhonius="Trachycephalus_venulosus", Physalaemus_nattereri="Eupemphix_nattereri")
m <- is.na(pool$key) & pool$name %in% names(syn); pool$key[m] <- syn[pool$name[m]]
cat("Anura in TT:", nrow(tt), "\n")
cat("pool:", nrow(pool), " matched:", sum(!is.na(pool$key)),
    " community matched:", sum(!is.na(pool$key[pool$in_community])), "/", sum(pool$in_community), "\n")
cat("unmatched:\n"); print(pool$name[is.na(pool$key)])
d <- tt[match(pool$key[!is.na(pool$key)], tt$key), ]
obs <- function(v, imp = NULL) {
  ok <- !is.na(d[[v]]); if (!is.null(imp)) ok <- ok & d[[imp]] == 0; mean(ok) }
print(round(c(
  body_length_obs = obs("BodyLength_mm", "ImputedLength"),
  body_length_any = obs("BodyLength_mm"),
  mass_obs = obs("BodyMass_g", "ImputedMass"),
  habitat_obs = obs("Ter", "ImputedHabitat"),
  majorhab_obs = obs("MajorHabitatSum", "ImputedMajorHabitat"),
  ecosystem_obs = obs("EcoFresh", "ImputedEcosystem"),
  acttime_obs = obs("Nocturnality", "ImputedActTime"),
  longevity = obs("MaxLongevity"),
  litter = obs("LitterSize"),
  directdev = obs("DirectDev"),
  larval = obs("LarvalStage"),
  diet = obs("DietBreadth"),
  rangesize = obs("RangeSize_km2"),
  poptrend = mean(!is.na(d$PopulationTrend) & d$PopulationTrend != ""),
  tree_imputed = mean(d$TreeImputed == 1)), 2))
cat("\nunique values:\n"); print(table(d$DirectDev, useNA="a")); print(table(d$PopulationTrend, useNA="a"))
print(summary(d$MaxLongevity)); print(summary(d$LitterSize))
print(head(unique(d$SourceLongevity), 5)); print(head(unique(d$SourceLitterSize), 5))
