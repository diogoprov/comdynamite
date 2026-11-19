# Robustness of species-level dynamics to model specification (v2b vs v4) and dependence on
# detectability. Detectability = empirical per-visit detection rate in stop-years with >= 1
# detection (model-free, unified protocol, all stops).
S <- "data-raw/persistence_proxy"
D <- "data-raw/naamp"
spp <- c("Pseudacris crucifer", "Anaxyrus fowleri", "Lithobates sylvaticus", "Acris crepitans", "Hyla cinerea",
         "Gastrophryne carolinensis", "Lithobates palustris", "Lithobates clamitans", "Lithobates catesbeianus",
         "Anaxyrus americanus")
get <- function(v, sp) {
  f <- file.path(S, "fits", paste0("pilot", v, "_", gsub(" ", "_", sp), ".rds"))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)[[1]]; s <- o$summary
  data.frame(species = sp, version = v, var = s$variable, mean = s$mean, q5 = s$q5, q95 = s$q95, rhat = s$rhat,
             ess = s$ess_bulk, div = sum(o$diagnostics$num_divergent), sw_p = o$ppc["p_value", "switches"])
}
res <- do.call(rbind, c(lapply(spp, get, v = "2b"), lapply(spp, get, v = "4")))
w <- function(vn, v) { x <- subset(res, var == vn & version == v); setNames(x$mean, x$species)[spp] }
runs <- read.csv(file.path(D, "Runs.csv")); names(runs) <- trimws(names(runs)); runs <- subset(runs, UnifiedProtocol == 1)
st <- merge(read.csv(file.path(D, "Stops.csv"))[, c("RunID", "StopNumber", "SiteID", "SkippedStop")],
            runs[, c("RunID", "SurveyYear")]); st <- st[st$SkippedStop == 0, ]
ct <- read.csv(file.path(D, "Counts.csv")); ct <- ct[ct$CallingIndex > 0, ]
detect <- sapply(spp, function(sp) {
  d <- paste(st$RunID, st$StopNumber) %in% with(ct[ct$Species == sp, ], paste(RunID, StopNumber))
  sy <- paste(st$SiteID, st$SurveyYear); occ <- tapply(d, sy, any)
  mean(d[occ[sy]]) })
lon <- readRDS(file.path(S, "species_table.rds")); lon <- setNames(lon$lon, lon$Species)[spp]
tab <- data.frame(species = spp, longevity = lon, detect = round(detect, 3),
  ext_v2 = w("ext_bar", "2b"), ext_v4 = w("ext_bar", "4"), col_v2 = w("col_bar", "2b"), col_v4 = w("col_bar", "4"),
  eps_v2 = w("eps_bar", "2b"), eps_v4 = w("eps_bar", "4"), tau_v2 = w("tau_bar", "2b"), tau_v4 = w("tau_bar", "4"))
cat("Convergence (max R-hat, min bulk ESS, divergences, PPC switches p):\n")
print(do.call(rbind, lapply(split(res, list(res$species, res$version), drop = TRUE), function(x)
  data.frame(species = x$species[1], version = x$version[1], max_rhat = round(max(x$rhat, na.rm = TRUE), 3),
             min_ess = round(min(x$ess, na.rm = TRUE)), div = x$div[1], ppc_sw_p = x$sw_p[1]))), row.names = FALSE)
cat("\nSpecies-level estimates (averaged over random effects):\n"); print(tab, digits = 3, row.names = FALSE)
sc <- function(a, b) { ok <- complete.cases(a, b); round(cor(a[ok], b[ok], method = "spearman"), 2) }
cat("\nSpearman v2b vs v4:  ext", sc(tab$ext_v2, tab$ext_v4), " col", sc(tab$col_v2, tab$col_v4),
    " eps", sc(tab$eps_v2, tab$eps_v4), " tau", sc(tab$tau_v2, tab$tau_v4), "\n")
cat("Spearman with detectability (v4):  ext", sc(tab$ext_v4, tab$detect), " eps", sc(tab$eps_v4, tab$detect),
    " tau", sc(tab$tau_v4, tab$detect), "\n")
cat("Spearman with detectability (v2b): ext", sc(tab$ext_v2, tab$detect), " eps", sc(tab$eps_v2, tab$detect),
    " tau", sc(tab$tau_v2, tab$detect), "\n")
cat("Descriptive only (n = 10), Spearman with longevity (v4): ext", sc(tab$ext_v4, tab$longevity),
    " eps", sc(tab$eps_v4, tab$longevity), " tau", sc(tab$tau_v4, tab$longevity), "\n")
saveRDS(list(res = res, tab = tab), file.path(S, "compare_versions.rds"))
