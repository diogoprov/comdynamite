# NAAMP species-level dynamics (v4 dynamic occupancy, averaged over random effects) against
# longevity, body size and detectability. PGLS (Pagel's lambda, ML) on the Portik et al. (2023)
# best tree; posterior uncertainty propagated by refitting on draws of the response.
suppressMessages({ library(ape); library(nlme) })
F <- "data-raw/persistence_proxy"
D <- "data-raw"
set.seed(2030)
sp_tab <- readRDS(file.path(F, "species_table.rds")); sp_tab <- subset(sp_tab, !is.na(lon))

# ---- species-level estimates and diagnostics ----
one <- function(sp) {
  f <- file.path(F, "fits", paste0("pilot4_", gsub(" ", "_", sp), ".rds"))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)[[1]]; s <- as.data.frame(o$summary)
  core <- s[!grepl("_m$|_bar$", s$variable), ]
  g <- o$gq
  summ <- function(v) {
    if (!is.null(g)) { x <- g[[v]]; c(median(x, na.rm = TRUE), quantile(x, c(.05, .95), na.rm = TRUE), mean(is.na(x))) }
    else { r <- s[s$variable == v, ]; c(r$mean, r$q5, r$q95, 0) } }
  est <- sapply(c("col_bar", "ext_bar", "eps_bar", "tau_bar"), summ)
  data.frame(species = sp, routes = if (!is.null(o$data)) o$data$R else NA,
    max_rhat = max(core$rhat, na.rm = TRUE), min_ess = min(core$ess_bulk, na.rm = TRUE),
    div = sum(o$diagnostics$num_divergent), ppc_switch_p = o$ppc["p_value", "switches"],
    col = est[1, 1], ext = est[1, 2], ext_q5 = est[2, 2], ext_q95 = est[3, 2],
    eps = est[1, 3], tau = est[1, 4], nan_share = max(est[4, ]))
}
tab <- do.call(rbind, lapply(sp_tab$Species, one))
logs <- list.files(file.path(F, "fits", "logs"), "^(v4_|batch_v4_|pilot4_)", full.names = TRUE)
routes <- sapply(tab$species, function(sp) {
  l <- unlist(lapply(logs, readLines)); l <- grep(paste0("^", sp, " : stops"), l, value = TRUE)
  if (length(l)) as.integer(sub(".*routes ([0-9]+).*", "\\1", l[1])) else NA })
tab$routes <- routes

# ---- covariates ----
k <- function(x) gsub(" ", "_", x)
tab$longevity <- sp_tab$lon[match(tab$species, sp_tab$Species)]
runs <- read.csv(file.path(D, "naamp/Runs.csv")); names(runs) <- trimws(names(runs)); runs <- subset(runs, UnifiedProtocol == 1)
st <- merge(read.csv(file.path(D, "naamp/Stops.csv"))[, c("RunID", "StopNumber", "SiteID", "SkippedStop")],
            runs[, c("RunID", "SurveyYear")]); st <- st[st$SkippedStop == 0, ]
ct <- read.csv(file.path(D, "naamp/Counts.csv")); ct <- ct[ct$CallingIndex > 0, ]
stkey <- paste(st$RunID, st$StopNumber); sy <- paste(st$SiteID, st$SurveyYear)
tab$detect <- sapply(tab$species, function(sp) {
  d <- stkey %in% with(ct[ct$Species == sp, ], paste(RunID, StopNumber)); occ <- tapply(d, sy, any); mean(d[occ[sy]]) })
tt <- read.csv(file.path(D, "tetrapodtraits/TetrapodTraits_v3.0.1.csv")); tt <- subset(tt, Order == "Anura")
tt$key <- k(tt$Scientific.Name)
rana <- function(x) { y <- sub("^Lithobates_", "Rana_", x)
  sub("^(Rana_)(sylvatic|areolat|catesbeian)us$", "\\1\\2a", y) }
ttk <- ifelse(k(tab$species) %in% tt$key, k(tab$species), rana(k(tab$species)))
i <- match(ttk, tt$key)
tab$svl <- tt$BodyLength_mm[i]; tab$svl_imputed <- tt$ImputedLength[i]

# ---- phylogeny ----
tree <- read.tree(file.path(D, "portik2023/TreePL-Rooted_Anura_bestTree.tre"))
tab$tip <- ifelse(k(tab$species) %in% tree$tip.label, k(tab$species), rana(k(tab$species)))
cat("not in tree:", tab$species[!tab$tip %in% tree$tip.label], "\n")

# ---- inclusion ----
tab$ok <- with(tab, routes >= 10 & max_rhat < 1.05 & nan_share < 0.05 & tip %in% tree$tip.label)
cat("\nAll fitted species:\n")
print(transform(tab[order(tab$ext), c("species", "routes", "max_rhat", "min_ess", "div", "ppc_switch_p", "col", "ext",
  "ext_q5", "ext_q95", "eps", "tau", "nan_share", "longevity", "svl", "detect", "ok")]), digits = 3, row.names = FALSE)
d <- subset(tab, ok); rownames(d) <- d$tip
d$log_e <- log(d$ext); d$log_lon <- log(d$longevity); d$log_svl <- log(d$svl); d$logit_det <- qlogis(d$detect)
sub <- keep.tip(tree, d$tip); sub$edge.length <- sub$edge.length / max(node.depth.edgelength(sub)); d <- d[sub$tip.label, ]
cat("\nn analysed:", nrow(d), "\n")
cat("Pearson among predictors:\n"); print(round(cor(d[, c("log_lon", "log_svl", "logit_det")]), 2))

fit <- function(f, dd = d, tr = sub) {
  # Pagel's lambda by profile likelihood on [0, 1] (unbounded corPagel optimisation is unstable at n ~ 20)
  grid <- seq(0, 1, by = 0.05)
  ms <- lapply(grid, function(l) tryCatch(gls(f, data = dd, correlation = corPagel(l, tr, form = ~tip, fixed = TRUE),
                                              method = "ML"), error = function(e) NULL))
  ll <- sapply(ms, function(m) if (is.null(m)) -Inf else as.numeric(logLik(m)))
  ms[[which.max(ll)]] }
lam <- function(m) round(coef(m$modelStruct$corStruct, unconstrained = FALSE), 2)
models <- list(
  "log e ~ log longevity" = log_e ~ log_lon,
  "log e ~ log longevity + logit detect" = log_e ~ log_lon + logit_det,
  "log e ~ log longevity + log SVL + logit detect" = log_e ~ log_lon + log_svl + logit_det,
  "eps ~ log longevity + logit detect" = eps ~ log_lon + logit_det,
  "tau ~ log longevity + logit detect" = tau ~ log_lon + logit_det,
  "log e ~ log SVL + logit detect" = log_e ~ log_svl + logit_det)
out <- list()
for (nm in names(models)) {
  m <- fit(models[[nm]]); tt_ <- summary(m)$tTable
  cat("\n==", nm, " (lambda =", lam(m), ")\n"); print(round(tt_, 3))
  out[[nm]] <- tt_ }

# ---- propagate posterior uncertainty in e (fits with saved draws; others normal on log scale from q5-q95) ----
draw_e <- function(sp) {
  o <- readRDS(file.path(F, "fits", paste0("pilot4_", gsub(" ", "_", sp), ".rds")))[[1]]
  if (!is.null(o$gq)) { x <- o$gq$ext_bar; x <- x[!is.na(x)]; sample(x, 200, replace = TRUE) }
  else { r <- d[d$species == sp, ]; exp(rnorm(200, log(r$ext), (log(r$ext_q95) - log(r$ext_q5)) / (2 * 1.645))) } }
E <- sapply(d$species, draw_e)
slopes <- t(sapply(seq_len(200), function(j) {
  dd <- d; dd$log_e <- log(E[j, ]); m <- fit(log_e ~ log_lon + logit_det, dd); coef(m)[c("log_lon", "logit_det")] }))
cat("\nPosterior propagation, log e ~ log longevity + logit detect (200 draws): median and 90% interval of slopes\n")
print(round(apply(slopes, 2, quantile, c(.05, .5, .95)), 3))
cat("share of draws with longevity slope < 0:", mean(slopes[, "log_lon"] < 0), "\n")
d2 <- subset(tab, routes >= 10 & max_rhat < 1.10 & nan_share < 0.10); rownames(d2) <- d2$tip
d2$log_e <- log(d2$ext); d2$log_lon <- log(d2$longevity); d2$logit_det <- qlogis(d2$detect)
sub2 <- keep.tip(tree, d2$tip); sub2$edge.length <- sub2$edge.length / max(node.depth.edgelength(sub2)); d2 <- d2[sub2$tip.label, ]
m2 <- fit(log_e ~ log_lon + logit_det, d2, sub2)
cat("\nSensitivity (R-hat < 1.10, NaN share < 0.10; n =", nrow(d2), ", lambda =", lam(m2), "):\n"); print(round(summary(m2)$tTable, 3))
saveRDS(list(table = tab, data = d, models = out, slopes = slopes), file.path(F, "pgls_persistence.rds"))
write.csv(tab, file.path(F, "species_dynamics_v4.csv"), row.names = FALSE)

# ---- abundance hypothesis: mean calling index (1-3) among detections as a local abundance proxy ----
ct_all <- read.csv(file.path(D, "naamp/Counts.csv")); ct_all <- ct_all[ct_all$CallingIndex > 0, ]
ct_all <- ct_all[paste(ct_all$RunID, ct_all$StopNumber) %in% stkey, ]
ci <- tapply(ct_all$CallingIndex, ct_all$Species, mean)
d$calling <- ci[d$species]; d2$calling <- ci[d2$species]
cat("\nCorrelation mean calling index vs logit detect:", round(cor(d$calling, d$logit_det), 2), "\n")
for (f in list(log_e ~ calling, log_e ~ calling + logit_det, log_e ~ calling + logit_det + log_lon)) {
  m <- fit(f); cat("\n==", deparse(f), " (n =", nrow(d), ", lambda =", lam(m), ")\n"); print(round(summary(m)$tTable, 3)) }
m <- fit(log_e ~ calling + logit_det, d2, sub2)
cat("\n== relaxed filter, log_e ~ calling + logit_det (n =", nrow(d2), ")\n"); print(round(summary(m)$tTable, 3))
tab$calling <- ci[tab$species]; write.csv(tab, file.path(F, "species_dynamics_v4.csv"), row.names = FALSE)
