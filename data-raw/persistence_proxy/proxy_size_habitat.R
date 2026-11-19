# Body size and habitat (TetrapodTraits v3.0.1, observed values) as candidate proxies of
# persistence: ES-sim against log DR (computed on the full Anura tree of Portik et al. 2023) for
# the anurans of the regional pool of Benicio et al. (2021); and the link of the same traits to
# occupancy (n sites in ATLANTIC AMPHIBIANS, Sao Paulo), the snapshot proxy of epsilon.
suppressMessages(devtools::load_all(quiet = TRUE))
source_es <- parse("data-raw/januario_essim.R")
for (e in source_es) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                         identical(e[[2]], as.name("es_sim"))) eval(e)
RNGkind("L'Ecuyer-CMRG"); set.seed(2028)
portik_dir <- "data-raw/portik2023"
best <- ape::read.tree(file.path(portik_dir, "TreePL-Rooted_Anura_bestTree.tre"))
boot <- ape::read.tree(file.path(portik_dir, "TreePL-Rooted_Anura_bootstraps.tre"))
anura <- function(t) ape::extract.clade(t, ape::getMRCA(t, c("Ascaphus_truei", "Rhinella_marina")))

tt <- read.csv("data-raw/tetrapodtraits/TetrapodTraits_v3.0.1.csv")
tt <- subset(tt, Order == "Anura"); tt$key <- gsub(" ", "_", tt$Scientific.Name)
pool <- read.csv("data-raw/benicio2021_regional_pool.csv")
syn <- c(Boana_albomarginata = "Hypsiboas_albomarginatus", Boana_albopunctata = "Hypsiboas_albopunctatus",
  Boana_bandeirantes = "Hypsiboas_bandeirantes", Boana_bischoffi = "Hypsiboas_bischoffi",
  Boana_caingua = "Hypsiboas_caingua", Boana_caipora = "Hypsiboas_caipora", Boana_faber = "Hypsiboas_faber",
  Boana_geographica = "Hypsiboas_geographicus", Boana_latistriata = "Hypsiboas_latistriatus",
  Boana_lundii = "Hypsiboas_lundii", Boana_pardalis = "Hypsiboas_pardalis",
  Boana_polytaenia = "Hypsiboas_polytaenius", Boana_prasina = "Hypsiboas_prasinus",
  Boana_pulchella = "Hypsiboas_pulchellus", Boana_punctata = "Hypsiboas_punctatus",
  Boana_raniceps = "Hypsiboas_raniceps", Boana_semilineata = "Hypsiboas_semilineatus",
  Trachycephalus_typhonius = "Trachycephalus_venulosus", Physalaemus_nattereri = "Eupemphix_nattereri")
pool$key <- NA
for (col in c("name", "portik_name", "atlantic_name")) {
  m <- is.na(pool$key) & pool[[col]] %in% tt$key; pool$key[m] <- pool[[col]][m] }
m <- is.na(pool$key) & pool$name %in% names(syn); pool$key[m] <- syn[pool$name[m]]
d <- cbind(pool, tt[match(pool$key, tt$key), c("BodyLength_mm", "ImputedLength", "Verticality",
  "Arb", "Ter", "Aqu", "Fos", "ImputedHabitat", "MajorHabitatSum", "ImputedMajorHabitat", "EcoFresh")])
tree_best <- anura(best)
d <- subset(d, portik_name %in% tree_best$tip.label & !is.na(BodyLength_mm) & ImputedLength == 0 &
               !is.na(Verticality) & ImputedHabitat == 0)
d$log_svl <- log(d$BodyLength_mm)
d$log_occ <- ifelse(d$n_sites > 0, log(d$n_sites), NA)
rownames(d) <- d$portik_name
cat("species analysed:", nrow(d), "\n")
traits <- c(log_svl = "log_svl", verticality = "Verticality", habitat_breadth = "MajorHabitatSum")

one_tree <- function(full) {
  t <- anura(full)
  dr <- dr_rates(t)[d$portik_name]
  sub <- ape::keep.tip(t, d$portik_name)
  out <- t(vapply(traits, function(v) {
    x <- stats::setNames(d[[v]], d$portik_name); ok <- !is.na(x)
    es_sim(dr[ok], x[ok], ape::keep.tip(sub, names(x)[ok]), n_null = 999)
  }, numeric(2)))
  data.frame(trait = names(traits), r = out[, "r"], p = out[, "p"])
}
res_best <- one_tree(best)
res_boot <- do.call(rbind, parallel::mclapply(seq_along(boot), function(i)
  cbind(tree = i, one_tree(boot[[i]])), mc.cores = 8))

# Traits against occupancy (snapshot epsilon): PGLS with Pagel's lambda on the best tree.
sub <- ape::keep.tip(tree_best, d$portik_name)
dd <- d[sub$tip.label, ]; dd <- dd[complete.cases(dd[, c("log_occ", "log_svl", "Verticality", "MajorHabitatSum")]), ]
sub <- ape::keep.tip(sub, rownames(dd)); dd <- dd[sub$tip.label, ]
sub_s <- sub; sub_s$edge.length <- sub_s$edge.length / max(ape::node.depth.edgelength(sub_s))
pgls <- nlme::gls(log_occ ~ log_svl + Verticality + MajorHabitatSum, data = dd,
                  correlation = ape::corPagel(0.5, sub_s, form = ~portik_name), method = "ML", control = nlme::glsControl(apVar = FALSE))
dr_best <- dr_rates(tree_best)[rownames(dd)]
pgls_dr_occ <- es_sim(dr_best, stats::setNames(dd$log_occ, rownames(dd)), sub)

cat("\nES-sim, best tree:\n"); print(res_best, digits = 3)
cat("\nES-sim, 100 bootstrap trees (median r, share p <= 0.05):\n")
print(aggregate(cbind(r, sig = p <= 0.05) ~ trait, res_boot, function(v) round(median(v), 3)))
print(aggregate(p ~ trait, res_boot, function(v) round(mean(v <= 0.05), 2)))
cat("\nPGLS log occupancy ~ traits (Pagel lambda):\n"); print(summary(pgls)$tTable, digits = 3)
cat("lambda:", coef(pgls$modelStruct$corStruct, unconstrained = FALSE), "\n")
cat("\nES-sim log DR vs log occupancy:\n"); print(pgls_dr_occ)
saveRDS(list(data = d, es_best = res_best, es_boot = res_boot, pgls = summary(pgls)$tTable,
             es_occ = pgls_dr_occ), "data-raw/persistence_proxy/proxy_size_habitat.rds")
