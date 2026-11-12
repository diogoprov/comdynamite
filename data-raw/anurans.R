# Anurans of urban ponds in two Brazilian cities (data/anurans.rda), from the raw files of the
# compendium of chapter 1 of the MSc dissertation of L. B. C. Menezes (unpublished). Set
# ANURANS_RAW to the folder `data/raw` of that compendium. The tree is the TreePL-dated tree of
# anurans of Portik et al. (2023, Mol. Biol. Evol. 40: msad109), kept in `data/raw/legacy` of the
# compendium, pruned to the species recorded.
raw <- Sys.getenv("ANURANS_RAW", file.path(
  "~/Library/CloudStorage/Box-Box/Meus Artigos/2026-Leandro Cap 1", "data", "raw"))
abundance <- read.delim(file.path(raw, "abundance.txt"), row.names = 1, check.names = FALSE)
environment <- read.delim(file.path(raw, "environment.txt"))
points <- read.delim(file.path(raw, "ponds_coord.txt"))
rownames(environment) <- environment$ponds
rownames(points) <- points$Ponds
sites <- rownames(abundance)
stopifnot(setequal(sites, rownames(environment)), setequal(sites, rownames(points)))
environment <- environment[sites, ]
points <- points[sites, ]

comm <- as.matrix(abundance)
storage.mode(comm) <- "integer"
comm <- comm[, colSums(comm) > 0]
city <- factor(unname(c(cg = "Campo Grande", sjc = "São José dos Campos")[sub(".*_", "", sites)]))

# hydrology: first principal component of log(x + 1) area, depth, floating vegetation and
# water temperature, standardised across ponds and rescaled to 0-100, as in the compendium;
# the sign is set so that larger, deeper ponds have higher scores
vars <- c("area", "prof", "float_veg", "temp_water")
measured <- as.matrix(environment[, vars])
pc <- stats::prcomp(scale(log1p(pmax(measured, 0))))
score <- pc$x[, 1]
if (cor(score, log1p(environment$area)) < 0) score <- -score
hydrology <- 100 * (score - min(score)) / diff(range(score))
env <- cbind(hydrology = hydrology, area = environment$area, depth = environment$prof,
             floating_vegetation = environment$float_veg,
             water_temperature = environment$temp_water, urban_cover = environment$urb_2022)
rownames(env) <- sites

coords <- as.matrix(points[, c("lon", "lat")])
dimnames(coords) <- list(sites, c("lon", "lat"))

# the tree: L. luctator was split from L. latrans (Magalhães et al. 2020), and the southeastern
# Brazilian populations once called L. latrans are L. luctator, so the tip of L. latrans stands for
# it; Phyllomedusa azurea is Pithecopus azureus after Duellman et al. (2016)
big <- ape::read.tree(file.path(raw, "legacy", "TreePL-Rooted_Anura_bestTree.tre"))
synonyms <- c(Leptodactylus_latrans = "Leptodactylus_luctator",
              Phyllomedusa_azurea = "Pithecopus_azureus")
renamed <- match(big$tip.label, names(synonyms))
big$tip.label[!is.na(renamed)] <- synonyms[renamed[!is.na(renamed)]]
stopifnot(all(colnames(comm) %in% big$tip.label))
tree <- ape::keep.tip(big, colnames(comm))
tree <- ape::ladderize(tree)
stopifnot(ape::is.ultrametric(tree, tol = 1e-6), ape::is.binary(tree))
comm <- comm[, tree$tip.label[order(tree$tip.label)]]

anurans <- list(comm = comm, env = env, city = city, coords = coords, tree = tree)
usethis::use_data(anurans, overwrite = TRUE, compress = "xz")
