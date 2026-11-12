# Anurans of urban ponds in two Brazilian cities (data/anurans.rda), from the raw files of the
# compendium of chapter 1 of the MSc dissertation of L. B. C. Menezes (unpublished). Set
# ANURANS_RAW to the folder `data/raw` of that compendium. The tree is the TreePL-dated tree of
# anurans of Portik et al. (2023, Mol. Phylogenet. Evol. 188: 107907), kept in `data/raw/legacy` of the
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

# regional pools: anurans that breed in lentic water and occur both in the ecoregion of each city
# (RESOLVE 2017; the union of Alto Paraná Atlantic forests and Serra do Mar coastal forests for São
# José dos Campos) and within 200 km of its ponds, from GBIF occurrences (accessed 2026-11-11), plus
# the species recorded in its ponds. The lists, queries and scripts are in data-raw/anurans_pool:
# >= 3 records in the ecoregion, >= 1 record within 200 km, and a lentic mode of reproduction in
# Nunes-de-Almeida et al. (2021), in da Silva et al. (2012) when absent from it, or by genus when
# at least three congeners are classified and agree. GBIF records of L. latrans in these regions are
# taken as L. luctator, and the introduced Lithobates catesbeianus is left out.
regional <- read.csv("data-raw/anurans_pool/pools.csv")
regional <- regional[with(regional, (passes_threshold & n_records_200km >= 1 &
                                       lentic_incl_genus_inferred %in% TRUE) | observed), ]
regional$species[regional$species == "Leptodactylus_latrans"] <- "Leptodactylus_luctator"
regional <- regional[regional$species != "Lithobates_catesbeianus", ]
regional <- regional[!duplicated(regional[c("city", "species")]), ]
regional$city <- c(`Campo Grande` = "Campo Grande",
                   `São José dos Campos` = "São José dos Campos")[regional$city]
stopifnot(all(colnames(comm) %in% regional$species))
species <- sort(unique(regional$species))
pool <- vapply(species, function(sp) as.character(city) %in% regional$city[regional$species == sp],
               logical(length(city))) * 1L
rownames(pool) <- sites
stopifnot(all(pool[, colnames(comm)][comm > 0] == 1))

# the tree: the TreePL-dated tree of Portik et al. (2023, Mol. Phylogenet. Evol.), with tips renamed to the names of GBIF
# where they differ (Ololygon, Pithecopus, Phyllomedusa); the tip of L. latrans stands for
# L. luctator, split from it by Magalhães et al. (2020). Species of the pools missing from the tree
# join it at the crown of their genus (or halfway along the branch of the only congener), and the
# polytomies are resolved at random with branches of length zero.
big <- ape::read.tree(file.path(raw, "legacy", "TreePL-Rooted_Anura_bestTree.tre"))
labels <- unique(regional[!is.na(regional$tree_label) & regional$tree_label != "" &
                            regional$tree_label != regional$species, c("tree_label", "species")])
labels$species[labels$tree_label == "Leptodactylus_latrans"] <- "Leptodactylus_luctator"
labels <- unique(labels)
renamed <- match(big$tip.label, labels$tree_label)
big$tip.label[!is.na(renamed)] <- labels$species[renamed[!is.na(renamed)]]
graft <- function(tree, species) {
  genus <- sub("_.*", "", species)
  congeners <- which(sub("_.*", "", tree$tip.label) == genus)
  if (!length(congeners)) stop("No congener of ", species, " in the tree.")
  depth <- ape::node.depth.edgelength(tree)
  height <- max(depth) - depth
  if (length(congeners) == 1) {
    edge <- tree$edge.length[tree$edge[, 2] == congeners]
    return(phytools::bind.tip(tree, species, edge.length = edge / 2, where = congeners,
                              position = edge / 2))
  }
  node <- ape::getMRCA(tree, congeners)
  phytools::bind.tip(tree, species, edge.length = height[node], where = node)
}
missing <- setdiff(species, big$tip.label)
for (sp in missing) big <- graft(big, sp)
set.seed(20261112)
tree <- ape::multi2di(ape::keep.tip(big, species), random = TRUE)
tree <- ape::ladderize(tree)
stopifnot(ape::is.ultrametric(tree, tol = 1e-6), setequal(tree$tip.label, species))
comm <- comm[, sort(colnames(comm))]

anurans <- list(comm = comm, env = env, city = city, coords = coords, pool = pool, tree = tree,
                grafted = missing)
usethis::use_data(anurans, overwrite = TRUE, compress = "xz")
