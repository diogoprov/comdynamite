# Species pools and trees for the case study of the DAMOCLES paper (dev/plano-artigo-damocles.md,
# step E1): the frogs of 14 protected areas of the Atlantic Forest of Sao Paulo of Benicio et al.
# (2021, Dryad doi:10.5061/dryad.sf7m0cg65) and the regional pool of the Serra do Mar and the
# interior (with its Cerrado) from ATLANTIC AMPHIBIANS (Vancine et al. 2018, Ecology 99:1692,
# doi:10.1002/ecy.2392), on the dated tree of Anura of Portik et al. (2023, OSF
# https://osf.io/3cbsh/; TreePL best tree and 100 dated bootstrap trees).
# 1. A table of names for the 248-species tree of Benicio et al. (their pyron version) and the
#    community matrix: name in that tree, name in the matrix, canonical name, name in Portik et
#    al., and what is done (match, rename or graft). The canonical name is the one in Portik et
#    al. when there is one, else a corrected name; the three unnamed lineages of Ischnocnema
#    guentheri (C1, C3, C4 of Benicio et al.) are Ischnocnema_sp1 to sp3. Written to
#    data-raw/benicio2021_taxonomy.csv; join the matrix to the trees by comm_name -> name.
# 2. The regional pool: anurans with a valid name in the sites of ATLANTIC AMPHIBIANS in the state
#    of Sao Paulo, which span the coast, the Serra do Mar and the interior with its Cerrado (the
#    14 areas lie within it), minus the introduced Rana catesbeiana, plus the 112 species of the
#    communities; I. guentheri is left out, its records taken as the lineages Ischnocnema_sp1 to
#    sp3. Names harmonised with Portik et al. (Ololygon in Scinax, Pithecopus in
#    Phyllomedusa, Megaelosia boticariana and M. massarti in Phantasmarana, gender of the
#    epithet). Written to data-raw/benicio2021_regional_pool.csv.
# 3. One master tree per tree of Portik et al., from which every pool is pruned: the tips renamed,
#    then
#    - taxa of the tree of Benicio et al. absent from Portik et al. grafted in the place they
#      hold in the tree of Benicio et al.: on the stem of the most recent common ancestor of their
#      sister clade, at the same relative height between that ancestor and its parent;
#    - species of the regional pool absent from both grafted at a point drawn uniformly along
#      the branches of the clade of their congeners in ATLANTIC AMPHIBIANS (for Ololygon, the
#      species ATLANTIC AMPHIBIANS places in it), a new draw in each tree.
# 4. The best tree and the 100 bootstrap trees pruned to three pools: the 112 species of the
#    communities, the 248 of the tree of Benicio et al. (species with sequences, not a regional
#    pool) and the regional pool.
# The data sets are not in the repository: the files of Dryad go in DATA_DIR/benicio2021/, those
# of OSF in DATA_DIR/portik2023/ (DATA_DIR defaults to data-raw) and ATLANTIC AMPHIBIANS in
# ATLANTIC_DIR. The trees of Portik et al. have no open licence, so the pruned trees are written to
# DATA_DIR/benicio2021_trees.rds and not committed either. phytools (GPL) is only used here.
suppressMessages({
  library(ape)
  library(phytools)
})
RNGkind("L'Ecuyer-CMRG")
set.seed(2021)
data_dir <- Sys.getenv("DATA_DIR", "data-raw")
atlantic_dir <- Sys.getenv("ATLANTIC_DIR", file.path(data_dir, "atlantic_amphibians"))
path <- function(...) file.path(data_dir, ...)
norm <- function(x) sub("_+$", "", gsub("'", "", x))

benicio <- read.tree(path("benicio2021", "tree_pyron_248_species"))
benicio$tip.label <- norm(benicio$tip.label)
comm <- read.csv(path("benicio2021", "frog_data.csv"), row.names = 1, check.names = FALSE)
colnames(comm) <- norm(colnames(comm))
portik <- read.tree(path("portik2023", "TreePL-Rooted_Anura_bestTree.tre"))
boot <- read.tree(path("portik2023", "TreePL-Rooted_Anura_bootstraps.tre"))
stopifnot(length(boot) == 100)
portik_tips <- norm(portik$tip.label)

# 1. Names of Benicio et al. --------------------------------------------------------------------
# Labels of the 248-species tree written differently in the matrix (same taxa).
tree_to_comm <- c(
  Adenomerasp_J = "Adenomera_sp_J", Aplastodiscussp_5 = "Aplastodiscus_sp_5",
  Aparasphenedon_bokermanni = "Aparasphenodon_bokermanni",
  Bokermannohyla_astartea = "Bokermanohyla_astartea",
  Bokermannohyla_circumdata = "Bokermanohyla_circumdata",
  Bokermannohyla_hylax = "Bokermanohyla_hylax",
  Dendrophryniscus_brevipollicatus = "Dendrophryniscus_brevipolicatus",
  Dendrophryniscus_sp = "Dendrophryniscus_sp.", Hypsiboas_bischoffi = "Hypsiboas_bischoffii",
  Hypsiboascaipora = "Hypsiboas_caipora", Ischnocnemaguentheri_C1 = "Ischnocnema_guentheri_C1",
  Ischnocnemaguentheri_C3 = "Ischnocnema_guentheri_C3",
  Ischnocnemaguentheri_C4 = "Ischnocnema_guentheri_C4",
  `Scinaxperereca-aff.hayii` = "Scinax_perereca__aff.hayii"
)
stopifnot(all(names(tree_to_comm) %in% benicio$tip.label),
          all(tree_to_comm %in% colnames(comm)))
# Canonical names (spelling, current genus, gender of the epithet); Hypsiboas is Boana.
renames <- c(
  Adenomera_sp_J = "Adenomera_sp_J", Aplastodiscus_sp_5 = "Aplastodiscus_sp_5",
  Aparasphenedon_bokermanni = "Aparasphenodon_bokermanni",
  Bokermanohyla_astartea = "Bokermannohyla_astartea",
  Bokermanohyla_circumdata = "Bokermannohyla_circumdata",
  Bokermanohyla_hylax = "Bokermannohyla_hylax",
  Bokermannohyla_izecksoni = "Bokermannohyla_izecksohni",
  Boana_alfarroi = "Boana_alfaroi", Boana_albomarginatus = "Boana_albomarginata",
  Boana_bischoffii = "Boana_bischoffi", Boana_calcaratus = "Boana_calcarata",
  Boana_fasciatus = "Boana_fasciata", Boana_leptolineatus = "Boana_leptolineata",
  Boana_marginatus = "Boana_marginata", Boana_multifasciatus = "Boana_multifasciata",
  Boana_polytaenius = "Boana_polytaenia", Boana_pulchellus = "Boana_pulchella",
  Boana_semiguttatus = "Boana_semiguttata", Boana_semilineatus = "Boana_semilineata",
  Chiasmocleis_cruci = "Chiasmocleis_crucis", Crossodactylus_schmdti = "Crossodactylus_schmidti",
  Dendrophryniscus_brevipolicatus = "Dendrophryniscus_brevipollicatus",
  Dendrophryniscus_sp. = "Dendrophryniscus_sp", Fritziana_aff._fissilis = "Fritziana_aff_fissilis",
  Hylodes_hayeri = "Hylodes_heyeri", Ischnocnema_henseli = "Ischnocnema_henselii",
  # Unnamed lineages of the I. guentheri complex, with deep divergences.
  Ischnocnema_guentheri_C1 = "Ischnocnema_sp1", Ischnocnema_guentheri_C3 = "Ischnocnema_sp2",
  Ischnocnema_guentheri_C4 = "Ischnocnema_sp3",
  Lysapus_laevis = "Lysapsus_laevis", Agalychnis_granulosa = "Hylomantis_granulosa",
  Megaelosia_boticariana = "Phantasmarana_boticariana",
  `Megaelosia_aff._boticariana` = "Phantasmarana_aff_boticariana",
  `Paratelmatobius_aff._cardosoi` = "Paratelmatobius_aff_cardosoi",
  Scinax_argyreonatus = "Scinax_argyreornatus",
  # "Scinax perereca aff. hayii" of the matrix is Scinax_perereca in the feng tree of Benicio et al.
  Scinax_perereca__aff.hayii = "Scinax_perereca"
)
notes <- c(
  Scinax_perereca = "labelled Scinax_perereca in the feng tree of Benicio et al.",
  Phantasmarana_boticariana = "Megaelosia boticariana is now in Phantasmarana",
  Hylomantis_granulosa = "Agalychnis granulosa is now in Hylomantis",
  Phantasmarana_aff_boticariana = "distinct from M. boticariana in the 248-species tree",
  Ischnocnema_sp1 = "Ischnocnema guentheri C1 of Benicio et al., unnamed lineage",
  Ischnocnema_sp2 = "Ischnocnema guentheri C3 of Benicio et al., unnamed lineage",
  Ischnocnema_sp3 = "Ischnocnema guentheri C4 of Benicio et al., unnamed lineage"
)

taxonomy <- data.frame(tree_name = benicio$tip.label)
taxonomy$comm_name <- ifelse(taxonomy$tree_name %in% names(tree_to_comm),
                             tree_to_comm[taxonomy$tree_name], taxonomy$tree_name)
taxonomy$comm_name[!taxonomy$comm_name %in% colnames(comm)] <- NA
stopifnot(setequal(stats::na.omit(taxonomy$comm_name), colnames(comm)))
original <- ifelse(is.na(taxonomy$comm_name), taxonomy$tree_name, taxonomy$comm_name)
original <- sub("^Hypsiboas_?", "Boana_", original)
taxonomy$name <- ifelse(original %in% names(renames), renames[original], original)
taxonomy$portik_name <- ifelse(taxonomy$name %in% portik_tips, taxonomy$name, NA)
taxonomy$action <- ifelse(is.na(taxonomy$portik_name), "graft",
                          ifelse(taxonomy$portik_name == taxonomy$tree_name, "match", "rename"))
taxonomy$note <- ifelse(taxonomy$name %in% names(notes), notes[taxonomy$name], "")
taxonomy$in_community <- !is.na(taxonomy$comm_name)
stopifnot(!anyDuplicated(taxonomy$name))
taxonomy <- taxonomy[order(taxonomy$name), ]
utils::write.csv(taxonomy, "data-raw/benicio2021_taxonomy.csv", row.names = FALSE)
cat("1. taxa:", nrow(taxonomy), "; in the communities:", sum(taxonomy$in_community), "\n")
print(table(action = taxonomy$action, community = taxonomy$in_community))
cat("   grafted:", paste(taxonomy$name[taxonomy$action == "graft"], collapse = ", "), "\n")

# 2. Regional pool -----------------------------------------------------------------------------
sites <- read.csv(file.path(atlantic_dir, "ATLANTIC_AMPHIBIANS_sites.csv"))
species <- read.csv(file.path(atlantic_dir, "ATLANTIC_AMPHIBIANS_species.csv"))
species <- species[species$order %in% "Anura" & !is.na(species$valid_name), ]
species$valid_name <- gsub(" ", "_", species$valid_name)
sp_sites <- sites$id[sites$state_abbreviation == "BR-SP"]
introduced <- "Lithobates_catesbeianus"
# Records of I. guentheri in Sao Paulo are taken as the unnamed lineages of the complex that
# Benicio et al. separate (Ischnocnema_sp1 to sp3), which enter the pool from the communities;
# keeping both would count the complex twice.
merged <- "Ischnocnema_guentheri"
genders <- function(e) {
  unique(c(e, sub("a$", "us", e), sub("us$", "a", e), sub("um$", "us", e), sub("ae$", "i", e),
           sub("i$", "ae", e), sub("is$", "e", e), sub("e$", "is", e)))
}
genera <- list(Ololygon = "Scinax", Pithecopus = "Phyllomedusa", Hypsiboas = "Boana",
               Megaelosia = "Phantasmarana")
harmonise <- function(name) {
  genus <- sub("_.*", "", name)
  epithet <- sub("^[^_]+_", "", name)
  hits <- intersect(as.vector(outer(c(genus, genera[[genus]]), genders(epithet), paste,
                                    sep = "_")), portik_tips)
  if (length(hits)) hits[1] else NA_character_
}
regional <- data.frame(atlantic_name = sort(unique(species$valid_name[species$id %in% sp_sites])))
regional <- regional[!regional$atlantic_name %in% c(introduced, merged), , drop = FALSE]
regional$n_sites <- as.vector(table(species$valid_name[species$id %in% sp_sites])[
  regional$atlantic_name])
regional$portik_name <- vapply(regional$atlantic_name, harmonise, "")
# Absent from Portik et al. but grafted by the tree of Benicio et al.
regional$name <- ifelse(!is.na(regional$portik_name), regional$portik_name,
  ifelse(regional$atlantic_name %in% taxonomy$name, regional$atlantic_name,
         sub("^Ololygon_", "Scinax_", regional$atlantic_name)))
regional$name[regional$name == "Scinax_atrata"] <- "Scinax_atratus"
regional$action <- ifelse(!is.na(regional$portik_name), "match",
  ifelse(regional$name %in% taxonomy$name, "graft_benicio", "graft_genus"))
# Clade for the grafts by genus: the congeners in all of ATLANTIC AMPHIBIANS.
atlantic_portik <- vapply(unique(species$valid_name), harmonise, "")
regional$graft_clade <- ""
for (i in which(regional$action == "graft_genus")) {
  genus <- sub("_.*", "", regional$atlantic_name[i])
  congeners <- atlantic_portik[startsWith(names(atlantic_portik), paste0(genus, "_"))]
  congeners <- sort(unique(stats::na.omit(congeners)))
  if (length(congeners) < 2) congeners <- portik_tips[startsWith(portik_tips, paste0(genus, "_"))]
  regional$graft_clade[i] <- paste(congeners, collapse = ";")
}
community <- taxonomy$name[taxonomy$in_community]
regional$in_community <- regional$name %in% community
extra <- setdiff(community, regional$name)
regional <- rbind(regional, data.frame(
  atlantic_name = NA, n_sites = 0, portik_name = taxonomy$portik_name[match(extra, taxonomy$name)],
  name = extra, action = ifelse(extra %in% portik_tips, "match", "graft_benicio"),
  graft_clade = "", in_community = TRUE))
stopifnot(!anyDuplicated(regional$name))
regional <- regional[order(regional$name), ]
utils::write.csv(regional, "data-raw/benicio2021_regional_pool.csv", row.names = FALSE)
cat("2. regional pool:", nrow(regional), "species (", sum(!is.na(regional$atlantic_name)),
    "from", length(sp_sites), "sites in Sao Paulo,", length(extra),
    "only in the communities )\n")
print(table(regional$action))
cat("   grafted by genus:", paste(regional$name[regional$action == "graft_genus"],
                                  collapse = ", "), "\n")

# 3. Master trees ------------------------------------------------------------------------------
# The tree of Benicio et al. with canonical names, the reference for the place of its taxa.
reference <- benicio
reference$tip.label <- taxonomy$name[match(reference$tip.label, taxonomy$tree_name)]
ref_height <- max(node.depth.edgelength(reference))
node_age <- function(tree, node) {
  depth <- node.depth.edgelength(tree)
  max(depth[seq_len(Ntip(tree))]) - depth[node]
}
# Place of `taxon` in the reference tree: its sister clade, and the relative height of its split
# between the ancestor of the sister clade and the parent of that ancestor.
graft_place <- function(taxon) {
  tip <- match(taxon, reference$tip.label)
  split <- reference$edge[reference$edge[, 2] == tip, 1]
  sister <- setdiff(reference$edge[reference$edge[, 1] == split, 2], tip)
  sister_tips <- unlist(lapply(sister, function(s) {
    if (s <= Ntip(reference)) reference$tip.label[s] else extract.clade(reference, s)$tip.label
  }))
  depth <- node.depth.edgelength(reference)
  split_age <- ref_height - depth[split]
  clade_age <- if (length(sister) == 1 && sister <= Ntip(reference)) 0 else
    ref_height - depth[sister[1]]
  parent <- reference$edge[reference$edge[, 2] == split, 1]
  parent_age <- if (length(parent)) ref_height - depth[parent] else split_age
  list(sister = sister_tips, fraction = (split_age - clade_age) / (parent_age - clade_age))
}
# Attach `taxon` on the stem of the ancestor of the tips of its sister clade present in `tree`.
graft_sister <- function(tree, taxon, place) {
  anchor <- intersect(place$sister, tree$tip.label)
  if (!length(anchor)) return(NULL)
  node <- if (length(anchor) == 1) match(anchor, tree$tip.label) else getMRCA(tree, anchor)
  stem <- tree$edge.length[tree$edge[, 2] == node]
  # Between the ancestor and its parent, clear of both.
  position <- min(max(place$fraction, 0.05), 0.95) * stem
  bind.tip(tree, taxon, edge.length = node_age(tree, node) + position, where = node,
           position = position)
}
# Attach `taxon` at a point drawn uniformly along the branches of the clade of `congeners`.
graft_uniform <- function(tree, taxon, congeners) {
  congeners <- intersect(congeners, tree$tip.label)
  crown <- if (length(congeners) == 1) match(congeners, tree$tip.label) else
    getMRCA(tree, congeners)
  edges <- which(tree$edge[, 2] %in% c(crown, getDescendants(tree, crown)))
  if (length(congeners) > 1) edges <- setdiff(edges, which(tree$edge[, 2] == crown))
  e <- edges[sample.int(length(edges), 1, prob = tree$edge.length[edges])]
  node <- tree$edge[e, 2]
  position <- stats::runif(1, 0, tree$edge.length[e])
  bind.tip(tree, taxon, edge.length = node_age(tree, node) + position, where = node,
           position = position)
}
benicio_grafts <- taxonomy$name[taxonomy$action == "graft"]
places <- lapply(stats::setNames(benicio_grafts, benicio_grafts), graft_place)
genus_grafts <- regional[regional$action == "graft_genus", ]
master <- function(tree) {
  tree$tip.label <- norm(tree$tip.label)
  tree <- drop.tip(tree, setdiff(tree$tip.label, c(taxonomy$portik_name, regional$portik_name,
    unlist(strsplit(genus_grafts$graft_clade, ";")))))
  left <- benicio_grafts
  while (length(left)) {
    done <- character(0)
    for (taxon in left) {
      grafted <- graft_sister(tree, taxon, places[[taxon]])
      if (!is.null(grafted)) {
        tree <- grafted
        done <- c(done, taxon)
      }
    }
    if (!length(done)) stop("no place for ", paste(left, collapse = ", "), call. = FALSE)
    left <- setdiff(left, done)
  }
  for (i in seq_len(nrow(genus_grafts))) {
    tree <- graft_uniform(tree, genus_grafts$name[i],
                          strsplit(genus_grafts$graft_clade[i], ";")[[1]])
  }
  tree
}

# 4. Pools --------------------------------------------------------------------------------------
pools <- list(community_112 = community, pool_248 = taxonomy$name, regional = regional$name)
prune <- function(tree, pool) {
  tree <- keep.tip(tree, pool)
  stopifnot(setequal(tree$tip.label, pool), is.binary(tree), is.ultrametric(tree, tol = 1e-6))
  tree
}
masters <- c(list(master(portik)), lapply(boot, master))
trees <- lapply(pools, function(pool) {
  list(best = prune(masters[[1]], pool),
       boot = structure(lapply(masters[-1], prune, pool = pool), class = "multiPhylo"))
})
for (p in names(trees)) {
  ages <- vapply(trees[[p]]$boot, function(t) max(node.depth.edgelength(t)), numeric(1))
  cat("4.", p, ":", Ntip(trees[[p]]$best), "tips; root age", round(max(node.depth.edgelength(
    trees[[p]]$best)), 1), "Myr (bootstraps", paste(round(range(ages), 1), collapse = " to "),
    ")\n")
}
# Agreement with the tree of Benicio et al.: Robinson-Foulds distance, normalised.
for (p in c("community_112", "pool_248")) {
  ref <- keep.tip(reference, pools[[p]])
  cat("   normalised RF distance to the tree of Benicio et al.,", p, ":",
      round(dist.topo(unroot(ref), unroot(trees[[p]]$best)) / (2 * (Ntip(ref) - 3)), 3), "\n")
}
saveRDS(trees, path("benicio2021_trees.rds"))
