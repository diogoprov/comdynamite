# Species pools and trees for the case study of the DAMOCLES paper (dev/plano-artigo-damocles.md,
# step E1): the frogs of 14 protected areas of the Atlantic Forest of Sao Paulo of Benicio et al.
# (2021, Dryad doi:10.5061/dryad.sf7m0cg65) on the dated tree of Anura of Portik et al. (2023,
# OSF https://osf.io/3cbsh/; TreePL best tree and 100 dated bootstrap trees).
# 1. A table of names: each taxon of the 248-species tree of Benicio et al. (their pyron
#    version) and of the community matrix, its name in the matrix, a canonical name, its name in
#    Portik et al. and what is done with it (match, rename or graft). Written to
#    data-raw/benicio2021_taxonomy.csv.
# 2. Taxa absent from Portik et al. (undescribed or "aff." taxa, the lineages of Ischnocnema
#    guentheri and a few described species) are grafted in the place they hold in the tree of
#    Benicio et al.: on the stem of the most recent common ancestor (in Portik et al.) of the
#    pool taxa of their sister clade in Benicio et al., at the same relative height between that
#    ancestor and its parent as in Benicio et al.
# 3. The best tree and each bootstrap tree are pruned to the 248-species pool, the missing taxa
#    grafted, and the result pruned to the 112 species of the communities.
# Neither data set is in the repository: the files of Dryad go in DATA_DIR/benicio2021/ and
# those of OSF in DATA_DIR/portik2023/ (DATA_DIR defaults to data-raw). The trees of Portik et al.
# have no open licence, so the pruned trees are written to DATA_DIR/benicio2021_trees.rds and are
# not committed either. phytools (GPL) is only used here.
suppressMessages({
  library(ape)
  library(phytools)
})
data_dir <- Sys.getenv("DATA_DIR", "data-raw")
path <- function(...) file.path(data_dir, ...)
norm <- function(x) sub("_+$", "", gsub("'", "", x))

benicio <- read.tree(path("benicio2021", "tree_pyron_248_species"))
benicio$tip.label <- norm(benicio$tip.label)
comm <- read.csv(path("benicio2021", "frog_data.csv"), row.names = 1, check.names = FALSE)
colnames(comm) <- norm(colnames(comm))
portik <- read.tree(path("portik2023", "TreePL-Rooted_Anura_bestTree.tre"))
portik$tip.label <- norm(portik$tip.label)
boot <- read.tree(path("portik2023", "TreePL-Rooted_Anura_bootstraps.tre"))
stopifnot(length(boot) == 100)

# 1. Names --------------------------------------------------------------------------------------
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
# Canonical names (spelling and current genus) that are in Portik et al.; Hypsiboas is Boana.
renames <- c(
  Aparasphenedon_bokermanni = "Aparasphenodon_bokermanni",
  Bokermanohyla_astartea = "Bokermannohyla_astartea",
  Bokermanohyla_circumdata = "Bokermannohyla_circumdata",
  Bokermanohyla_hylax = "Bokermannohyla_hylax",
  Boana_alfarroi = "Boana_alfaroi", Boana_albomarginatus = "Boana_albomarginata",
  Boana_bischoffii = "Boana_bischoffi", Boana_calcaratus = "Boana_calcarata",
  Boana_fasciatus = "Boana_fasciata", Boana_leptolineatus = "Boana_leptolineata",
  Boana_marginatus = "Boana_marginata", Boana_multifasciatus = "Boana_multifasciata",
  Boana_polytaenius = "Boana_polytaenia", Boana_pulchellus = "Boana_pulchella",
  Boana_semiguttatus = "Boana_semiguttata", Boana_semilineatus = "Boana_semilineata",
  Chiasmocleis_cruci = "Chiasmocleis_crucis", Crossodactylus_schmdti = "Crossodactylus_schmidti",
  Dendrophryniscus_brevipolicatus = "Dendrophryniscus_brevipollicatus",
  Hylodes_hayeri = "Hylodes_heyeri", Ischnocnema_henseli = "Ischnocnema_henselii",
  Lysapus_laevis = "Lysapsus_laevis", Agalychnis_granulosa = "Hylomantis_granulosa",
  Megaelosia_boticariana = "Phantasmarana_boticariana",
  Scinax_argyreonatus = "Scinax_argyreornatus",
  # "Scinax perereca aff. hayii" of the matrix is Scinax_perereca in the feng tree of Benicio et al.
  Scinax_perereca__aff.hayii = "Scinax_perereca"
)
notes <- c(
  Scinax_perereca__aff.hayii = "labelled Scinax_perereca in the feng tree of Benicio et al.",
  Megaelosia_boticariana = "Megaelosia boticariana is now in Phantasmarana",
  Agalychnis_granulosa = "Agalychnis granulosa is now in Hylomantis",
  `Megaelosia_aff._boticariana` = "distinct from M. boticariana in the 248-species tree",
  Ischnocnema_guentheri_C1 = "lineage of the I. guentheri complex; check against current taxonomy",
  Ischnocnema_guentheri_C3 = "lineage of the I. guentheri complex; check against current taxonomy",
  Ischnocnema_guentheri_C4 = "lineage of the I. guentheri complex; check against current taxonomy"
)

taxonomy <- data.frame(tree_name = benicio$tip.label)
taxonomy$comm_name <- ifelse(taxonomy$tree_name %in% names(tree_to_comm),
                             tree_to_comm[taxonomy$tree_name], taxonomy$tree_name)
taxonomy$comm_name[!taxonomy$comm_name %in% colnames(comm)] <- NA
stopifnot(setequal(stats::na.omit(taxonomy$comm_name), colnames(comm)))
taxonomy$name <- ifelse(is.na(taxonomy$comm_name), taxonomy$tree_name, taxonomy$comm_name)
taxonomy$name <- sub("^Hypsiboas_?", "Boana_", taxonomy$name)
portik_name <- ifelse(taxonomy$name %in% names(renames), renames[taxonomy$name], taxonomy$name)
taxonomy$portik_name <- ifelse(portik_name %in% portik$tip.label, portik_name, NA)
taxonomy$action <- ifelse(is.na(taxonomy$portik_name), "graft",
                          ifelse(taxonomy$portik_name == taxonomy$tree_name, "match", "rename"))
key <- ifelse(taxonomy$name %in% names(notes), taxonomy$name, taxonomy$tree_name)
taxonomy$note <- ifelse(key %in% names(notes), notes[key], "")
taxonomy$in_community <- !is.na(taxonomy$comm_name)
stopifnot(!anyDuplicated(taxonomy$name), !anyDuplicated(stats::na.omit(taxonomy$portik_name)))
taxonomy <- taxonomy[order(taxonomy$name), ]
utils::write.csv(taxonomy, "data-raw/benicio2021_taxonomy.csv", row.names = FALSE)
cat("1. taxa:", nrow(taxonomy), "; in the communities:", sum(taxonomy$in_community), "\n")
print(table(action = taxonomy$action, community = taxonomy$in_community))
cat("   grafted:", paste(taxonomy$name[taxonomy$action == "graft"], collapse = ", "), "\n")

# 2. Grafting ----------------------------------------------------------------------------------
# The tree of Benicio et al. with canonical names, the reference for the place of grafted taxa.
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
  if (length(sister) == 1 && sister <= Ntip(reference)) {
    clade_age <- 0
    above <- split_age
  } else {
    clade <- if (length(sister) == 1) sister else split
    clade_age <- ref_height - depth[clade]
    above <- if (clade == split) 0 else split_age - clade_age
  }
  parent <- reference$edge[reference$edge[, 2] == split, 1]
  parent_age <- if (length(parent)) ref_height - depth[parent] else split_age
  list(sister = sister_tips, fraction = if (above == 0) 0 else
    above / (parent_age - clade_age))
}
# Attach `taxon` to `tree` on the stem of the ancestor of the pool taxa of its sister clade.
graft <- function(tree, taxon, place) {
  anchor <- intersect(place$sister, tree$tip.label)
  if (!length(anchor)) return(NULL)
  node <- if (length(anchor) == 1) match(anchor, tree$tip.label) else getMRCA(tree, anchor)
  age <- node_age(tree, node)
  parent <- tree$edge[tree$edge[, 2] == node, 1]
  stem <- tree$edge.length[tree$edge[, 2] == node]
  # Between the ancestor and its parent, clear of both.
  position <- min(max(place$fraction, 0.05), 0.95) * stem
  bind.tip(tree, taxon, edge.length = age + position, where = node, position = position)
}
graft_all <- function(tree, taxa) {
  places <- lapply(stats::setNames(taxa, taxa), graft_place)
  left <- taxa
  while (length(left)) {
    done <- character(0)
    for (taxon in left) {
      grafted <- graft(tree, taxon, places[[taxon]])
      if (!is.null(grafted)) {
        tree <- grafted
        done <- c(done, taxon)
      }
    }
    if (!length(done)) stop("no place for ", paste(left, collapse = ", "), call. = FALSE)
    left <- setdiff(left, done)
  }
  tree
}

# 3. Pools --------------------------------------------------------------------------------------
pools <- list(pool_248 = taxonomy$name, community_112 = taxonomy$name[taxonomy$in_community])
# Grafted on the 248-species pool, which holds the sister clades of all grafted taxa, then
# pruned to the pool.
build <- function(tree, pool) {
  tree$tip.label <- norm(tree$tip.label)
  tree <- keep.tip(tree, stats::na.omit(taxonomy$portik_name))
  tree$tip.label <- taxonomy$name[match(tree$tip.label, taxonomy$portik_name)]
  tree <- keep.tip(graft_all(tree, taxonomy$name[taxonomy$action == "graft"]), pool)
  stopifnot(setequal(tree$tip.label, pool), is.binary(tree), is.ultrametric(tree, tol = 1e-6))
  tree
}
trees <- lapply(pools, function(pool) {
  list(best = build(portik, pool),
       boot = structure(lapply(boot, build, pool = pool), class = "multiPhylo"))
})
for (p in names(trees)) {
  ages <- vapply(trees[[p]]$boot, function(t) max(node.depth.edgelength(t)), numeric(1))
  cat("3.", p, ":", Ntip(trees[[p]]$best), "tips; root age", round(max(node.depth.edgelength(
    trees[[p]]$best)), 1), "Myr (bootstraps", paste(round(range(ages), 1), collapse = " to "),
    ")\n")
}
# Agreement with the tree of Benicio et al.: Robinson-Foulds distance, normalised.
for (p in names(trees)) {
  ref <- keep.tip(reference, pools[[p]])
  cat("   normalised RF distance to the tree of Benicio et al.,", p, ":",
      round(dist.topo(unroot(ref), unroot(trees[[p]]$best)) / (2 * (Ntip(ref) - 3)), 3), "\n")
}
saveRDS(trees, path("benicio2021_trees.rds"))
