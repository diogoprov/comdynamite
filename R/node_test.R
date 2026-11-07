#' Node contrasts of a phylogeny
#'
#' Codes each node of a tree as a contrast between its two descendant clades, as
#' in Leibold et al. (2010). One clade starts at -1 and the other at 1, and the
#' value is halved at each node on the way to the tips. Species outside the node
#' get zero.
#'
#' With `branch_lengths = TRUE`, values follow section III of the supplement of Leibold
#' et al. (2010): each branch between a species and the node counts its length divided by
#' the number of species below it, and the sum is divided by the total length of that side
#' of the node. Species that share history get less weight, and each side still sums
#' to 1. The supplement proposes this coding for trees in which all species are equally
#' distant from the root.
#'
#' @param tree A `phylo` object with dichotomous nodes.
#' @param cherries Keep nodes that join two tips? The article tests them (its Figure 3),
#'   but the function `CreateNodeContrasts` of the authors leaves them out.
#' @param branch_lengths Use the branch lengths of `tree`?
#' @return Matrix with one row per tip and one column per node, named by the node
#'   numbers of `tree`.
#' @note Adapted from the MATLAB function `CreateNodeContrasts` in the supplement of
#'   Leibold et al. (2010).
#' @references
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @examples
#' node_contrasts(ape::read.tree(text = "((a,b),(c,(d,e)));"))
#'
#' # Figure S2 of the supplement
#' tree <- ape::read.tree(text = "((((A:0.2,B:0.2):0.65,C:0.85):0.15,D:1):0.4,(E:0.5,F:0.5):0.9);")
#' round(node_contrasts(tree, branch_lengths = TRUE), 3)
#' @export
node_contrasts <- function(tree, cherries = TRUE, branch_lengths = FALSE) {
  stopifnot(inherits(tree, "phylo"))
  if (branch_lengths && is.null(tree$edge.length)) {
    stop("`tree` has no branch lengths", call. = FALSE)
  }
  n_tips <- length(tree$tip.label)
  children <- split(tree$edge[, 2], tree$edge[, 1])
  if (any(lengths(children) != 2)) stop("the tree must be dichotomous and rooted", call. = FALSE)
  # value of each tip below a node, starting from the value of the node
  spread <- function(node, value) {
    if (node <= n_tips) return(stats::setNames(value, node))
    unlist(lapply(children[[as.character(node)]], spread, value = value / 2),
           use.names = TRUE)
  }
  edge_length <- numeric(max(tree$edge))
  if (branch_lengths) edge_length[tree$edge[, 2]] <- tree$edge.length
  n_below <- function(node) {
    if (node <= n_tips) return(1)
    sum(vapply(children[[as.character(node)]], n_below, numeric(1)))
  }
  # branch lengths on the way to each tip, shared among the species below each branch
  shared <- function(node, value) {
    value <- value + edge_length[node] / n_below(node)
    if (node <= n_tips) return(stats::setNames(value, node))
    unlist(lapply(children[[as.character(node)]], shared, value = value), use.names = TRUE)
  }
  side_values <- function(node, sign) {
    if (!branch_lengths) return(spread(node, sign))
    values <- shared(node, 0)
    # a clade with branches of length zero has equal weights
    if (sum(values) == 0) values[] <- 1
    sign * values / sum(values)
  }
  informative <- Filter(function(node) cherries || any(children[[node]] > n_tips), names(children))
  contrasts <- vapply(informative, function(node) {
    pair <- children[[node]]
    out <- numeric(n_tips)
    for (side in 1:2) {
      values <- side_values(pair[side], c(-1, 1)[side])
      out[as.integer(sub(".*\\.", "", names(values)))] <- values
    }
    out
  }, numeric(n_tips))
  matrix(contrasts, n_tips, dimnames = list(tree$tip.label, informative))
}

#' Test clades against the environment, node by node
#'
#' Implements the regression of Leibold et al. (2010). For each node, sites without
#' species of that node are removed, the contrast is standardised with species
#' weights and summed within sites, and the result is regressed on the predictors
#' standardised with site weights, using site richness as weights. Each node shows
#' whether the two clades that split there differ in their response to the predictors.
#'
#' The test statistic is the explained over the residual sum of squares. It is
#' tested by permuting sites, and by replacing the species of the node with a random
#' set of the same size from all species. Following Peres-Neto et al. (2012) and ter Braak
#' et al. (2012), the p-value to report is the larger of the two tests. With many nodes these p-values need a correction, and
#' `p_adjusted` gives the false discovery rate that the article uses. This correction
#' needs many permutations, as the smallest p-value is 1 / (`n_perm` + 1). The article used
#' 9999, which is advisable with many nodes.
#'
#' The article permutes sites and species at once, which is `p_joint`. This tests the
#' null hypothesis that neither the environment nor the phylogeny matters, so it is not
#' a test of their association. In simulations where the environment acts but species
#' responses have no phylogenetic signal, `p_joint` rejected 76% of the nodes at the 5%
#' level, and `p_max` 4.8%.
#'
#' With few species at a node, permutations of species have little power, so `p_max` is
#' conservative: with the [daphniids] of the article, of 13 species, no predictor passes
#' it, whereas sites alone recover most of the environmental results of its Figure 3 (see
#' the article on validation).
#'
#' Each predictor gets a semipartial correlation, the square root of the loss of R2
#' when it is removed, with the sign of its slope when it is alone in the model.
#'
#' Biogeographic events follow the article: a column of `biogeography` has 1 for sites
#' in one area, -1 for sites in the other and 0 for sites in neither. Each event is
#' tested only in its sites, controlling for `env`. Other contrasts are allowed, for
#' example a column without zeros, which compares one area with all the others. Events
#' with few sites in one area give unstable values. [biogeography_matrix()] builds the
#' matrix from the area of each site.
#'
#' Free permutations treat sites as exchangeable. When sites of an area are neighbours and
#' the community has spatial structure of its own, they reject too often: 45% of tests at
#' the 5% level in a simulation with areas placed at random, against 5% with
#' [torus_permutations()] (see the article on analysing clades). `site_permutations` accepts
#' a restricted null model. For events, each row moves the areas across the sites, which
#' keep their community and predictors. In a landscape with edges, [patch_permutations()]
#' gives the permutations of each event for `area_permutations`.
#'
#' R2, adjusted R2 and the statistic equal those of the MATLAB function
#' `MetacommunityPhylogenetics` of the authors, checked in Octave. Four parts of that
#' function differ from this one:
#' \itemize{
#'   \item In the original, the species permutations shuffle only the first columns
#'     of the community matrix, as many as there are species in the node, so they
#'     do not draw from all species. Here they do.
#'   \item The original keeps the total of the observed data after permuting
#'     species. Here it is recalculated.
#'   \item The semipartial values of the original use the raw predictors, without an
#'     intercept, and multiply by the absolute slope. The article uses standardised
#'     predictors and divides by it, to keep only the sign, which is what we do.
#'   \item The original compares the permuted statistics as they are. The one of each
#'     predictor divides a proportion by a sum of squares, so it depends on the scale of the
#'     response, and with one predictor it does not rank draws as the global one does. Neither
#'     takes into account that draws of species occupy different numbers of sites. Here
#'     draws are ranked on F statistics. `statistic` is still the one of the original, and
#'     p-values of sites alone are the same, but those that use species differ.
#' }
#'
#' As in the original, R2 sums squares without site weights, although the fit
#' uses them, so it is not the R2 of a weighted regression and can exceed 1. The adjusted
#' R2 uses the degrees of freedom of an unweighted model. `r2_weighted` sums squares with the
#' site weights, stays between 0 and 1, and is the value to report as a proportion.
#'
#' Nodes for which a predictor has no variance among the remaining sites, with too
#' few sites, or with a clade without individuals get `NA`. The original stops with an
#' error.
#'
#' Sites are matched by name when `comm`, `env` and `biogeography` all have row names, and
#' otherwise by position.
#'
#' @param comm Matrix with one row per site and one column per species. Species are
#'   matched to the rows of `nodes` by name, in any order, and those in only one of the
#'   two are dropped with a warning.
#' @param nodes Matrix of contrasts with one row per species, from
#'   [node_contrasts()].
#' @param env Predictor, or matrix or data frame of predictors, with one row per
#'   site.
#' @param biogeography Optional matrix of biogeographic events with one row per site
#'   and values of 1, -1 and 0.
#' @param n_perm Number of permutations of sites, of species and of both. With
#'   `site_permutations`, only of species.
#' @param site_permutations Optional matrix of site indices with one permutation per row,
#'   for example from [torus_permutations()], to use instead of free permutations of sites.
#' @param area_permutations Optional list with one matrix per event, from
#'   [patch_permutations()], to use for the events instead of `site_permutations` or free
#'   permutations. Sites must be in the order of `comm`.
#' @param cores Number of cores, with nodes run in parallel. Not available on Windows,
#'   where one core is used. For reproducible p-values with more than one core, call
#'   `RNGkind("L'Ecuyer-CMRG")` before `set.seed()`. Results still depend on the number of
#'   cores.
#' @return List of class `mc_node_test` with `nodes`, a data frame with one row per
#'   node with `node`, `n_species`, `n_sites`, `r2`, `r2_adjusted`, `r2_weighted`, `statistic`,
#'   `p_site`, `p_species`, `p_max`, `p_joint`, `n_draws`, the valid draws of sites, since a
#'   draw without variance does not count, and `p_adjusted`; `semipartial`, a matrix of
#'   nodes by predictors; `p_partial_max`, with the larger of the two p-values for each of
#'   them; `p_partial_site`, with the p-values of sites alone; `p_partial_joint`, with
#'   those of the permutation of the article; `biogeography`, a list with the same
#'   four matrices for the events, or `NULL`; and `permutations`, the numbers of
#'   permutations and whether those of sites were free or restricted.
#' @note Adapted from the MATLAB function `MetacommunityPhylogenetics` in the supplement of
#'   Leibold et al. (2010), with the differences listed in Details.
#' @references
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#'
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of
#' spatial contingency and environmental filtering on metacommunity phylogenetics.
#' *Ecology*, 93, S14--S30. \doi{10.1890/11-0494.1}
#'
#' ter Braak, C. J. F., Cormont, A., & Dray, S. (2012). Improved testing of species
#' traits-environment relationships in the fourth-corner problem. *Ecology*, 93,
#' 1525--1526.
#' @seealso [biogeography_matrix()] codes the events and [patch_permutations()] moves them in
#'   a landscape with edges. [plot_node_test()] and
#'   [plot_node_heatmap()] show the results. [pcc_test()] and
#'   [phylo_beta_test()] test the whole tree at once.
#' @examples
#' # daphniids of Leibold et al. (2010)
#' nodes <- node_contrasts(daphniids$tree)
#' set.seed(1)
#' out <- node_test(daphniids$comm, nodes, daphniids$env, n_perm = 49)
#' out
#' round(out$semipartial, 2)
#' plot_node_test(out, daphniids$tree, test = "site")
#' @export
node_test <- function(comm, nodes, env, biogeography = NULL, n_perm = 999, cores = 1,
                      site_permutations = NULL, area_permutations = NULL) {
  if (inherits(nodes, "phylo")) {
    stop("`nodes` is a tree: use `node_contrasts()` on it first", call. = FALSE)
  }
  comm <- check_comm(comm)
  nodes <- as.matrix(nodes)
  env <- check_env(env, comm)
  stopifnot(n_perm >= 1, cores >= 1)
  check_permutations(site_permutations, comm)
  if (!is.null(biogeography)) biogeography <- check_biogeography(biogeography, comm)
  check_area_permutations(area_permutations, biogeography, comm)
  species <- match_species(comm, rownames(nodes), "`nodes`")
  comm <- comm[, species, drop = FALSE]
  nodes <- nodes[species, , drop = FALSE]
  node_names <- colnames(nodes) %||% seq_len(ncol(nodes))
  env_names <- colnames(env) %||% paste0("predictor", seq_len(ncol(env)))
  # forking is not available on Windows
  workers <- if (.Platform$OS.type == "windows") 1L else cores
  # nodes run in batches, one per worker, so the progress bar can follow them
  per_node <- function(run, name) {
    columns <- seq_len(ncol(nodes))
    batches <- unname(split(columns, (columns - 1) %/% workers))
    results <- progress_map(batches, function(batch) {
      # moves the parent generator on, so each batch gets its own streams
      if (workers > 1) stats::runif(1)
      parallel::mclapply(batch, function(i) run(nodes[, i]), mc.cores = workers)
    }, name)
    results <- unlist(results, recursive = FALSE)
    failed <- vapply(results, inherits, logical(1), "try-error")
    if (any(failed)) {
      stop("a node failed in a parallel worker: ", conditionMessage(attr(results[[which(failed)[1]]],
                                                                         "condition")), call. = FALSE)
    }
    results
  }
  matrix_of <- function(results, field, names) {
    matrix(unlist(lapply(results, `[[`, field)), length(results), byrow = TRUE,
           dimnames = list(node_names, names))
  }
  # one draw of the null model of sites: rows of the community or areas moved across sites
  shuffle_sites <- function(r) {
    rows <- if (is.null(site_permutations)) sample.int(nrow(comm)) else site_permutations[r, ]
    list(comm = comm[rows, , drop = FALSE], env = env)
  }
  shuffle_areas <- function(event, permutations) {
    function(r) {
      moved <- event[permutations[r, ]]
      sites <- moved != 0
      list(comm = comm[sites, , drop = FALSE],
           env = cbind(moved[sites], env[sites, , drop = FALSE]))
    }
  }
  n_site <- if (is.null(site_permutations)) n_perm else nrow(site_permutations)
  results <- per_node(function(contrast) {
    one_node(comm, contrast, env, n_perm, shuffle_sites, n_site)
  }, "Nodes")
  table <- data.frame(node = node_names, do.call(rbind, lapply(results, `[[`, "summary")))
  table$p_adjusted <- stats::p.adjust(table$p_max, "BH")
  structure(list(nodes = table, semipartial = matrix_of(results, "semipartial", env_names),
                 p_partial_max = matrix_of(results, "p_partial_max", env_names),
                 p_partial_site = matrix_of(results, "p_partial_site", env_names),
                 p_partial_joint = matrix_of(results, "p_partial_joint", env_names),
                 biogeography = if (!is.null(biogeography)) {
                   events <- colnames(biogeography) %||% paste0("event", seq_len(ncol(biogeography)))
                   by_event <- lapply(seq_len(ncol(biogeography)), function(e) {
                     sites <- biogeography[, e] != 0
                     predictors <- cbind(biogeography[sites, e], env[sites, , drop = FALSE])
                     local_comm <- comm[sites, , drop = FALSE]
                     permutations <- area_permutations[[e]] %||% site_permutations
                     shuffle <- if (is.null(permutations)) {
                       function(r) {
                         list(comm = local_comm[sample.int(nrow(local_comm)), , drop = FALSE],
                              env = predictors)
                       }
                     } else {
                       shuffle_areas(biogeography[, e], permutations)
                     }
                     n_area <- if (is.null(permutations)) n_perm else nrow(permutations)
                     per_node(function(contrast) {
                       one_node(local_comm, contrast, predictors, n_perm, shuffle, n_area)
                     }, paste("Nodes,", events[e]))
                   })
                   first <- function(field) {
                     vapply(by_event, function(event) {
                       vapply(event, function(node) node[[field]][1], numeric(1))
                     }, numeric(ncol(nodes)))
                   }
                   shaped <- function(x) {
                     matrix(x, ncol(nodes), length(events), dimnames = list(node_names, events))
                   }
                   list(semipartial = shaped(first("semipartial")),
                        p_partial_max = shaped(first("p_partial_max")),
                        p_partial_site = shaped(first("p_partial_site")),
                        p_partial_joint = shaped(first("p_partial_joint")))
                 },
                 permutations = permutation_record(n_perm, site_permutations,
                                                   area_permutations)),
            class = "mc_node_test")
}

#' Matrix of biogeographic events
#'
#' Codes the events of [node_test()] from the area of each site, for example the refugium,
#' basin or island where it lies. Each event gets 1 for sites on one side, -1 for sites on
#' the other and 0 for sites on neither.
#'
#' Events should come from a hypothesis stated in advance. Leibold et al. (2010) contrasted
#' areas on the two sides of historical barriers, and Leibold & Chase (2018, chapter 9) a
#' published biogeographic hypothesis. A named list in `contrast` gives this design, and
#' nested events are written as nested groups of areas. `contrast = "others"` compares each
#' area with the rest, which does not need that hypothesis, but asks a different question
#' and is better seen as exploratory.
#'
#' A site on both sides of an event gets 0.
#'
#' The sites of an area are often neighbours, so their number overstates the information
#' in an event. See `site_permutations` in [node_test()] and [patch_permutations()].
#'
#' @param area Area of each site, with `NA` for sites outside all areas, or a matrix or data
#'   frame with one row per site and one column per area, with 1 or `TRUE` for the sites of
#'   each area. Polygons can give it with `sf::st_join()` or `sf::st_intersects()`.
#' @param contrast Named list of events, each a list with the areas of the two sides, or
#'   `"others"` for one event per area.
#' @param outside For `contrast = "others"`, sites outside all areas are left out of each
#'   event with `"zero"` or count as the other side with `"other"`.
#' @param min_sites A warning is given for events with fewer sites on one side. This is a
#'   rule of thumb, not a value from the literature.
#' @return Matrix with one row per site and one column per event, with values of 1, -1 and 0.
#' @references
#' Leibold, M. A., & Chase, J. M. (2018). *Metacommunity Ecology*. Monographs in Population
#' Biology 59. Princeton University Press. \doi{10.2307/j.ctt1wf4d24}
#'
#' Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
#' phylogenetics: separating the roles of environmental filters and historical
#' biogeography. *Ecology Letters*, 13, 1290--1299.
#' \doi{10.1111/j.1461-0248.2010.01523.x}
#' @examples
#' area <- rep(c("north", "centre", "south", NA), c(8, 6, 6, 10))
#'
#' # a barrier between the north and the rest, and a later one between centre and south
#' events <- list(barrier = list("north", c("centre", "south")),
#'                river = list("centre", "south"))
#' colSums(biogeography_matrix(area, events) == 1)
#'
#' # each area against the others
#' head(biogeography_matrix(area, "others"))
#' @export
biogeography_matrix <- function(area, contrast = "others", outside = c("zero", "other"),
                                min_sites = 5) {
  outside <- match.arg(outside)
  member <- area_membership(area)
  inside <- if (outside == "zero") rowSums(member) > 0 else rep(TRUE, nrow(member))
  sides <- function(first, second) {
    in_any <- function(areas) rowSums(member[, areas, drop = FALSE]) > 0
    in_any(first) * 1 - in_any(second) * 1
  }
  events <- if (identical(contrast, "others")) {
    vapply(colnames(member), function(a) ifelse(member[, a], 1, ifelse(inside, -1, 0)),
           numeric(nrow(member)))
  } else {
    check_events(contrast, colnames(member))
    vapply(contrast, function(event) sides(event[[1]], event[[2]]), numeric(nrow(member)))
  }
  events <- matrix(events, nrow(member), dimnames = list(rownames(member), colnames(events)))
  smaller <- pmin(colSums(events == 1), colSums(events == -1))
  if (any(smaller < min_sites)) {
    warning("fewer than ", min_sites, " sites on one side of: ",
            paste(names(smaller)[smaller < min_sites], collapse = ", "), call. = FALSE)
  }
  events
}

# one logical column per area
area_membership <- function(area) {
  if (is.matrix(area) || is.data.frame(area)) {
    member <- as.matrix(area)
    if (is.null(colnames(member)) || !all(member %in% c(0, 1))) {
      stop("`area` as a matrix needs names of areas and values of 0 and 1", call. = FALSE)
    }
    return(member == 1)
  }
  areas <- sort(unique(as.character(area)))
  if (length(areas) < 2) stop("`area` needs at least two areas", call. = FALSE)
  member <- outer(as.character(area), areas, `==`)
  member[is.na(member)] <- FALSE
  dimnames(member) <- list(names(area), areas)
  member
}

check_events <- function(contrast, areas) {
  well_formed <- is.list(contrast) && !is.null(names(contrast)) &&
    all(lengths(contrast) == 2)
  if (!well_formed) {
    stop("`contrast` must be \"others\" or a named list of events, each with two sides, as ",
         "`list(barrier = list(\"north\", c(\"centre\", \"south\")))`", call. = FALSE)
  }
  unknown <- setdiff(unlist(contrast), areas)
  if (length(unknown) > 0) {
    stop("areas in `contrast` not found in `area`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
}

# events with one row per site, in the order of `comm`
check_biogeography <- function(biogeography, comm) {
  biogeography <- check_sites(as.matrix(biogeography), comm, "`biogeography`")
  if (!all(biogeography %in% c(-1, 0, 1))) {
    stop("`biogeography` needs values of 1, -1 and 0. `biogeography_matrix()` builds it from ",
         "the area of each site", call. = FALSE)
  }
  biogeography
}

check_area_permutations <- function(permutations, biogeography, comm) {
  if (is.null(permutations)) return(invisible())
  if (is.null(biogeography)) {
    stop("`area_permutations` needs `biogeography`", call. = FALSE)
  }
  if (!is.list(permutations) || length(permutations) != ncol(biogeography)) {
    stop("`area_permutations` needs one matrix per event, as given by `patch_permutations()`",
         call. = FALSE)
  }
  for (event in permutations) {
    check_permutations(event, comm, "`area_permutations`")
    if (!is.null(colnames(event)) && !is.null(rownames(comm)) &&
          !identical(colnames(event), rownames(comm))) {
      stop("sites of `area_permutations` are not in the order of `comm`", call. = FALSE)
    }
  }
}

#' @export
print.mc_node_test <- function(x, ...) {
  tested <- x$nodes[!is.na(x$nodes$p_max), ]
  cat("<mc_node_test>", nrow(x$nodes), "nodes,", nrow(tested), "tested\n")
  permutation_line(x$permutations)
  cat("nodes with p (max) <= 0.05:", sum(tested$p_max <= 0.05), "before correction,",
      sum(tested$p_adjusted <= 0.05), "after\n")
  untested <- nrow(x$nodes) - nrow(tested)
  if (untested > 0) {
    cat("nodes not tested:", untested, "(a clade without individuals, too few sites or a",
        "predictor without variance in the sites of the node)\n")
  }
  invisible(x)
}

one_node <- function(comm, contrast, env, n_perm, shuffle, n_site) {
  members <- which(contrast != 0)
  values <- contrast[members]
  observed <- node_fit(comm[, members, drop = FALSE], values, env)
  blank <- rep(NA_real_, ncol(env))
  if (is.na(observed$statistic)) {
    return(list(summary = data.frame(n_species = length(members), n_sites = observed$n_sites,
                                     r2 = NA_real_, r2_adjusted = NA_real_,
                                     r2_weighted = NA_real_, statistic = NA_real_,
                                     p_site = NA_real_, p_species = NA_real_, p_max = NA_real_,
                                     p_joint = NA_real_, n_draws = NA_integer_),
                semipartial = blank, p_partial_max = blank, p_partial_site = blank,
                p_partial_joint = blank))
  }
  null_of <- function(n, draw) {
    t(vapply(seq_len(n), function(r) {
      fit <- draw(r)
      c(fit$f, fit$partial_f)
    }, numeric(ncol(env) + 1)))
  }
  by_site <- null_of(n_site, function(r) {
    moved <- shuffle(r)
    node_fit(moved$comm[, members, drop = FALSE], values, moved$env)
  })
  by_taxa <- null_of(n_perm, function(r) {
    node_fit(comm[, sample.int(ncol(comm), length(members)), drop = FALSE], values, env)
  })
  by_both <- null_of(n_site, function(r) {
    moved <- shuffle(r)
    node_fit(moved$comm[, sample.int(ncol(comm), length(members)), drop = FALSE], values,
             moved$env)
  })
  p_site <- permutation_p(by_site[, 1], observed$f)
  p_species <- permutation_p(by_taxa[, 1], observed$f)
  p_joint <- permutation_p(by_both[, 1], observed$f)
  partial_p <- function(null) {
    vapply(seq_len(ncol(env)), function(j) permutation_p(null[, j + 1], observed$partial_f[j]),
           numeric(1))
  }
  partial_site <- partial_p(by_site)
  partial_taxa <- partial_p(by_taxa)
  partial_joint <- partial_p(by_both)
  list(summary = data.frame(n_species = length(members), n_sites = observed$n_sites,
                            r2 = observed$r2, r2_adjusted = observed$r2_adjusted,
                            r2_weighted = observed$r2_weighted,
                            statistic = observed$statistic, p_site = p_site, p_species = p_species,
                            p_max = max(p_site, p_species), p_joint = p_joint,
                            n_draws = sum(!is.na(by_site[, 1]))),
       semipartial = observed$semipartial, p_partial_max = pmax(partial_site, partial_taxa),
       p_partial_site = partial_site, p_partial_joint = partial_joint)
}

# statistics of one node; all NA when the model cannot be fitted. `statistic` is that of
# the original code. Draws are ranked on F, which is free of the scale of the response and
# takes into account that draws can occupy different numbers of sites
node_fit <- function(incidence, values, env) {
  occupied <- rowSums(incidence) > 0
  incidence <- incidence[occupied, , drop = FALSE]
  env <- env[occupied, , drop = FALSE]
  n <- nrow(incidence)
  p <- ncol(env)
  failed <- list(statistic = NA_real_, f = NA_real_, partial_f = rep(NA_real_, p), n_sites = n)
  site_weights <- rowSums(incidence)
  species_weights <- colSums(incidence)
  both_clades <- any(species_weights[values > 0] > 0) && any(species_weights[values < 0] > 0)
  if (!both_clades || n < p + 2 || any(apply(env, 2, stats::var) == 0)) return(failed)
  standard <- weighted_standardise(matrix(values), species_weights)
  response <- drop(incidence %*% standard)
  predictors <- weighted_standardise(env, site_weights)
  fit <- weighted_fit(predictors, response, site_weights)
  if (is.null(fit)) return(failed)
  total <- sum(response^2)
  explained <- sum(fit^2)
  residual <- sum((response - fit)^2)
  r2 <- explained / total
  partial_r2 <- vapply(seq_len(p), function(j) {
    if (p == 1) return(r2)
    r2 - sum(weighted_fit(predictors[, -j, drop = FALSE], response, site_weights)^2) / total
  }, numeric(1))
  slopes <- colSums(site_weights * predictors * response)
  mean_residual <- residual / (n - p - 1)
  list(statistic = explained / residual, f = explained / p / mean_residual,
       partial_f = partial_r2 * total / mean_residual, r2 = r2,
       r2_weighted = sum(site_weights * fit^2) / sum(site_weights * response^2),
       r2_adjusted = 1 - (n - 1) / (n - p - 1) * (1 - r2),
       semipartial = sign(slopes) * sqrt(pmax(partial_r2, 0)), n_sites = n)
}

weighted_fit <- function(predictors, response, weights) {
  cross <- crossprod(predictors, weights * predictors)
  if (rcond(cross) < .Machine$double.eps) return(NULL)
  drop(predictors %*% solve(cross, crossprod(predictors, weights * response)))
}
