simulated <- function(lambda, n_sites = 50, n_species = 14) {
  tree <- ape::rcoal(n_species)
  traits <- species_traits(tree, lambda = lambda, optima_range = c(-1, 101),
                           rmax_range = c(0, 30))
  env <- runif(n_sites, 0, 100)
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  keep <- rowSums(comm) > 0
  list(comm = comm[keep, ], tree = tree, env = env[keep])
}

test_that("contrasts halve at each node, as worked out by hand", {
  tree <- ape::read.tree(text = "((a,b),(c,(d,e)));")
  expect_equal(ncol(node_contrasts(tree)), 4)
  expect_equal(unname(node_contrasts(tree)[, "7"]), c(-1, 1, 0, 0, 0))
  x <- node_contrasts(tree, cherries = FALSE)
  # nodes: root, (a,b) is left out, (c,(d,e)), (d,e) is left out
  expect_equal(ncol(x), 2)
  expect_equal(x[, 1], c(a = -0.5, b = -0.5, c = 0.5, d = 0.25, e = 0.25))
  expect_equal(x[, 2], c(a = 0, b = 0, c = -1, d = 0.5, e = 0.5))
  expect_error(node_contrasts(ape::read.tree(text = "(a,b,c);")), "dichotomous")
})

test_that("statistics equal the original function run in Octave", {
  # expected values from MetacommunityPhylogenetics in Octave 11.3, see
  # data-raw/validation_node.R
  set.seed(11)
  n <- 40
  k <- 12
  inc <- matrix(rbinom(n * k, 1, 0.45), n)
  inc[rowSums(inc) == 0, 1] <- 1
  env <- cbind(runif(n, 10, 30), rnorm(n, 5, 2))
  nodes <- cbind(c(-.5, -.5, .5, .25, .25, rep(0, 7)),
                 c(rep(0, 6), -.5, -.5, .5, .25, .125, .125),
                 c(-.25, -.25, -.25, -.125, -.125, .5, .25, .25, .125, .125, .0625, .0625))
  expected <- rbind(c(0.00783339, -0.05229792, 0.00799555),
                    c(0.03382408, -0.02138597, 0.03483741),
                    c(0.07671692, 0.02680973, 0.07970519))
  fits <- lapply(1:3, function(i) {
    m <- which(nodes[, i] != 0)
    node_fit(inc[, m, drop = FALSE], nodes[m, i], env)
  })
  mine <- t(vapply(fits, function(f) c(f$r2, f$r2_adjusted, f$statistic), numeric(3)))
  expect_equal(mine, expected, tolerance = 1e-6)
})

test_that("with one predictor, the semipartial squared is the R2", {
  set.seed(1)
  d <- simulated(1)
  x <- node_contrasts(d$tree)
  out <- node_test(d$comm, x, d$env, n_perm = 3)
  expect_equal(unname(out$semipartial[, 1]^2), out$nodes$r2)
  flipped <- node_test(d$comm, x, -d$env, n_perm = 3)
  expect_equal(flipped$semipartial, -out$semipartial)
  expect_equal(out$nodes$p_adjusted, p.adjust(out$nodes$p_max, "BH"))
})

test_that("several cores give the same statistics", {
  x <- node_contrasts(daphniids$tree)
  one <- node_test(daphniids$comm, x, daphniids$env, n_perm = 2)
  two <- node_test(daphniids$comm, x, daphniids$env, n_perm = 2, cores = 2)
  expect_equal(two$nodes$r2, one$nodes$r2)
  expect_equal(two$semipartial, one$semipartial)
  expect_error(node_test(daphniids$comm, x, daphniids$env, n_perm = 2, cores = 0))
})

test_that("the joint test breaks both links", {
  set.seed(4)
  x <- node_contrasts(daphniids$tree)
  out <- node_test(daphniids$comm, x, daphniids$env, n_perm = 9)
  expect_true(all(out$nodes$p_joint > 0, na.rm = TRUE))
  expect_equal(dim(out$p_partial_joint), dim(out$p_partial_max))
})

test_that("the daphniids run from the data of the package", {
  set.seed(1)
  out <- node_test(daphniids$comm, node_contrasts(daphniids$tree), daphniids$env,
                   biogeography = daphniids$biogeography[, 1:2], n_perm = 2)
  expect_equal(dim(out$p_partial_site), c(12, 9))
  expect_true(all(out$p_partial_max >= out$p_partial_site, na.rm = TRUE))
  expect_equal(dim(out$biogeography$p_partial_site), c(12, 2))
})

test_that("biogeographic events use only their sites and control for the environment", {
  set.seed(7)
  d <- simulated(1, n_sites = 70)
  x <- node_contrasts(d$tree)
  n <- nrow(d$comm)
  event <- sample(c(-1, 0, 1), n, replace = TRUE)
  out <- node_test(d$comm, x, d$env, biogeography = cbind(north_south = event), n_perm = 3)
  expect_identical(colnames(out$biogeography$semipartial), "north_south")
  expect_equal(dim(out$biogeography$p_partial_max), c(ncol(x), 1))
  sites <- event != 0
  m <- which(x[, 1] != 0)
  direct <- node_fit(d$comm[sites, names(m), drop = FALSE], x[m, 1],
                     cbind(event[sites], d$env[sites]))
  expect_equal(unname(out$biogeography$semipartial[1, 1]), direct$semipartial[1])
  unnamed <- node_test(d$comm, x, d$env, biogeography = cbind(event), n_perm = 2)
  expect_null(node_test(d$comm, x, d$env, n_perm = 2)$biogeography)
  expect_identical(colnames(unname(unnamed$biogeography$semipartial)), NULL)
  expect_error(node_test(d$comm, x, d$env, biogeography = cbind(event * 2), n_perm = 2))
})

test_that("the test detects clades that differ in their response", {
  set.seed(2)
  d <- simulated(1, n_sites = 80)
  out <- node_test(d$comm, node_contrasts(d$tree), d$env, n_perm = 99)
  expect_s3_class(out, "mc_node_test")
  expect_equal(nrow(out$nodes), ncol(node_contrasts(d$tree)))
  expect_lt(min(out$nodes$p_site, na.rm = TRUE), 0.05)
  expect_equal(out$nodes$p_max, pmax(out$nodes$p_site, out$nodes$p_species))
  expect_output(print(out), "tested")
})

test_that("species permutations draw from all species", {
  # the first species carry no signal: the original would shuffle only them
  set.seed(3)
  comm <- matrix(rbinom(60 * 12, 1, 0.5), 60, dimnames = list(NULL, paste0("s", 1:12)))
  env <- seq(0, 1, length.out = 60)
  comm[, 11] <- as.integer(env > 0.5)
  comm[, 12] <- as.integer(env <= 0.5)
  nodes <- matrix(0, 12, 1, dimnames = list(colnames(comm), "n"))
  nodes[11:12, 1] <- c(-1, 1)
  out <- node_test(comm, nodes, env, n_perm = 99)
  expect_lte(out$nodes$p_max, 0.02)
})

test_that("several predictors give one partial value each", {
  set.seed(4)
  d <- simulated(1)
  env <- cbind(temperature = d$env, noise = runif(nrow(d$comm)))
  out <- node_test(d$comm, node_contrasts(d$tree), env, n_perm = 5)
  expect_identical(colnames(out$semipartial), colnames(env))
  expect_equal(dim(out$p_partial_max), dim(out$semipartial))
  unnamed <- node_test(d$comm, node_contrasts(d$tree), unname(env), n_perm = 2)
  expect_identical(colnames(unnamed$semipartial), c("predictor1", "predictor2"))
})

test_that("nodes that cannot be fitted give NA", {
  set.seed(5)
  comm <- matrix(rbinom(30 * 4, 1, 0.6), 30, dimnames = list(NULL, paste0("s", 1:4)))
  comm[, 3:4] <- 0
  comm[1:2, 3] <- 1
  comm[1:2, 4] <- c(1, 0)
  env <- cbind(c(1, 1, runif(28)), c(rep(1, 2), rep(2, 28)))
  nodes <- cbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  rownames(nodes) <- colnames(comm)
  out <- node_test(comm, nodes, env, n_perm = 2)
  expect_true(is.na(out$nodes$p_max[2]))
  expect_false(is.na(out$nodes$r2[1]))
  # one species only: no contrast within the node
  single <- cbind(c(1, 1, 0, 0))
  rownames(single) <- colnames(comm)
  expect_true(is.na(node_test(comm, single, runif(30), n_perm = 2)$nodes$r2))
  # collinear predictors
  expect_true(is.na(node_fit(comm[, 1:2], c(1, -1), cbind(1:30, 2 * (1:30)))$statistic))
  expect_null(weighted_fit(cbind(1:5, 1:5), 1:5, rep(1, 5)))
})

test_that("invalid data are rejected", {
  set.seed(6)
  d <- simulated(1)
  x <- node_contrasts(d$tree)
  expect_warning(node_test(d$comm[, -1], x, d$env, n_perm = 2), "not in `comm` dropped")
  expect_error(node_test(d$comm, x, d$env[-1], n_perm = 5))
})

test_that("branch lengths give the values of Figure S2 of the supplement", {
  tree <- ape::read.tree(text = "((((A:0.2,B:0.2):0.65,C:0.85):0.15,D:1):0.4,(E:0.5,F:0.5):0.9);")
  x <- node_contrasts(tree, branch_lengths = TRUE)
  expect_equal(abs(x[c("A", "B", "C", "D", "E", "F"), "7"]),
               c(A = 0.196, B = 0.196, C = 0.290, D = 0.319, E = 0.5, F = 0.5), tolerance = 2e-3)
  expect_equal(abs(x[c("A", "B", "C", "D"), "8"]), c(A = 0.281, B = 0.281, C = 0.439, D = 1),
               tolerance = 2e-3)
  expect_equal(unname(abs(x[c("A", "B", "C"), "9"])), c(0.5, 0.5, 1))
  sums <- apply(x, 2, function(column) c(sum(column[column < 0]), sum(column[column > 0])))
  expect_true(all(abs(sums) - 1 < 1e-12))
  expect_equal(dimnames(x), dimnames(node_contrasts(tree)))
  expect_error(node_contrasts(ape::read.tree(text = "((a,b),c);"), branch_lengths = TRUE),
               "no branch lengths")
})

test_that("the events of the daphniids come out of the areas", {
  observed <- daphniids$biogeography
  area <- apply(observed, 1, paste, collapse = " ")
  events <- lapply(colnames(observed), function(e) {
    list(unique(area[observed[, e] == 1]), unique(area[observed[, e] == -1]))
  })
  names(events) <- colnames(observed)
  built <- biogeography_matrix(area, events)
  expect_equal(unname(built), unname(observed))
  expect_equal(colnames(built), colnames(observed))
})

test_that("each area can face the others, with sites outside left out or not", {
  area <- c(a = "x", b = "x", c = "y", d = "z", e = NA)
  zero <- biogeography_matrix(area, "others", min_sites = 1)
  expect_equal(zero[, "x"], c(a = 1, b = 1, c = -1, d = -1, e = 0))
  other <- biogeography_matrix(area, "others", outside = "other", min_sites = 1)
  expect_equal(other[, "y"], c(a = -1, b = -1, c = 1, d = -1, e = -1))
})

test_that("a matrix of areas gives the same events, with overlaps as zero", {
  area <- rep(c("x", "y", NA), c(3, 3, 2))
  member <- cbind(x = area %in% "x", y = area %in% "y")
  expect_equal(biogeography_matrix(member, min_sites = 1),
               biogeography_matrix(area, min_sites = 1))
  member[1, "y"] <- TRUE
  events <- list(e = list("x", "y"))
  expect_equal(biogeography_matrix(member, events, min_sites = 1)[, "e"],
               c(0, 1, 1, -1, -1, -1, 0, 0))
})

test_that("mistakes in areas and events are explained", {
  area <- rep(c("x", "y"), each = 6)
  expect_error(biogeography_matrix(rep("x", 6)), "two areas")
  expect_error(biogeography_matrix(matrix(1:4, 2)), "names of areas")
  expect_error(biogeography_matrix(area, "pairs"), "named list")
  expect_error(biogeography_matrix(area, list(e = list("x"))), "named list")
  expect_error(biogeography_matrix(area, list(e = list("x", "w"))), "not found in `area`: w")
  expect_warning(biogeography_matrix(area, min_sites = 7), "fewer than 7 sites on one side of: x, y")
})

test_that("restricted permutations reach sites and events", {
  set.seed(1)
  nodes <- node_contrasts(daphniids$tree)
  rows <- which(rowSums(daphniids$comm) > 0)[1:60]
  comm <- daphniids$comm[rows, ]
  env <- daphniids$env[rows, 1:2]
  events <- cbind(e = rep(c(1, -1, 0), each = 20))
  perms <- t(replicate(9, sample(60)))
  out <- node_test(comm, nodes, env, events, n_perm = 9, site_permutations = perms)
  tested <- !is.na(out$nodes$p_site)
  expect_true(all(out$nodes$p_site[tested] >= 0.1))
  expect_true(any(!is.na(out$biogeography$p_partial_site)))
  # identical rows leave the data as they are, so the observed value is never beaten
  same <- matrix(seq_len(60), 3, 60, byrow = TRUE)
  fixed <- node_test(comm, nodes, env, events, n_perm = 3, site_permutations = same)
  expect_equal(fixed$nodes$p_site[tested], rep(1, sum(tested)))
  expect_error(node_test(comm, nodes, env, n_perm = 3, site_permutations = same[, -1]))
})

test_that("with one predictor, the global and partial tests rank draws alike", {
  set.seed(7)
  nodes <- node_contrasts(daphniids$tree)
  out <- node_test(daphniids$comm, nodes, daphniids$env[, "PHEQ", drop = FALSE], n_perm = 19)
  expect_equal(unname(out$p_partial_max[, 1]), out$nodes$p_max)
  expect_equal(unname(out$p_partial_joint[, 1]), out$nodes$p_joint)
})

test_that("the partial statistic does not depend on the scale of the response", {
  set.seed(8)
  comm <- matrix(rbinom(240, 1, 0.5), 60)
  env <- cbind(runif(60), runif(60))
  one <- node_fit(comm, c(1, 1, -1, -1), env)
  two <- node_fit(cbind(comm, comm), rep(c(1, 1, -1, -1), 2), env)
  expect_equal(one$partial_f, two$partial_f)
  expect_equal(one$f, two$f)
})

test_that("a node with an absent clade is not tested", {
  comm <- matrix(rbinom(120, 1, 0.5), 30, dimnames = list(NULL, letters[1:4]))
  comm[, "d"] <- 0
  env <- runif(30)
  expect_true(is.na(node_fit(comm[, 3:4], c(1, -1), cbind(env))$statistic))
  expect_false(is.na(node_fit(comm[, 1:2], c(1, -1), cbind(env))$statistic))
})

test_that("mistakes in the arguments of node_test() come early and explained", {
  nodes <- node_contrasts(daphniids$tree)
  expect_error(node_test(daphniids$comm, daphniids$tree, daphniids$env), "use `node_contrasts")
  expect_error(node_test(daphniids$comm, nodes, daphniids$env, 999), "1 rows and `comm`")
  expect_error(node_test(daphniids$comm, nodes, daphniids$env, cbind(rep("a", 314))),
               "biogeography_matrix")
})

test_that("each parallel call gets its own streams", {
  skip_on_os("windows")
  kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(kind[1]))
  set.seed(9)
  nodes <- node_contrasts(daphniids$tree)
  env <- daphniids$env[, 1:2]
  one <- node_test(daphniids$comm, nodes, env, n_perm = 9, cores = 2)$nodes$p_site
  two <- node_test(daphniids$comm, nodes, env, n_perm = 9, cores = 2)$nodes$p_site
  expect_false(identical(one, two))
})

test_that("a test without valid draws has no p-value", {
  expect_true(is.na(permutation_p(c(NA_real_, NA_real_), 1)))
})

test_that("patch_permutations() keeps the number and size of patches", {
  set.seed(1)
  coords <- expand.grid(x = 1:12, y = 1:8)
  event <- rep(0, 96)
  event[coords$x <= 3 & coords$y <= 3] <- 1
  event[coords$x == 12 & coords$y <= 2] <- 1
  event[coords$x >= 7 & coords$x <= 9 & coords$y >= 6] <- -1
  drawn <- patch_permutations(coords, cbind(a = event), n_perm = 30)
  expect_named(drawn, "a")
  expect_equal(dim(drawn$a), c(30, 96))
  expect_true(all(apply(drawn$a, 1, function(row) all(sort(row) == 1:96))))
  moved <- apply(drawn$a, 1, function(row) event[row])
  expect_true(all(colSums(moved == 1) == 11))
  expect_true(all(colSums(moved == -1) == 9))
  # places change among permutations
  expect_gt(nrow(unique(t(moved))), 20)
  # a patch of one site is kept
  lone <- rep(0, 96)
  lone[c(1, 50)] <- c(1, -1)
  expect_true(all(apply(patch_permutations(coords, cbind(lone), 5)$lone, 1,
                        function(row) sum(lone[row] != 0)) == 2))
})

test_that("patches are found with the distance between sites", {
  coords <- cbind(c(1, 2, 3, 10, 11), 0)
  d <- as.matrix(dist(coords))
  expect_equal(lengths(area_patches(rep(TRUE, 5), d, 1.5)), c(3, 2))
  expect_length(area_patches(rep(TRUE, 5), d, 8), 1)
  expect_length(area_patches(rep(FALSE, 5), d, 1.5), 0)
})

test_that("a side without other sites takes what is left", {
  set.seed(2)
  coords <- expand.grid(x = 1:6, y = 1:6)
  event <- ifelse(coords$x <= 2 & coords$y <= 2, 1, -1)
  drawn <- patch_permutations(coords, cbind(event), n_perm = 10, gap = 1.5)$event
  expect_true(all(apply(drawn, 1, function(row) sum(event[row] == 1)) == 4))
})

test_that("patch_permutations() rejects wrong input", {
  coords <- expand.grid(1:3, 1:3)
  event <- cbind(rep(c(1, 0, -1), 3))
  expect_error(patch_permutations(coords[-1, ], event), "rows")
  expect_error(patch_permutations(cbind(coords, 1), event), "two numeric")
  expect_error(patch_permutations(coords, event * 2), "values of 1")
  expect_error(patch_permutations(data.frame(a = letters[1:9], b = 1:9), event), "not numeric")
})

test_that("node_test() uses area_permutations for events only", {
  set.seed(3)
  tree <- ape::rcoal(6)
  nodes <- node_contrasts(tree)
  coords <- expand.grid(x = 1:6, y = 1:5)
  comm <- matrix(rpois(180, 3), 30, dimnames = list(paste0("s", 1:30), tree$tip.label))
  env <- cbind(a = rnorm(30))
  events <- cbind(e = ifelse(coords$x <= 2, 1, ifelse(coords$x >= 5, -1, 0)))
  rownames(events) <- rownames(comm)
  drawn <- patch_permutations(coords, events, n_perm = 7)
  out <- node_test(comm, nodes, env, events, n_perm = 9, area_permutations = drawn)
  expect_equal(out$permutations$areas, 7)
  expect_equal(out$permutations$null, "free")
  expect_output(print(out), "7 of areas")
  expect_gte(min(out$biogeography$p_partial_site, na.rm = TRUE), 1 / 8)
  expect_error(node_test(comm, nodes, env, n_perm = 3, area_permutations = drawn),
               "needs `biogeography`")
  expect_error(node_test(comm, nodes, env, events, n_perm = 3, area_permutations = drawn[[1]]),
               "one matrix per event")
  expect_error(node_test(comm, nodes, env, events, n_perm = 3,
                         area_permutations = list(drawn[[1]][, -1])), "`area_permutations` needs")
  rownames(comm) <- rev(rownames(comm))
  rownames(env) <- rownames(comm)
  expect_error(node_test(comm, nodes, env, events, n_perm = 3, area_permutations = drawn),
               "order of `comm`")
})

test_that("node_test() has a weighted R2 and counts draws", {
  set.seed(5)
  nodes <- node_contrasts(daphniids$tree)
  out <- node_test(daphniids$comm, nodes, daphniids$env[, 1:2], n_perm = 5)
  tested <- out$nodes[!is.na(out$nodes$r2), ]
  expect_true(all(tested$r2_weighted >= 0 & tested$r2_weighted <= 1))
  expect_true(all(tested$n_draws <= 5))
  expect_equal(nrow(as.data.frame(out)), ncol(nodes))
  expect_false(is.unsorted(summary(out)$p_adjusted))
})

test_that("clades with branches of length zero have equal weights", {
  tree <- ape::read.tree(text = "((a:0,b:0):0,(c:1,d:1):1);")
  contrasts <- node_contrasts(tree, branch_lengths = TRUE)
  expect_equal(unname(contrasts[c("a", "b"), 1]), c(-0.5, -0.5))
})

test_that("a failure in a worker stops node_test()", {
  skip_on_os("windows")
  nodes <- node_contrasts(daphniids$tree)
  local_mocked_bindings(one_node = function(...) stop("broken"))
  expect_error(suppressWarnings(node_test(daphniids$comm, nodes, daphniids$env[, 1:2],
                                          n_perm = 2, cores = 2)), "broken")
})

test_that("result methods give tables and plots", {
  set.seed(6)
  tree <- ape::rcoal(6)
  comm <- matrix(rpois(180, 2) + 1, 30, dimnames = list(NULL, tree$tip.label))
  env <- rnorm(30)
  pcc <- pcc_test(comm, tree, env, n_perm = 3)
  beta <- phylo_beta_test(comm, tree, env, n_perm = 3)
  expect_equal(as.data.frame(pcc), summary(pcc))
  expect_equal(as.data.frame(beta)$p_max, summary(beta)$p_max)
  vector <- phylo_vector(comm, tree, env, n_perm = 3, level = 0.5, n_sim = 3)
  expect_equal(as.data.frame(vector), summary(vector))
  traits <- cbind(t1 = rnorm(6))
  rownames(traits) <- tree$tip.label
  part <- trait_partition(comm, traits, env, rnorm(30), n_perm = 3)
  expect_equal(as.data.frame(part), summary(part))
  expect_silent(plot(part))
  expect_warning(plot(part, test = "taxa"), "use `\"species\"`")
  out <- node_test(comm, node_contrasts(tree), env, n_perm = 3)
  expect_silent(plot(out, tree, level = 0.5))
})

test_that("the heatmap handles repeated names and empty results", {
  set.seed(7)
  tree <- ape::rcoal(5)
  comm <- matrix(rpois(100, 3), 20, dimnames = list(NULL, tree$tip.label))
  env <- cbind(same = rnorm(20))
  events <- cbind(same = rep(c(1, -1), 10))
  out <- node_test(comm, node_contrasts(tree), env, events, n_perm = 3)
  expect_equal(ncol(plot_node_heatmap(out, level = 1)), 2)
  out$semipartial[] <- NA
  out$biogeography$semipartial[] <- NA
  expect_error(plot_node_heatmap(out), "no node was tested")
})
