# the class of every column of the tables returned

classes <- function(x) vapply(x, function(v) class(v)[1], character(1))

individuals_data <- function(factors = FALSE) {
  set.seed(1)
  ind <- data.frame(site = rep(paste0("s", 1:6), each = 12),
                    species = rep(c("a", "b", "c"), 24), stringsAsFactors = factors)
  ind$trait <- c(a = 1, b = 2, c = 3)[as.character(ind$species)] + stats::rnorm(72)
  ind$abundance <- 5
  ind
}

test_that("the tables of the functions for individuals have fixed column classes", {
  ind <- individuals_data()
  t1 <- tstats(ind, n_perm = 49)
  expect_identical(classes(t1), c(
    site = "character", trait = "character", n_populations = "integer",
    T_IP.IC = "numeric", T_IC.IR = "numeric", T_PC.PR = "numeric",
    ses_T_IP.IC = "numeric", p_T_IP.IC = "numeric", ses_T_IC.IR = "numeric",
    p_T_IC.IR = "numeric", ses_T_PC.PR = "numeric", p_T_PC.PR = "numeric"))
  s <- tstats_summary(t1)
  expect_identical(classes(s$traits), c(
    trait = "character", statistic = "character", n_sites = "integer",
    percent_low = "numeric", percent_high = "numeric", mean_ses = "numeric",
    p_mean_low = "numeric", p_mean_high = "numeric", sites_low = "character",
    sites_high = "character"))
  expect_identical(classes(s$sites), c(
    site = "character", trait = "character", statistic = "character", observed = "numeric",
    ses = "numeric", p_low = "numeric", p_high = "numeric", deviation = "character"))
  expect_identical(classes(null_intervals(t1)), c(
    site = "character", trait = "character", statistic = "character", observed = "numeric",
    null_median = "numeric", null_low = "numeric", null_high = "numeric",
    inner_low = "numeric", inner_high = "numeric", ses = "numeric", p_low = "numeric",
    p_high = "numeric", deviation = "character"))
  expect_identical(classes(community_index(ind, stats::var, n_perm = 9)), c(
    site = "character", trait = "character", index = "character", observed = "numeric",
    ses = "numeric", p = "numeric"))
  f <- trait_flex(ind)
  expect_identical(classes(f$decomposition), c(
    trait = "character", component = "character", explained = "numeric",
    residual = "numeric", total = "numeric", p = "numeric"))
  expect_identical(classes(f$sites), c(
    site = "character", trait = "character", specific = "numeric", fixed = "numeric",
    intraspecific = "numeric"))
  expect_identical(classes(f$model), c(
    site = "character", trait = "character", component = "character",
    response = "numeric", fitted = "numeric", residual = "numeric"))
  nv <- suppressWarnings(nested_variance(ind, "trait", c("site", "species")))
  expect_identical(classes(nv), c(level = "character", trait = "numeric"))
  expect_identical(classes(trait_densities(ind)), c(
    site = "character", species = "character", level = "character", weight = "numeric",
    trait = "numeric", density = "numeric"))
  expect_identical(classes(trait_rarefaction(ind, n = 2:3, reps = 2)), c(
    effort = "integer", statistic = "character", mean = "numeric", low = "numeric",
    high = "numeric", reps = "integer"))
  expect_identical(classes(as.data.frame(trait_gradient(ind, abundance = "abundance"))), c(
    species = "character", n_sites = "integer", n_alone = "integer", mean = "numeric",
    beta = "numeric", alpha = "numeric", breadth = "numeric", slope = "numeric",
    slope_se = "numeric", intercept = "numeric", min_site = "numeric",
    max_site = "numeric"))
})

test_that("factor labels come back as character columns", {
  ind <- individuals_data(factors = TRUE)
  expect_identical(class(tstats(ind, n_perm = 0)$site), "character")
  expect_identical(class(community_index(ind, stats::var, n_perm = 0)$site), "character")
  expect_identical(class(trait_densities(ind)$species), "character")
})

test_that("as.data.frame() of node_test() has fixed column classes", {
  set.seed(1)
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  comm <- matrix(stats::rpois(40, 3), 10, 4,
                 dimnames = list(paste0("s", 1:10), c("a", "b", "c", "d")))
  out <- node_test(comm, node_contrasts(tree), stats::rnorm(10), n_perm = 3)
  expect_identical(classes(as.data.frame(out)), c(
    node = "character", n_species = "integer", n_sites = "integer", r2 = "numeric",
    r2_adjusted = "numeric", r2_original = "numeric",
    r2_adjusted_original = "numeric", r2_weighted = "numeric", statistic = "numeric",
    p_site = "numeric", p_species = "numeric", p_max = "numeric", p_joint = "numeric",
    n_draws = "integer", p_adjusted = "numeric"))
})

test_that("list columns stop with a message that names them", {
  ind <- individuals_data()
  listed <- ind
  listed$trait <- as.list(listed$trait)
  expect_error(tstats(listed, n_perm = 0), "columns in `traits` must be numeric")
  expect_error(trait_gradient(listed), "`trait` must be numeric")
  comm <- matrix(1:12, 4, 3, dimnames = list(paste0("s", 1:4), c("a", "b", "c")))
  table <- as.data.frame(comm)
  table$b <- as.list(table$b)
  expect_error(check_comm(table), "columns of `comm` that are not numeric: b")
  env <- data.frame(x = 1:4, y = I(as.list(4:1)), row.names = paste0("s", 1:4))
  expect_error(check_env(env, comm), "columns of `env` that are not numeric: y")
})

test_that("columns of class units give the same results as plain columns", {
  skip_if_not_installed("units")
  ind <- individuals_data()
  measured <- ind
  measured$trait <- units::set_units(measured$trait, "mm")
  expect_equal(tstats(measured, n_perm = 0), tstats(ind, n_perm = 0))
  expect_equal(as.data.frame(trait_gradient(measured)), as.data.frame(trait_gradient(ind)))
  comm <- matrix(1:12, 4, 3, dimnames = list(paste0("s", 1:4), c("a", "b", "c")))
  env <- data.frame(x = c(1, 4, 2, 8), row.names = paste0("s", 1:4))
  env_units <- env
  env_units$x <- units::set_units(env_units$x, "degC")
  expect_equal(check_env(env_units, comm), check_env(env, comm))
})
