toy_individuals <- function(seed = 1, n_sites = 6, n = 120) {
  set.seed(seed)
  d <- data.frame(site = rep(sprintf("s%02d", seq_len(n_sites)), length.out = n),
                  species = sample(c("a", "b", "c", "d"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  d$trait <- c(a = 10, b = 20, c = 30, d = 40)[d$species] + rnorm(n)
  d$other <- rnorm(n)
  d
}

test_that("tstats_summary() lists deviating sites and skips untestable ones", {
  d <- toy_individuals()
  d <- rbind(d, data.frame(site = "solo", species = "a", trait = 10 + rnorm(8),
                           other = rnorm(8)))
  set.seed(2)
  res <- tstats(d, traits = c("trait", "other"), nperm = 49)
  s <- tstats_summary(res)
  expect_named(s, c("traits", "sites"))
  expect_equal(nrow(s$traits), 6)
  solo <- s$sites[s$sites$site == "solo" & s$sites$statistic == "T_IP.IC", ]
  expect_true(all(is.na(solo$p_low) & is.na(solo$p_high) & is.na(solo$deviation)))
  ip <- s$traits[s$traits$statistic == "T_IP.IC" & s$traits$trait == "trait", ]
  expect_equal(ip$n_sites, 6)
  expect_false(grepl("solo", ip$sites_low))
  # site lists agree with the per-site table
  for (i in seq_len(nrow(s$traits))) {
    row <- s$traits[i, ]
    rows <- s$sites[s$sites$trait == row$trait & s$sites$statistic == row$statistic, ]
    low <- rows$site[which(rows$deviation == "low")]
    expect_identical(row$sites_low, paste(low, collapse = " "))
    expect_equal(row$percent_low, 100 * length(low) / row$n_sites)
  }
  # T_IP.IC is far below the null of permutation within sites when species differ
  expect_equal(ip$percent_low, 100)
  expect_lt(ip$p_mean_low, 0.05)
  expect_error(tstats_summary(tstats(d, nperm = 0)), "nperm")
})

test_that("trait_flex() splits the community mean as Leps et al. (2011)", {
  d <- toy_individuals()
  env <- data.frame(site = rev(sprintf("s%02d", 1:6)), moisture = 6:1)
  d$trait <- d$trait + env$moisture[match(d$site, env$site)]
  res <- trait_flex(d, formula = ~moisture, env = env)
  dec <- res$decomposition
  expect_equal(dec$component, c("turnover", "intraspecific", "covariation", "total"))
  expect_equal(sum(dec$total[1:3]), dec$total[4])
  expect_equal(dec$explained + dec$residual, dec$total)
  # hand computation of the specific and fixed means of one site
  s1 <- d[d$site == "s01", ]
  pop_mean <- tapply(s1$trait, s1$species, mean)
  n <- table(s1$species)[names(pop_mean)]
  fixed <- tapply(d$trait, d$species, mean)[names(pop_mean)]
  site1 <- res$sites[res$sites$site == "s01", ]
  expect_equal(site1$specific, sum(n * pop_mean) / sum(n))
  expect_equal(site1$fixed, sum(n * fixed) / sum(n))
  # env is matched by name, not by order
  shuffled <- trait_flex(d, formula = ~moisture, env = env[c(3, 1, 6, 2, 5, 4), ])
  expect_equal(shuffled$decomposition, dec)
  expect_lt(dec$p[dec$component == "intraspecific"], 0.05)
  # without a model nothing is explained
  plain <- trait_flex(d)$decomposition
  expect_true(all(plain$explained == 0) && all(is.na(plain$p)))
})

test_that("trait_flex() uses abundances and handles missing values consistently", {
  d <- toy_individuals()
  d$abundance <- c(a = 1, b = 5, c = 2, d = 10)[d$species]
  site1 <- trait_flex(d, abundance = "abundance")$sites[1, ]
  s1 <- d[d$site == "s01", ]
  pm <- tapply(s1$trait, s1$species, mean)
  w <- c(a = 1, b = 5, c = 2, d = 10)[names(pm)]
  expect_equal(site1$specific, sum(w * pm) / sum(w))
  # a population with all values missing is dropped from both community means
  gone <- d$site == "s02" & d$species == "a"
  d$trait[gone] <- NA
  res <- trait_flex(d)$sites
  s2 <- d[d$site == "s02" & !is.na(d$trait), ]
  n <- table(s2$species)
  fixed <- tapply(d$trait, d$species, mean, na.rm = TRUE)[names(n)]
  expect_equal(res$fixed[res$site == "s02"], sum(n * fixed) / sum(n))
  expect_error(trait_flex(d, formula = ~moisture), "env")
  expect_error(trait_flex(d, formula = ~moisture, env = data.frame(site = "s01", moisture = 1)),
               "lacks")
  expect_error(trait_flex(d, traits = "nope"), "no column")
})

test_that("nested_variance() checks the nesting and sums to one", {
  skip_if_not_installed("nlme")
  d <- toy_individuals(n_sites = 8, n = 240)
  d$region <- ifelse(d$site %in% sprintf("s%02d", 1:4), "north", "south")
  d$population <- paste(d$site, d$species)
  res <- suppressWarnings(nested_variance(d, c("trait", "other"),
                                          c("region", "site", "population")))
  expect_equal(res$level, c("region", "site", "population", "within"))
  expect_equal(colSums(res[c("trait", "other")]), c(trait = 1, other = 1))
  expect_gt(res$trait[res$level == "population"], 0.9)
  expect_error(nested_variance(d, "trait", c("region", "population", "site")),
               "does not split")
  expect_warning(nested_variance(d, "trait", c("region", "species")), "crossed")
  expect_error(nested_variance(d, "trait", c("region", "nope")), "no column")
})

test_that("nested_variance() agrees with cati::partvar()", {
  skip_if_not_installed("nlme")
  skip_if_not_installed("cati")
  d <- toy_individuals(n_sites = 8, n = 240)
  d$region <- ifelse(d$site %in% sprintf("s%02d", 1:4), "north", "south")
  ours <- nested_variance(d, c("trait", "other"), c("region", "site"))
  theirs <- suppressMessages(utils::capture.output(
    pv <- cati::partvar(traits = d[c("trait", "other")],
                        factors = data.frame(region = d$region, site = d$site))))
  expect_equal(unname(as.matrix(ours[c("trait", "other")])), unname(unclass(pv)[, 1:2]),
               tolerance = 1e-6)
})

test_that("community_index() computes observed values and keeps the null models apart", {
  d <- toy_individuals()
  res <- community_index(d, stats::var, null = "local", nperm = 9)
  expect_equal(res$observed, as.vector(tapply(d$trait, d$site, stats::var)[unique(d$site)]))
  # permuting within sites cannot change an index that ignores species
  expect_true(all(abs(attr(res, "null") - rep(res$observed, each = 9)) < 1e-9))
  expect_true(all(is.na(res$ses)))
  # an index of species sees the permutation
  pop_var <- function(x, species) stats::var(tapply(x, species, mean))
  sp <- community_index(d, pop_var, null = "local", nperm = 19)
  # permuting within sites mixes the species, so the variance of their means falls
  expect_true(all(sp$ses > 0))
  # regional draws stay within the values of the data
  reg <- community_index(d, mean, null = "regional.ind", nperm = 19)
  expect_true(all(attr(reg, "null") >= min(d$trait) & attr(reg, "null") <= max(d$trait)))
  # the order of the rows does not matter for observed values
  shuffled <- community_index(d[sample(nrow(d)), ], pop_var, null = "local", nperm = 0)
  expect_equal(shuffled$observed[match(sp$site, shuffled$site)], sp$observed)
  expect_named(shuffled, c("site", "trait", "observed"))
})

test_that("community_index() handles populations, several traits and checks the index", {
  d <- toy_individuals()
  prab <- community_index(d, length, null = "regional.pop.prab", nperm = 5)
  expect_equal(prab$observed,
               as.vector(tapply(d$species, d$site, function(s) length(unique(s)))[unique(d$site)]))
  pop <- community_index(d, stats::var, null = "regional.pop", nperm = 5)
  means <- unique(stats::ave(d$trait, d$site, d$species))
  expect_true(all(round(attr(pop, "null"), 8) >= 0))
  two <- community_index(d, mean, traits = c("trait", "other"), null = "regional.ind", nperm = 5)
  expect_equal(unique(two$trait), c("trait", "other"))
  expect_equal(dim(attr(two, "null")), c(5, 12))
  # multivariate: rows move together, so the correlation within a site is kept by "local"
  d$twin <- d$trait * 2
  mv <- community_index(d, function(x) stats::cor(x[, 1], x[, 2]), traits = c("trait", "twin"),
                        null = "local", nperm = 5, multivariate = TRUE)
  expect_equal(unique(mv$trait), "multivariate")
  expect_true(all(abs(attr(mv, "null") - 1) < 1e-9))
  expect_error(community_index(d, range), "single number")
  expect_error(community_index(d, mean, traits = "nope"), "no column")
})

test_that("community_index() and trait_flex() agree with cati where cati is right", {
  skip_if_not_installed("cati")
  d <- toy_individuals()
  d <- d[order(d$site), ]
  ours <- community_index(d, stats::var, traits = c("trait", "other"), nperm = 0)
  theirs <- cati::ComIndex(as.matrix(d[c("trait", "other")]),
                           index = "tapply(x, ind.plot, var)", sp = factor(d$species),
                           ind.plot = factor(d$site), nullmodels = "regional.ind", nperm = 2,
                           printprogress = FALSE)
  expect_equal(unname(cbind(ours$observed[ours$trait == "trait"],
                            ours$observed[ours$trait == "other"])),
               unname(theirs$obs[[1]]))
  # decompCTRE() looks for the variables of its formula outside the function, so it is
  # compared without a model
  dc <- cati::decompCTRE(traits = d["trait"], sp = factor(d$species),
                         ind.plot = factor(d$site), printprogress = FALSE)
  ours <- trait_flex(d)$decomposition
  expect_equal(ours$total, unname(unlist(dc$trait$SumSq[1, ])), tolerance = 1e-8)
})

test_that("the regional nulls of tstats() keep the observed denominators", {
  set.seed(3)
  d <- data.frame(site = rep(sprintf("s%d", 1:4), each = 15),
                  species = sample(c("a", "b", "c"), 60, replace = TRUE))
  d$trait <- rnorm(60)
  res <- tstats(d, nperm = 2000)
  null <- attr(res, "null")
  # each site draws without replacement from the pool, so its variance is unbiased for var(z)
  expect_lt(abs(mean(null[, , "T_IC.IR"]) - 1), 0.02)
})

test_that("tstats_summary() warns when nperm is too small for the level", {
  d <- toy_individuals()
  expect_warning(tstats_summary(tstats(d, nperm = 19)), "nperm >= 40")
  expect_no_warning(tstats_summary(tstats(d, nperm = 40)))
})

test_that("trait_flex() guards the model against missing values, names and no intercept", {
  d <- toy_individuals()
  env <- data.frame(site = sprintf("s%02d", 1:6), moisture = c(1, 3, 2, NA, 5, 4))
  expect_warning(res <- trait_flex(d, formula = ~moisture, env = env), "s04")
  dec <- res$decomposition
  expect_equal(dec$explained + dec$residual, dec$total)
  kept <- res$sites[res$sites$site != "s04", ]
  expect_equal(dec$total[dec$component == "total"],
               sum((kept$specific - mean(kept$specific))^2))
  # an environmental variable called y does not clash with the response
  env$y <- c(1, 3, 2, 6, 5, 4)
  env$moisture <- env$y
  expect_equal(suppressWarnings(trait_flex(d, formula = ~y, env = env))$decomposition,
               suppressWarnings(trait_flex(d, formula = ~moisture, env = env))$decomposition)
  expect_error(trait_flex(d, formula = ~ y - 1, env = env), "intercept")
  d$abundance <- seq_len(nrow(d))
  expect_error(trait_flex(d, abundance = "abundance"), "same for all")
})

test_that("nested_variance() rejects a factor with one level", {
  skip_if_not_installed("nlme")
  d <- toy_individuals()
  d$region <- "all"
  expect_error(nested_variance(d, "trait", c("region", "site")), "single level")
})

test_that("the functions handle untestable statistics, missing values and bad input", {
  # one species per site: T_PC.PR cannot be computed anywhere
  d <- data.frame(site = rep(c("s1", "s2", "s3"), each = 10),
                  species = rep(c("a", "b", "c"), each = 10))
  set.seed(4)
  d$trait <- rnorm(30)
  s <- suppressWarnings(tstats_summary(tstats(d, nperm = 39)))
  pc <- s$traits[s$traits$statistic == "T_PC.PR", ]
  expect_equal(pc$n_sites, 0)
  expect_true(is.na(pc$percent_low) && is.na(pc$p_mean_low))
  d <- toy_individuals()
  env <- data.frame(site = sprintf("s%02d", 1:6), moisture = 1:6)
  expect_error(trait_flex(d, formula = ~ depth, env = env), "no column depth")
  d$abundance <- NA_real_
  expect_error(trait_flex(d, abundance = "abundance"), "missing values")
  d$label <- "x"
  expect_error(trait_flex(d, traits = "label"), "numeric")
  expect_error(community_index(d, mean, traits = "label"), "numeric")
  d$trait[d$site == "s01"] <- NA
  res <- community_index(d, mean, nperm = 0)
  expect_true(is.na(res$observed[res$site == "s01"]))
})
