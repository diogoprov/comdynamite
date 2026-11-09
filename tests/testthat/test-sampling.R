toy_sampling <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(site = rep(c("s1", "s2"), each = 50),
                  species = sample(c("a", "b", "c"), 100, TRUE, prob = c(6, 3, 1)))
  d$trait <- c(a = 10, b = 12, c = 15)[d$species] + rnorm(100)
  d
}

test_that("subsample_individuals() applies n, fraction, cover and min_n", {
  d <- toy_sampling()
  sizes <- table(paste(d$site, d$species))
  sub <- subsample_individuals(d, n = 5)
  expect_equal(as.vector(table(paste(sub$site, sub$species))[names(sizes)]),
               as.vector(pmin(sizes, 5)))
  expect_true(all(rownames(sub) %in% rownames(d)))
  expect_false(is.unsorted(as.integer(rownames(sub))))
  # cover: species from the most abundant down until they reach the share
  covered <- subsample_individuals(d, cover = 0.6)
  for (s in c("s1", "s2")) {
    counts <- sort(table(d$species[d$site == s]), decreasing = TRUE)
    share <- cumsum(counts) / sum(counts)
    expected <- names(counts)[c(0, share[-length(share)]) < 0.6]
    expect_setequal(unique(covered$species[covered$site == s]), expected)
  }
  expect_equal(nrow(subsample_individuals(d, cover = 1)), nrow(d))
  # min_n removes small populations
  big <- subsample_individuals(d, min_n = 15)
  expect_true(all(table(paste(big$site, big$species)) >= 15))
  set.seed(2)
  half <- subsample_individuals(d, fraction = 0.5)
  expect_lt(abs(nrow(half) - 50), 20)
  # abundance column instead of counts
  d$abundance <- c(a = 1, b = 100, c = 1)[d$species]
  expect_equal(unique(subsample_individuals(d, cover = 0.5, abundance = "abundance")$species),
               "b")
  d$abundance[1] <- 999
  expect_error(subsample_individuals(d, cover = 0.5, abundance = "abundance"), "one value")
  expect_error(subsample_individuals(d, site = "nope"), "no column")
})

test_that("sample_individuals() measures only the dominant species with cover", {
  skip_if_not_installed("treeshapesim")
  set.seed(4)
  tree <- treeshapesim::sim_tree_shape(10)
  sim <- simulate_metacommunity(landscape_implicit(5), species_traits(tree),
                                competition_matrix(tree, "none"), n_steps = 50,
                                seed_steps = 1:5)
  all_species <- sample_individuals(sim, n = 3)
  dominant <- sample_individuals(sim, n = 3, cover = 0.5)
  pops <- unique(all_species[c("site", "species", "abundance")])
  for (s in unique(pops$site)) {
    p <- pops[pops$site == s, ]
    p <- p[order(-p$abundance), ]
    share <- cumsum(p$abundance) / sum(p$abundance)
    expected <- p$species[c(0, share[-length(share)]) < 0.5]
    expect_setequal(unique(dominant$species[dominant$site == s]), expected)
  }
  few <- sample_individuals(sim, n = 3, min_n = 3)
  expect_true(all(table(paste(few$site, few$species)) == 3))
})

test_that("trait_rarefaction() summarises subsamples and plots them", {
  d <- toy_sampling()
  set.seed(5)
  rare <- trait_rarefaction(d, n = c(2, 5), reps = 4)
  expect_s3_class(rare, "trait_rarefaction")
  expect_equal(rare$statistic, rep(c("T_IP.IC", "T_IC.IR", "T_PC.PR"), 2))
  values <- attr(rare, "values")
  expect_equal(rare$mean[1:3], unname(colMeans(values[["2"]])))
  expect_true(all(rare$low <= rare$mean & rare$mean <= rare$high))
  own <- trait_rarefaction(d, function(x) c(mean = mean(x$trait), n = nrow(x)), n = 3,
                           reps = 2)
  expect_equal(own$mean[own$statistic == "n"],
               sum(pmin(table(paste(d$site, d$species)), 3)))
  expect_error(trait_rarefaction(d, function(x) mean(x$trait), n = 2, reps = 1), "named")
  default_n <- trait_rarefaction(d, function(x) c(n = nrow(x)), reps = 1)
  expect_equal(min(default_n$effort), 2)
  expect_equal(attr(default_n, "by"), "n")
  # over cover: deterministic without n, one subsample per level
  by_cover <- trait_rarefaction(d, function(x) c(species = length(unique(x$species))),
                                cover = c(0.5, 1))
  expect_equal(attr(by_cover, "by"), "cover")
  expect_equal(by_cover$reps, c(1, 1))
  expect_equal(by_cover$mean[2], 3)
  expect_lt(by_cover$mean[1], 3)
  skip_if_not_installed("ggplot2")
  p <- plot(rare)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(ggplot2::layer_data(p, 2)), 6)
  expect_length(p$layers, 3)
  # the published rule for cover is drawn as a dashed line
  pc <- plot(by_cover)
  rule <- ggplot2::layer_data(pc, 1)
  expect_equal(unique(rule$xintercept), 0.8)
  expect_equal(unique(rule$linetype), "dashed")
  expect_match(unique(ggplot2::layer_data(pc, 2)$label), "Pakeman")
  expect_length(plot(by_cover, reference = FALSE)$layers, 3)
  own_rule <- plot(rare, reference = c(protocol = 4))
  expect_equal(unique(ggplot2::layer_data(own_rule, 1)$xintercept), 4)
})

test_that("sampling helpers handle edge inputs", {
  ind <- data.frame(site = rep(c("s1", "s2"), each = 6), species = rep(c("a", "b"), 6),
                    trait = 1:12)
  ind$ab <- ifelse(ind$site == "s1", 0, ifelse(ind$species == "a", 5, 1))
  sub <- subsample_individuals(ind, cover = 0.8, abundance = "ab")
  expect_true(all(sub$site == "s2"))
  expect_error(trait_rarefaction(ind, cover = 1.5))
  r <- trait_rarefaction(ind, n = 2:3, reps = 2)
  attr(r, "by") <- NULL
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(r), "ggplot")
})

test_that("trait_rarefaction over cover with a fixed n keeps the replicates", {
  set.seed(6)
  ind <- data.frame(site = rep(c("s1", "s2"), each = 20),
                    species = sample(c("a", "b", "c"), 40, replace = TRUE), trait = rnorm(40))
  r <- trait_rarefaction(ind, statistic = function(x) c(m = mean(x$trait)), cover = c(0.5, 1),
                         n = 2, reps = 3)
  expect_equal(attr(r, "by"), "cover")
  expect_equal(nrow(attr(r, "values")[[1]]), 3)
  expect_equal(r$reps, c(3, 3))
  expect_error(trait_rarefaction(ind, cover = 0.5, n = 1.5), "round")
})

test_that("the default statistic of trait_rarefaction runs no permutations", {
  set.seed(9)
  ind <- data.frame(site = rep(c("s1", "s2"), each = 20),
                    species = sample(c("a", "b", "c"), 40, replace = TRUE), trait = rnorm(40))
  calls <- 0
  trace_fun <- function(...) calls <<- calls + 1
  local_mocked_bindings(
    tstats = function(..., nperm) {
      expect_equal(nperm, 0)
      calls <<- calls + 1
      data.frame(T_IP.IC = 1, T_IC.IR = 1, T_PC.PR = 1)
    }
  )
  trait_rarefaction(ind, n = 2, reps = 2)
  expect_equal(calls, 2)
})
