# A genetic architecture of a few diploid, biallelic, additive loci (`loci` of
# evolution_control()). Each population carries the frequency of allele 1 at each locus, as an
# array of sites x species x loci, and its genetic mean and variance follow from them under
# Hardy-Weinberg proportions and linkage equilibrium. The helpers at the top are functions of
# the frequencies; those at the bottom (`loci_step_*`) change the state frame of mc_dynamics().
# Selection changes the frequencies by the formula of Wright (1937), and the loss of variance
# to stabilising selection at loci with two alleles is that of Barton (1986); both are in the
# references of evolution_control().

#' stops unless `loci` is a list with `n`, `effect`, `mutation` and, optionally, `origin`;
#' returns it with `origin` filled in
#' @noRd
check_loci <- function(loci) {
  known <- c("n", "effect", "mutation", "origin")
  ok <- is.list(loci) && !is.null(names(loci)) && all(names(loci) %in% known) &&
    all(known[1:3] %in% names(loci)) && !anyDuplicated(names(loci))
  if (!ok) {
    stop("`loci` must be `NULL` or a list with `n`, `effect` and `mutation` (and, optionally, ",
         "`origin`)", call. = FALSE)
  }
  check_count(loci$n, "loci$n")
  check_number(loci$effect, "loci$effect", lower = 0, open = c(TRUE, FALSE))
  check_number(loci$mutation, "loci$mutation", lower = 0, upper = 0.5)
  if (is.null(loci$origin)) loci$origin <- 0
  check_number(loci$origin, "loci$origin")
  loci[known]
}

#' stops when `loci` meets a setting of evolution_control() that it cannot go with, and warns
#' when `genetic_var` was given (`given`), since the variance comes from the loci
#' @noRd
check_loci_settings <- function(given, env_var, mutation_var, variance_selection,
                                immigrant_selection) {
  if (!is.null(mutation_var)) {
    stop("`loci` cannot go with `mutation_var`: the variance comes from the allele ",
         "frequencies, and mutation from `loci$mutation`", call. = FALSE)
  }
  if (!isFALSE(variance_selection)) {
    stop("`loci` cannot go with `variance_selection`: selection changes the allele ",
         "frequencies, and the variance with them", call. = FALSE)
  }
  if (!isFALSE(immigrant_selection)) {
    stop("`loci` cannot go with `immigrant_selection = TRUE`", call. = FALSE)
  }
  if (is_covariance(env_var)) {
    stop("`loci` are for one trait: `env_var` cannot be a covariance matrix", call. = FALSE)
  }
  if (given) {
    warning("`genetic_var` is ignored with `loci`: the genetic variance comes from the allele ",
            "frequencies", call. = FALSE)
  }
  invisible(NULL)
}

#' the allele frequencies of the regional pool of each species (species x loci), which give
#' its `optimum`: the first loci fixed for allele 1, at most one locus at an intermediate
#' frequency and the others fixed for allele 0
#' @noRd
loci_initial <- function(optimum, loci) {
  top <- loci$origin + 2 * loci$effect * loci$n
  full <- (optimum - loci$origin) / (2 * loci$effect)
  if (any(full < -1e-8 | full > loci$n + 1e-8)) {
    stop("with `loci`, the optimum of every species must be between `origin` (", loci$origin,
         ") and `origin + 2 * effect * n` (", top, ")", call. = FALSE)
  }
  full[abs(full - round(full)) < 1e-8] <- round(full[abs(full - round(full)) < 1e-8])
  p <- outer(full, seq_len(loci$n) - 1, "-")
  matrix(pmin(pmax(p, 0), 1), length(optimum), loci$n)
}

#' the genetic mean of each population (sites x species) from its allele frequencies `p`
#' (sites x species x loci): origin + 2 * effect * sum of the frequencies; NA where absent
#' @noRd
loci_mean <- function(p, loci) {
  loci$origin + 2 * loci$effect * rowSums(p, dims = 2)
}

#' the genetic variance of each population (sites x species) from its allele frequencies `p`,
#' with Hardy-Weinberg proportions and linkage equilibrium: 2 * effect^2 * sum of p (1 - p)
#' @noRd
loci_variance <- function(p, loci) {
  2 * loci$effect^2 * rowSums(p * (1 - p), dims = 2)
}

#' the allele frequencies `p` after Gaussian selection in the populations `cells`, by the
#' formula of Wright (1937), dp = p (1 - p) / 2 * dlog(W) / dp, with the mean fitness
#' log(W) = -log(V) / 2 - d^2 / (2 V), where `d` (sites x species) is the distance of the mean
#' phenotype from the optimum and `width` (V = sigma^2 + G + VE) the width of selection:
#' dp = p (1 - p) (a d / V - a^2 (1 - 2 p) / (2 V) (1 - d^2 / V)), over `per` steps per
#' generation. The first term changes the mean (summed over loci, the response of Lande 1976)
#' and the second, which comes from G in V, removes variance near the optimum (Barton 1986)
#' @noRd
loci_select <- function(p, cells, d, width, per, loci) {
  n_loci <- dim(p)[3]
  a <- loci$effect
  idx <- rep(cells, n_loci)
  mean_term <- a * rep(d / width, n_loci)
  var_term <- a^2 * (1 - 2 * p) / 2 * rep((1 - d^2 / width) / width, n_loci)
  out <- p + p * (1 - p) * (mean_term - var_term) / per
  p[idx] <- pmin(pmax(out[idx], 0), 1)
  p
}

#' the allele frequencies `p` after a fraction `f` of a generation (a number, or sites x
#' species) of mutation between the two alleles, each at the rate `mu` per generation
#' @noRd
loci_mutate <- function(p, mu, f) {
  p + rep(f * mu, length.out = length(p)) * (1 - 2 * p)
}

#' drift of the allele frequencies `p` in the populations `cells`: a binomial sample of
#' `size` (sites x species, rounded, at least 1) alleles at each locus that is not fixed
#' @noRd
loci_drift <- function(p, cells, size) {
  n_loci <- dim(p)[3]
  idx <- which(rep(cells, n_loci) & !is.na(p) & p > 0 & p < 1)
  if (!length(idx)) return(p)
  draws <- pmax(round(rep(size, n_loci)[idx]), 1)
  p[idx] <- stats::rbinom(length(idx), draws, p[idx]) / draws
  p
}

#' the allele frequencies of two groups merged in the populations `cells`: `w` (sites x
#' species) of `p1` and `1 - w` of `p2`; a group without frequencies (NA) counts as absent
#' @noRd
loci_mix <- function(p1, p2, w, cells) {
  idx <- rep(cells, dim(p1)[3])
  a <- ifelse(is.na(p1), 0, p1)
  b <- ifelse(is.na(p2), 0, p2)
  out <- rep(w, dim(p1)[3]) * a + rep(1 - w, dim(p1)[3]) * b
  p1[idx] <- out[idx]
  p1
}

#' the allele frequencies of the immigrants into each population: the mean of those of their
#' sources, weighted by the emigrants `em` that each sends, locus by locus
#' @noRd
loci_immigrants <- function(em, p, disp_matrix) {
  for (l in seq_len(dim(p)[3])) {
    x <- matrix(p[, , l], nrow(em), ncol(em))
    p[, , l] <- immigrant_mean(em, x, disp_matrix, x)
  }
  p
}

#' the allele frequencies at the start, in the state `st`: every population present carries
#' those of the pool of its species (`pool_p`), and no species is frozen
#' @noRd
loci_setup <- function(st) {
  ev <- st$ev
  dims <- c(st$m, st$k_all, ev$loci$n)
  st$pool_p <- array(rep(ev$loci_pool, each = st$m), dims)
  st$lp <- array(NA_real_, dims)
  at_start <- rep(st$start_pop, ev$loci$n)
  st$lp[at_start] <- st$pool_p[at_start]
  st$pool_mean <- loci_mean(st$pool_p, ev$loci)
  st$pool_g <- loci_variance(st$pool_p, ev$loci)
  st$z <- loci_mean(st$lp, ev$loci)
  st$vg <- loci_variance(st$lp, ev$loci)
  st$fixed_z <- rep(FALSE, st$k_all)
  invisible(st)
}

#' the allele frequencies after the seeds of this step, which carry those of the pool and mix
#' with residents by their numbers; changes `lp`, `z` and `vg` in the state `st`
#' @noRd
loci_step_seed <- function(st) {
  alive <- st$alive
  got <- st$added_m > 0
  la <- loci_mix(st$lp[, alive, , drop = FALSE], st$pool_p[, alive, , drop = FALSE],
                 st$n / pmax(st$n + st$added_m, 1), got)
  st$lp[, alive, ] <- la
  st$z[, alive] <- loci_mean(la, st$ev$loci)
  st$vg[, alive] <- loci_variance(la, st$ev$loci)
  invisible(st)
}

#' the allele frequencies of the alive species at the start of a step (`la`, those of the
#' pool where a population is absent) and the genetic means that follow (`za`), in `st`
#' @noRd
loci_step_start <- function(st) {
  la <- st$lp[, st$alive, , drop = FALSE]
  pool <- st$pool_p[, st$alive, , drop = FALSE]
  la[is.na(la)] <- pool[is.na(la)]
  st$la <- la
  st$za <- loci_mean(la, st$ev$loci)
  invisible(st)
}

#' writes the allele frequencies `la` of the alive species to the state `st`, with the genetic
#' means (`za`), variances (`ga`, `vg`) and phenotypes (`pa`, unless given) that follow
#' @noRd
loci_update <- function(st, la, pa = NULL) {
  ev <- st$ev
  st$la <- la
  st$za <- loci_mean(la, ev$loci)
  st$ga <- loci_variance(la, ev$loci)
  st$vg[, st$alive] <- st$ga
  st$pa <- if (!is.null(pa)) pa else {
    st$za + sweep(matrix(st$e_now - ev$ref, st$m, st$k), 2, ev$b[st$alive], "*")
  }
  invisible(st)
}

#' the adults that survive a step of the lottery are a sample of the adults: drift of their
#' allele frequencies (`l_surv`), with their mean and variance (`z_surv`, `g_surv`), in `st`
#' @noRd
loci_step_survive <- function(st) {
  ev <- st$ev
  cells <- st$surv > 0 & st$x > 0 & is.finite(ev$nu)
  st$l_surv <- loci_drift(st$la, cells, 2 * ev$nu / ifelse(cells, st$x, 1))
  st$z_surv <- loci_mean(st$l_surv, ev$loci)
  st$g_surv <- loci_variance(st$l_surv, ev$loci)
  invisible(st)
}

#' selection, mutation and drift of the allele frequencies in the populations present after
#' reproduction; with a lottery, mutation and drift act on the recruits instead
#' @noRd
loci_step_select <- function(st) {
  ev <- st$ev
  width <- st$ga + matrix(ev$sigma[st$alive]^2 + ev$ve[st$alive], st$m, st$k, byrow = TRUE)
  live <- st$n > 0
  la <- loci_select(st$la, live, st$e_now - st$pa, width, ev$per, ev$loci)
  if (!st$lottery) {
    la <- loci_mutate(la, ev$loci$mutation, 1 / ev$per)
    if (is.finite(ev$nu)) la <- loci_drift(la, live, 2 * ev$nu * st$n * ev$per)
  }
  la[rep(!live, ev$loci$n)] <- NA_real_
  loci_update(st, la)
}

#' gene flow: residents and immigrants mix their allele frequencies by their numbers, and
#' populations founded by immigrants take a binomial sample of theirs (at recruitment, with a
#' lottery); changes the evolutionary state in `st`
#' @noRd
loci_step_flow <- function(st) {
  ev <- st$ev
  res <- st$res
  imm <- st$imm_m
  both <- res > 0 & imm > 0
  founded <- res == 0 & imm > 0
  li <- loci_immigrants(st$em_to, st$la, st$disp_matrix)
  la <- loci_mix(st$la, li, res / pmax(res + imm, 1), both | founded)
  pa_now <- st$pa
  pa_now[both] <- ((res * st$pa + imm * st$p_imm) / (res + imm))[both]
  pa_now[founded] <- st$p_imm[founded]
  if (st$lottery) {
    st$no_drift[founded] <- !ev$founder
  } else if (ev$founder) {
    la <- loci_drift(la, founded, 2 * imm)
  }
  loci_update(st, la, pa_now)
}

#' allele frequencies after a lottery: those of the survivors and of the recruits, which are a
#' sample of the offspring, weighted by their numbers, with mutation in the recruits
#' @noRd
loci_step_recruit <- function(st) {
  ev <- st$ev
  rec <- st$rec
  la <- st$la
  if (is.finite(ev$nu)) la <- loci_drift(la, rec > 0 & !st$no_drift, 2 * ev$nu * rec)
  has <- st$n_new > 0
  la <- loci_mix(st$l_surv, la, st$surv / pmax(st$n_new, 1), has)
  la <- loci_mutate(la, ev$loci$mutation, rec / pmax(st$n_new, 1))
  la[rep(!has, ev$loci$n)] <- NA_real_
  loci_update(st, la)
}

#' writes the allele frequencies of the alive species back to `lp` (NA where the population is
#' gone), in the state `st`
#' @noRd
loci_step_store <- function(st) {
  la <- st$la
  la[rep(st$n == 0, st$ev$loci$n)] <- NA_real_
  st$lp[, st$alive, ] <- la
  invisible(st)
}

#' the allele frequencies of the output: sites x species x loci, NA where absent
#' @noRd
loci_output <- function(lp, absent, sites, species) {
  lp[rep(absent, dim(lp)[3])] <- NA_real_
  dimnames(lp) <- list(sites, species, paste0("locus", seq_len(dim(lp)[3])))
  lp
}
