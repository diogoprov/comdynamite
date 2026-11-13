# Evolution of several traits at once, one per environmental axis of the niche, with genetic
# and environmental covariance matrices. The univariate engine in mc_step() is left as it was;
# these helpers do the same steps on arrays of sites x species x traits.

# is this a setting of evolution for several traits?
evolution_is_mv <- function(evolution) {
  is.list(evolution$traits) && !is.data.frame(evolution$traits)
}

# a covariance matrix for each species from one number (times the identity), one matrix, or a
# list with one matrix per species
evolution_matrices <- function(x, q, species, name) {
  as_one <- function(v) {
    v <- if (length(v) == 1 && is.null(dim(v))) diag(v, q) else as.matrix(v)
    if (!is.numeric(v) || !identical(dim(v), c(q, q)) || anyNA(v) ||
        max(abs(v - t(v))) > 1e-8 * max(1, abs(v))) {
      stop("`", name, "` needs symmetric ", q, " x ", q, " matrices (one trait per axis)",
           call. = FALSE)
    }
    if (min(eigen(v, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 * max(1, abs(v))) {
      stop("`", name, "` must be positive semidefinite", call. = FALSE)
    }
    dimnames(v) <- NULL
    v
  }
  if (is.list(x)) {
    if (!is.null(names(x))) {
      if (!all(species %in% names(x))) {
        stop("`", name, "` lacks matrices for ", species_list(setdiff(species, names(x))),
             call. = FALSE)
      }
      x <- x[species]
    }
    if (length(x) != length(species)) {
      stop("`", name, "` as a list needs one matrix per species (", length(species), ")",
           call. = FALSE)
    }
    return(lapply(x, as_one))
  }
  rep(list(as_one(x)), length(species))
}

# a matrix W with W' P W = I on the range of P, so that distances between rows of z %*% W are
# Mahalanobis distances under P (directions without variance are left out)
evolution_whiten <- function(p) {
  e <- eigen(p, symmetric = TRUE)
  keep <- e$values > 1e-10 * max(1, e$values[1])
  e$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(e$values[keep]), sum(keep))
}

# Mahalanobis distances under the P of species j between the rows of z (populations x traits)
mv_distance <- function(ev, j, z) {
  as.matrix(stats::dist(z %*% ev$whiten[[j]]))
}

# any matrix with L %*% t(L) = v, also for singular v
evolution_root <- function(v) {
  e <- eigen(v, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(v))
}

# the state of multivariate evolution inside mc_dynamics()
evolution_state_mv <- function(evolution, m, k_all, species) {
  traits <- lapply(evolution$traits, as_traits)
  q <- length(traits)
  axes <- names(evolution$traits) %||% paste0("axis", seq_len(q))
  if (anyDuplicated(axes) || any(axes %in% c("site", "species", "individual", "abundance",
                                              "step", ""))) {
    stop("the axes of `env` need distinct names other than `site`, `species`, ",
         "`individual`, `abundance` and `step`", call. = FALSE)
  }
  if (!is.null(species) && all(vapply(traits, function(t) all(species %in% rownames(t)),
                                      logical(1)))) {
    traits <- lapply(traits, function(t) t[species, , drop = FALSE])
  }
  if (any(vapply(traits, nrow, integer(1)) != k_all)) {
    stop("`traits` of `evolution` need one row per species (", k_all, ") on every axis",
         call. = FALSE)
  }
  env_fun <- if (is.function(evolution$env)) evolution$env else {
    fixed_env <- evolution$env
    function(step) fixed_env
  }
  env <- as.matrix(env_fun(1))
  if (!is.numeric(env) || nrow(env) != m || ncol(env) != q || anyNA(env)) {
    stop("`env` of `evolution` needs one row per site (", m, ") and one column per trait (", q,
         ")", call. = FALSE)
  }
  dimnames(env) <- NULL
  optimum <- vapply(traits, function(t) t$optimum, numeric(k_all))
  sigma <- vapply(traits, function(t) t$sigma, numeric(k_all))
  optimum <- matrix(optimum, k_all, q)
  sigma <- matrix(sigma, k_all, q)
  r_axis <- matrix(vapply(traits, function(t) t$r_max, numeric(k_all)), k_all, q)
  r_max <- Reduce(`*`, lapply(traits, function(t) t$r_max))
  sp <- rownames(traits[[1]])
  G <- evolution_matrices(evolution$genetic_var, q, sp, "genetic_var")
  E <- evolution_matrices(evolution$env_var, q, sp, "env_var")
  b <- evolution$plasticity
  if (is.matrix(b) && !is.null(rownames(b)) && all(sp %in% rownames(b))) {
    b <- b[sp, , drop = FALSE]
  }
  b <- if (is.matrix(b)) b else if (length(b) == 1) matrix(b, k_all, q) else
    if (length(b) == q) matrix(b, k_all, q, byrow = TRUE) else NULL
  if (!is.numeric(b) || !identical(dim(b), c(k_all, q)) || anyNA(b)) {
    stop("`plasticity` needs one value, one per trait (", q, ") or a matrix of species x traits",
         call. = FALSE)
  }
  dimnames(b) <- NULL
  # selection on the axes with a finite niche breadth: the response of each population is
  # K (E - p), K = G[, S] (Omega + P)[S, S]^-1, and growth falls with the Mahalanobis distance
  sel <- lapply(seq_len(k_all), function(j) {
    om <- sigma[j, ]^2
    s <- which(is.finite(om))
    p <- G[[j]] + E[[j]]
    if (!length(s)) {
      return(list(s = s, A = matrix(0, 0, 0), K = matrix(0, q, 0), c = r_max[j]))
    }
    Q <- diag(om[s], length(s)) + p[s, s, drop = FALSE]
    if (any(diag(Q) <= 0)) {
      stop("a niche breadth `sigma` of 0 needs phenotypic variance on that axis", call. = FALSE)
    }
    A <- solve(Q)
    # without covariances among the selected axes, growth is computed axis by axis as in
    # growth_rates(), so that G = E = 0 gives the same numbers to the last bit
    diagonal <- all(Q[row(Q) != col(Q)] == 0)
    logc <- 0.5 * (sum(log(om[s])) - as.numeric(determinant(Q)$modulus))
    list(s = s, A = A, K = G[[j]][, s, drop = FALSE] %*% A, c = r_max[j] * exp(logc),
         diagonal = diagonal, width = om + diag(p), om = om, r_axis = r_axis[j, ])
  })
  env_now <- if (is.function(evolution$env)) {
    function(step) {
      e <- as.matrix(env_fun(step))
      dimnames(e) <- NULL
      e
    }
  } else {
    function(step) env
  }
  list(mv = TRUE, q = q, axes = axes, env = env_now, ref = colMeans(env),
       optimum = optimum, sigma = sigma, r_max = r_max, r_axis = r_axis, G = G, E = E, b = b,
       sel = sel, whiten = Map(function(g, e) evolution_whiten(g + e), G, E),
       L = lapply(G, evolution_root), fixed = vapply(G, function(g) all(g == 0), logical(1)),
       nu = evolution$ne_ratio, per = evolution$steps_per_generation,
       founder = evolution$founder)
}

# phenotypes: genetic means plus the plastic response to the environment e (sites x traits)
mv_phenotype <- function(ev, za, e, cols) {
  pa <- za
  for (a in seq_len(ev$q)) pa[, , a] <- za[, , a] + outer(e[, a] - ev$ref[a], ev$b[cols, a])
  pa
}

# growth rates (sites x species) of phenotypes pa under the multivariate Gaussian niche, with
# the variance of phenotypes around their mean
mv_growth <- function(ev, pa, e, cols) {
  r <- matrix(0, dim(pa)[1], length(cols))
  for (i in seq_along(cols)) {
    s <- ev$sel[[cols[i]]]
    if (!length(s$s)) {
      r[, i] <- s$c
      next
    }
    if (s$diagonal) {
      ri <- 1
      for (a in seq_len(ev$q)) {
        f <- if (is.finite(s$om[a])) {
          exp(-(e[, a] - pa[, i, a])^2 / (2 * s$width[a])) *
            (s$r_axis[a] * sqrt(s$om[a] / s$width[a]))
        } else {
          s$r_axis[a]
        }
        ri <- ri * f
      }
      r[, i] <- ri
      next
    }
    d <- e[, s$s, drop = FALSE] - matrix(pa[, i, s$s], dim(pa)[1])
    r[, i] <- s$c * exp(-0.5 * rowSums((d %*% s$A) * d))
  }
  r
}

# selection and drift in the populations present after reproduction (n > 0)
mv_select <- function(ev, za, pa, e, n, cols) {
  for (i in seq_along(cols)) {
    j <- cols[i]
    live <- n[, i] > 0
    if (!any(live) || ev$fixed[j]) next
    s <- ev$sel[[j]]
    z <- matrix(za[live, i, ], sum(live))
    if (length(s$s)) {
      d <- e[live, s$s, drop = FALSE] - matrix(pa[live, i, s$s], sum(live))
      z <- z + (d %*% t(s$K)) / ev$per
    }
    if (is.finite(ev$nu)) {
      noise <- matrix(stats::rnorm(sum(live) * ev$q), sum(live)) %*% t(ev$L[[j]])
      z <- z + noise * sqrt(1 / (ev$per * ev$nu * n[live, i]))
    }
    za[live, i, ] <- z
  }
  za[rep(n == 0, ev$q)] <- NA_real_
  za
}

# means of founders in `cells` (sites x species of the alive columns `cols`): k founders drawn
# from sources with means `from`, with variance `between` (per trait) among the sources
mv_found <- function(ev, target, from, k, cells, cols, between = NULL) {
  for (i in which(colSums(cells) > 0)) {
    rows <- which(cells[, i])
    base <- matrix(from[rows, i, ], length(rows))
    j <- cols[i]
    if (ev$founder && !ev$fixed[j]) {
      kk <- k[rows, i]
      draw <- matrix(stats::rnorm(length(rows) * ev$q), length(rows)) %*% t(ev$L[[j]])
      if (!is.null(between)) {
        extra <- matrix(between[rows, i, ], length(rows))
        draw <- draw + matrix(stats::rnorm(length(rows) * ev$q), length(rows)) * sqrt(extra)
      }
      base <- base + draw / sqrt(kk)
    }
    target[rows, i, ] <- base
  }
  target
}

# gene flow: residents and immigrants mixed by their numbers, immigrants with the means of
# their sources; populations without residents are founded by a sample of their immigrants
mv_flow <- function(ev, za, pa, n, em, imm, disp_matrix, cols) {
  imm <- imm + matrix(0L, nrow(n), ncol(n))
  res <- n - em
  zi <- between <- p_imm <- za
  for (a in seq_len(ev$q)) {
    zi[, , a] <- immigrant_mean(em, za[, , a], disp_matrix, za[, , a])
    # variance among sources; the residue of the subtraction (cancellation for large means)
    # is set to 0
    v <- immigrant_mean(em, za[, , a]^2, disp_matrix, za[, , a]^2) - zi[, , a]^2
    v[v < 1e-10 * pmax(1, zi[, , a]^2)] <- 0
    between[, , a] <- v
    p_imm[, , a] <- immigrant_mean(em, pa[, , a], disp_matrix, pa[, , a])
  }
  both <- res > 0 & imm > 0
  w_res <- ifelse(both, res / (res + imm), 1)
  for (a in seq_len(ev$q)) {
    za[, , a] <- ifelse(both, w_res * za[, , a] + (1 - w_res) * zi[, , a], za[, , a])
    pa[, , a] <- ifelse(both, w_res * pa[, , a] + (1 - w_res) * p_imm[, , a], pa[, , a])
  }
  founded <- res == 0 & imm > 0
  if (any(founded)) {
    za <- mv_found(ev, za, zi, imm, founded, cols, between)
    for (a in seq_len(ev$q)) pa[, , a][founded] <- p_imm[, , a][founded]
  }
  list(za = za, pa = pa)
}

#' Genetic and environmental covariance matrices for several traits
#'
#' Builds one additive genetic covariance matrix \eqn{G}{G} and one environmental covariance
#' matrix \eqn{E}{E} per species, for the evolution of several traits (one per niche axis) in
#' [evolution_control()], under three families of relationship between them:
#'
#' * `"proportional"`: \eqn{G = h^2 P}{G = h2 P} and \eqn{E = (1 - h^2) P}{E = (1 - h2) P},
#'   both proportional to a phenotypic matrix \eqn{P}{P} (`P`) that is shared by all species up
#'   to a scale factor. Under neutral evolution with such matrices, the covariance among
#'   population means is expected to be proportional to that within populations (Lande 1979),
#'   the expectation tested with phenotypic matrices by Ackermann & Cheverud (2004) and whose
#'   error rates were studied by Prôa et al. (2013).
#' * `"independent"`: \eqn{G}{G} and \eqn{E}{E} have no common structure: each species draws
#'   two random orientations (eigenvectors) and two sets of variances, with total variances
#'   \eqn{h^2}{h2} and \eqn{1 - h^2}{1 - h2} times that of `P`.
#' * `"shared"`: \eqn{G}{G} and \eqn{E}{E} share the principal components (eigenvectors) of `P`
#'   but differ in the variance along them, which each species draws around \eqn{h^2}{h2} and
#'   \eqn{1 - h^2}{1 - h2} times the eigenvalues of `P`.
#'
#' `divergence` sets how much species differ: the standard deviation of the log of the scale
#' factor of `P` (`"proportional"`) or of the eigenvalues (`"independent"` and `"shared"`).
#'
#' @param species Names of the species, or their number.
#' @param n_traits Number of traits (axes). Ignored when `P` is given.
#' @param family One of `"proportional"`, `"independent"` and `"shared"`.
#' @param P Phenotypic covariance matrix, symmetric and positive definite. By default, unit
#'   variances with a common correlation `correlation`.
#' @param correlation Correlation between all pairs of traits in the default `P`.
#' @param heritability Share \eqn{h^2}{h2} of the phenotypic variance that is additive genetic,
#'   between 0 and 1.
#' @param divergence Standard deviation (on the log scale) of the differences among species.
#' @return A list with `genetic_var` and `env_var`, each a named list of one matrix per species,
#'   ready for `evolution_control(genetic_var = m$genetic_var, env_var = m$env_var)`.
#' @references
#'   Ackermann, R. R., & Cheverud, J. M. (2004). Detecting genetic drift versus selection in
#'   human evolution. *Proceedings of the National Academy of Sciences*, 101, 17946--17951.
#'   \doi{10.1073/pnas.0405919102}
#'
#'   Lande, R. (1979). Quantitative genetic analysis of multivariate evolution, applied to
#'   brain:body size allometry. *Evolution*, 33, 402--416.
#'   \doi{10.1111/j.1558-5646.1979.tb04694.x}
#'
#'   Prôa, M., O'Higgins, P., & Monteiro, L. R. (2013). Type I error rates for testing genetic
#'   drift with phenotypic covariance matrices: a simulation study. *Evolution*, 67, 185--195.
#'   \doi{10.1111/j.1558-5646.2012.01746.x}
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' m <- genetic_matrices(c("sp1", "sp2"), n_traits = 2, family = "shared",
#'                       correlation = 0.6, heritability = 0.4)
#' m$genetic_var$sp1
#' ev <- evolution_control(genetic_var = m$genetic_var, env_var = m$env_var)
#' @export
genetic_matrices <- function(species, n_traits = 2,
                             family = c("proportional", "independent", "shared"), P = NULL,
                             correlation = 0, heritability = 0.5, divergence = 0) {
  family <- match.arg(family)
  if (is.numeric(species) && length(species) == 1) species <- paste0("sp", seq_len(species))
  species <- as.character(species)
  stopifnot(length(species) >= 1, !anyDuplicated(species),
            is.numeric(heritability), length(heritability) == 1, heritability >= 0,
            heritability <= 1, is.numeric(divergence), length(divergence) == 1,
            divergence >= 0)
  if (is.null(P)) {
    stopifnot(is.numeric(n_traits), length(n_traits) == 1, n_traits >= 1,
              n_traits == round(n_traits), is.numeric(correlation), length(correlation) == 1,
              correlation > -1 / max(n_traits - 1, 1), correlation < 1)
    P <- matrix(correlation, n_traits, n_traits)
    diag(P) <- 1
  }
  labels <- dimnames(as.matrix(P))
  P <- evolution_matrices(P, nrow(as.matrix(P)), "P", "P")[[1]]
  q <- nrow(P)
  eig <- eigen(P, symmetric = TRUE)
  if (min(eig$values) <= 0) stop("`P` must be positive definite", call. = FALSE)
  build <- function(vectors, values) {
    m <- vectors %*% diag(values, q) %*% t(vectors)
    (m + t(m)) / 2
  }
  rotation <- function() qr.Q(qr(matrix(stats::rnorm(q * q), q)))
  spread <- function() exp(stats::rnorm(q, 0, divergence))
  one <- function() {
    switch(family,
      proportional = {
        s <- exp(stats::rnorm(1, 0, divergence))
        list(G = heritability * s * P, E = (1 - heritability) * s * P)
      },
      independent = {
        total <- sum(eig$values)
        vg <- eig$values[sample.int(q)] * spread()
        ve <- eig$values[sample.int(q)] * spread()
        list(G = build(rotation(), heritability * total * vg / sum(vg)),
             E = build(rotation(), (1 - heritability) * total * ve / sum(ve)))
      },
      shared = list(G = build(eig$vectors, heritability * eig$values * spread()),
                    E = build(eig$vectors, (1 - heritability) * eig$values * spread()))
    )
  }
  mats <- lapply(species, function(s) one())
  name <- function(m) {
    dimnames(m) <- labels
    m
  }
  list(genetic_var = stats::setNames(lapply(mats, function(m) name(m$G)), species),
       env_var = stats::setNames(lapply(mats, function(m) name(m$E)), species))
}
