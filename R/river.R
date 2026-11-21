#' Random dendritic river network
#'
#' Grows a branching river network site by site, from the outlet upwards, for
#' [landscape_river()]. Each new site joins the network upstream of an existing one: with
#' probability `branching` it starts a new branch at a site that is not a headwater (a
#' confluence), and otherwise it extends a headwater, so that low values give long, little
#' branched channels and high values bushy networks. It is a simple random growth model, not an
#' optimal channel network (Rinaldo et al. 2014); networks from other sources, such as the
#' OCNet package or a real river, can be passed to [landscape_river()] as a data frame with
#' the downstream neighbour and the length of the reach of each site.
#'
#' The coordinates place each site at its distance from the outlet along the channel, in a
#' direction inside a wedge of angles that each confluence splits among its branches in
#' proportion to their number of headwaters. Straight-line distances are therefore shorter
#' than distances along the channel, and sites on neighbouring branches can be close in a
#' straight line but far apart along the channel, as in real catchments.
#'
#' @param n_sites Number of sites (at least 2).
#' @param branching Probability that a new site starts a branch rather than extending a
#'   headwater (between 0 and 1).
#' @param reach_length Mean length of the reach between a site and its downstream neighbour.
#'   Lengths are drawn from an exponential distribution with this mean, or are all equal with
#'   `fixed_length = TRUE`.
#' @param fixed_length Logical: equal reach lengths.
#' @param max_branches Largest number of branches that join at a confluence (2 for binary
#'   confluences).
#' @return A data frame with one row per site, the outlet first: `site`, `downstream` (row of
#'   the downstream neighbour, `NA` at the outlet), `length` (of the reach to the downstream
#'   neighbour, 0 at the outlet), `distance` (along the channel to the outlet), `magnitude`
#'   (Shreve magnitude, the number of headwaters upstream of the site, itself included),
#'   `order` (Strahler order), `headwater` (logical) and the coordinates `x` and `y`.
#' @references Rinaldo, A., Rigon, R., Banavar, J. R., Maritan, A., & Rodriguez-Iturbe, I.
#'   (2014). Evolution and selection of river networks: Statics, dynamics, and complexity.
#'   *Proceedings of the National Academy of Sciences*, 111, 2417--2424.
#'   \doi{10.1073/pnas.1322700111}
#' @family landscapes
#' @examples
#' set.seed(1)
#' net <- river_network(30)
#' head(net)
#' plot(net$x, net$y, asp = 1, pch = 19, cex = 0.6)
#' down <- !is.na(net$downstream)
#' segments(net$x[down], net$y[down], net$x[net$downstream[down]], net$y[net$downstream[down]])
#' @export
river_network <- function(n_sites, branching = 0.3, reach_length = 1, fixed_length = FALSE,
                          max_branches = 2) {
  stopifnot(length(n_sites) == 1, n_sites >= 2, length(branching) == 1, branching >= 0,
            branching <= 1, length(reach_length) == 1, reach_length > 0,
            is.logical(fixed_length), length(max_branches) == 1, max_branches >= 2)
  n_sites <- as.integer(n_sites)
  downstream <- c(NA_integer_, rep(0L, n_sites - 1))
  children <- integer(n_sites)
  for (i in 2:n_sites) {
    open <- seq_len(i - 1)
    tips <- open[children[open] == 0]
    forks <- open[children[open] > 0 & children[open] < max_branches]
    pool <- if (length(forks) && (!length(tips) || stats::runif(1) < branching)) forks else tips
    to <- pool[sample.int(length(pool), 1)]
    downstream[i] <- to
    children[to] <- children[to] + 1L
  }
  len <- if (fixed_length) rep(reach_length, n_sites) else stats::rexp(n_sites, 1 / reach_length)
  len[1] <- 0
  river_table(downstream, len)
}

#' Landscape on a river network
#'
#' Builds a landscape whose sites lie on a dendritic river network, with dispersal along the
#' channel that can differ upstream and downstream, as in the riverine metacommunities of
#' Carrara et al. (2012) and Altermatt (2013). The probability that an emigrant from site s
#' reaches site t is proportional to
#' \deqn{\exp\left(-\frac{d_{\downarrow}}{s_{\downarrow}} - \frac{d_{\uparrow}}{s_{\uparrow}}
#' \right),}{exp(-d_down / scale_down - d_up / scale_up),}
#' where \eqn{d_\downarrow}{d_down} is the distance it travels downstream (from s to the
#' confluence where the paths of s and t meet) and \eqn{d_\uparrow}{d_up} the distance it then
#' travels upstream (from the confluence to t). A `scale_up` smaller than `scale_down` gives
#' the drift of passive dispersers, and equal scales a symmetric kernel along the channel. A
#' share `overland` of emigrants instead moves in a straight line, with an exponential kernel
#' of the straight-line distance (scale `overland_scale`), for organisms with an aerial or
#' terrestrial stage that can cross between branches.
#'
#' @param network A data frame with one row per site and the columns `downstream` (row of the
#'   downstream neighbour, `NA` for the single outlet) and `length` (of the reach to the
#'   downstream neighbour), and optionally `x` and `y`; for instance from [river_network()].
#'   Without `x` and `y`, coordinates are laid out as in [river_network()].
#' @param env Numeric vector with the environment of each site, or `NULL` (the default) for
#'   the distance from the outlet along the channel rescaled to `env_range`, a gradient from
#'   the lowlands to the headwaters.
#' @param scale_down,scale_up Scales of the kernel for the distance travelled downstream and
#'   upstream, in the units of `length`.
#' @param overland Share of emigrants that move in a straight line rather than along the
#'   channel (between 0 and 1).
#' @param overland_scale Scale of the exponential kernel of straight-line dispersal.
#' @param env_range Range of the environment when `env` is `NULL`.
#' @return An `mc_landscape` object, as from [landscape_matrix()], whose `dist` holds the
#'   distances along the channel and which also keeps `network` (the network completed as in
#'   [river_network()]) and `downstream` (sites x sites matrix of the distance travelled
#'   downstream from the site in the row to the site in the column; its transpose is the
#'   distance travelled upstream).
#' @references Altermatt, F. (2013). Diversity in riverine metacommunities: a network
#'   perspective. *Aquatic Ecology*, 47, 365--377. \doi{10.1007/s10452-013-9450-3}
#'
#'   Carrara, F., Altermatt, F., Rodriguez-Iturbe, I., & Rinaldo, A. (2012). Dendritic
#'   connectivity controls biodiversity patterns in experimental metacommunities.
#'   *Proceedings of the National Academy of Sciences*, 109, 5761--5766.
#'   \doi{10.1073/pnas.1119651109}
#' @family landscapes
#' @examples
#' set.seed(1)
#' land <- landscape_river(river_network(40), scale_down = 3, scale_up = 1)
#' land
#' # emigrants from a headwater go mostly downstream
#' head(order(-land$dispersal[which.max(land$network$distance), ]))
#' traits <- data.frame(optimum = c(20, 50, 80), sigma = 20, r_max = 4,
#'                      row.names = paste0("sp", 1:3))
#' alpha <- matrix(0.05, 3, 3, dimnames = list(rownames(traits), rownames(traits)))
#' sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.1, n_steps = 100)
#' @export
landscape_river <- function(network = river_network(50), env = NULL, scale_down = 2,
                            scale_up = scale_down, overland = 0, overland_scale = 1,
                            env_range = c(0, 100)) {
  stopifnot(is.data.frame(network), all(c("downstream", "length") %in% names(network)),
            length(scale_down) == 1, scale_down > 0, length(scale_up) == 1, scale_up > 0,
            length(overland) == 1, overland >= 0, overland <= 1, length(overland_scale) == 1,
            overland_scale > 0, length(env_range) == 2, env_range[1] < env_range[2])
  xy <- if (all(c("x", "y") %in% names(network))) network[, c("x", "y")]
  net <- river_table(network$downstream, network$length)
  if (!is.null(xy)) {
    stopifnot(is.numeric(xy$x), is.numeric(xy$y), all(is.finite(c(xy$x, xy$y))))
    net$x <- xy$x
    net$y <- xy$y
  }
  m <- nrow(net)
  if (is.null(env)) {
    env <- net$distance / max(net$distance) * diff(env_range) + env_range[1]
  }
  stopifnot(is.numeric(env), length(env) == m, all(is.finite(env)))
  down <- river_downstream(net)
  w <- -down / scale_down - t(down) / scale_up
  diag(w) <- -Inf
  channel <- exp(w - apply(w, 1, max))
  channel <- channel / rowSums(channel)
  coords <- data.frame(x = net$x, y = net$y)
  dispersal <- if (overland > 0) {
    straight <- dispersal_kernel(as.matrix(stats::dist(coords)), overland_scale)
    (1 - overland) * channel + overland * straight
  } else channel
  out <- landscape_matrix(dispersal, env, coords = coords, dist = down + t(down))
  out$network <- net
  out$downstream <- down
  out$settings <- list(scale_down = scale_down, scale_up = scale_up, overland = overland,
                       overland_scale = overland_scale)
  out
}

# the network table of river_network() from the downstream neighbour and reach length of each
# site; checks that the sites form a single tree draining to one outlet
river_table <- function(downstream, length) {
  m <- base::length(downstream)
  downstream <- as.integer(downstream)
  length <- as.numeric(length)
  if (m < 2 || base::length(length) != m) {
    stop("the network needs at least 2 sites, each with `downstream` and `length`",
         call. = FALSE)
  }
  outlet <- which(is.na(downstream))
  if (base::length(outlet) != 1) {
    stop("the network needs exactly one outlet (one `NA` in `downstream`)", call. = FALSE)
  }
  if (any(downstream[-outlet] < 1 | downstream[-outlet] > m | downstream[-outlet] ==
          seq_len(m)[-outlet])) {
    stop("`downstream` must give the row of another site", call. = FALSE)
  }
  if (any(!is.finite(length[-outlet]) | length[-outlet] <= 0)) {
    stop("reach lengths must be positive", call. = FALSE)
  }
  length[outlet] <- 0
  # sites from the outlet upwards: every site after its downstream neighbour
  children <- split(seq_len(m)[-outlet], factor(downstream[-outlet], levels = seq_len(m)))
  ord <- outlet
  i <- 1
  while (i <= base::length(ord)) {
    ord <- c(ord, children[[ord[i]]])
    i <- i + 1
  }
  if (base::length(ord) != m) {
    stop("`downstream` must form a tree draining to the outlet, without loops", call. = FALSE)
  }
  distance <- numeric(m)
  for (s in ord[-1]) distance[s] <- distance[downstream[s]] + length[s]
  n_children <- lengths(children)
  magnitude <- as.numeric(n_children == 0)
  order <- ifelse(n_children == 0, 1L, 0L)
  for (s in rev(ord)) {
    ch <- children[[s]]
    if (!base::length(ch)) next
    magnitude[s] <- sum(magnitude[ch])
    top <- max(order[ch])
    order[s] <- top + as.integer(sum(order[ch] == top) > 1)
  }
  # angular layout: each site takes the middle of its wedge, split among branches by magnitude
  lo <- hi <- numeric(m)
  lo[outlet] <- pi / 6
  hi[outlet] <- 5 * pi / 6
  for (s in ord) {
    ch <- children[[s]]
    if (!base::length(ch)) next
    cuts <- lo[s] + (hi[s] - lo[s]) * c(0, cumsum(magnitude[ch])) / sum(magnitude[ch])
    lo[ch] <- cuts[-base::length(cuts)]
    hi[ch] <- cuts[-1]
  }
  angle <- (lo + hi) / 2
  data.frame(site = paste0("site", seq_len(m)), downstream = downstream, length = length,
             distance = distance, magnitude = magnitude, order = as.integer(order),
             headwater = n_children == 0, x = distance * cos(angle),
             y = distance * sin(angle), stringsAsFactors = FALSE)
}

# distance travelled downstream from each site (rows) to each other site (columns): from the
# site to the confluence where their paths meet
river_downstream <- function(net) {
  m <- nrow(net)
  outlet <- which(is.na(net$downstream))
  children <- split(seq_len(m)[-outlet], factor(net$downstream[-outlet], levels = seq_len(m)))
  ord <- outlet
  i <- 1
  while (i <= length(ord)) {
    ord <- c(ord, children[[ord[i]]])
    i <- i + 1
  }
  below <- vector("list", m)
  for (s in rev(ord)) below[[s]] <- c(s, unlist(below[children[[s]]], use.names = FALSE))
  meet <- matrix(0, m, m)
  # sites from the outlet upwards, so the deepest common site is written last
  for (s in ord) meet[below[[s]], below[[s]]] <- net$distance[s]
  net$distance - meet
}
