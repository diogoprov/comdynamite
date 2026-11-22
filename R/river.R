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
#' @param n_sites Number of sites: one whole number, 2 or more.
#' @param branching Probability that a new site starts a branch rather than extending a
#'   headwater: one number between 0 and 1.
#' @param reach_length Mean length of the reach between a site and its downstream neighbour:
#'   one number greater than 0.
#'   Lengths are drawn from an exponential distribution with this mean, or are all equal with
#'   `fixed_length = TRUE`.
#' @param fixed_length Equal reach lengths? `TRUE` or `FALSE`.
#' @param max_branches Largest number of branches that join at a confluence (2 for binary
#'   confluences): one whole number, 2 or more.
#' @return A data frame with one row per site, the outlet first: `site`, `downstream` (row of
#'   the downstream neighbour, `NA` at the outlet), `length` (of the reach to the downstream
#'   neighbour, 0 at the outlet), `distance` (along the channel to the outlet), `magnitude`
#'   (Shreve magnitude, the number of headwaters upstream of the site, itself included),
#'   `upstream_length` (the length of channel upstream of the site, its own reach included, a
#'   proxy for the catchment area, which describes the size of a river better than Strahler
#'   order; Altermatt 2013), `order` (Strahler order), `headwater` (logical) and the coordinates `x`
#'   and `y`.
#' @references Altermatt, F. (2013). Diversity in riverine metacommunities: a network
#'   perspective. *Aquatic Ecology*, 47, 365--377. \doi{10.1007/s10452-013-9450-3}
#'
#'   Rinaldo, A., Rigon, R., Banavar, J. R., Maritan, A., & Rodriguez-Iturbe, I.
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
  check_count(n_sites, "n_sites", min = 2)
  check_number(branching, "branching", 0, 1)
  check_number(reach_length, "reach_length", lower = 0, open = TRUE)
  check_flag(fixed_length, "fixed_length")
  check_count(max_branches, "max_branches", min = 2)
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
#'   upstream, in the units of `length`: each one number greater than 0, or `Inf` for no cost
#'   of travel in that direction.
#' @param overland Share of emigrants that move in a straight line rather than along the
#'   channel: one number between 0 and 1.
#' @param overland_scale Scale of the exponential kernel of straight-line dispersal: one number
#'   greater than 0.
#' @param env_range Range of the environment when `env` is `NULL`: two finite numbers, the
#'   first smaller.
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
  if (!is.data.frame(network) || !all(c("downstream", "length") %in% names(network))) {
    stop("`network` must be a data frame with columns `downstream` and `length`",
         call. = FALSE)
  }
  # an infinite scale removes the cost of travel in that direction
  check_number(scale_down, "scale_down", lower = 0, open = c(TRUE, FALSE),
               finite = FALSE)
  check_number(scale_up, "scale_up", lower = 0, open = c(TRUE, FALSE),
               finite = FALSE)
  check_number(overland, "overland", 0, 1)
  check_number(overland_scale, "overland_scale", lower = 0, open = TRUE)
  check_env_range(env_range)
  net <- river_given_coords(river_table(network$downstream, network$length), network)
  env <- river_site_env(env, net, env_range)
  down <- river_downstream(net)
  channel <- river_channel(down, scale_down, scale_up)
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

#' the network table `net` with the coordinates given in the columns `x` and `y` of `network`,
#' checked, in place of the layout of river_table(); `net` as it is without those columns
#' @noRd
river_given_coords <- function(net, network) {
  if (!all(c("x", "y") %in% names(network))) return(net)
  xy <- network[, c("x", "y")]
  if (!is.numeric(xy$x) || !is.numeric(xy$y) || !all(is.finite(c(xy$x, xy$y)))) {
    stop("columns `x` and `y` of `network` must be finite numbers", call. = FALSE)
  }
  net$x <- xy$x
  net$y <- xy$y
  net
}

#' the environment of each site of a river: `env` checked, or, when NULL, a gradient over
#' `env_range` with the distance to the outlet
#' @noRd
river_site_env <- function(env, net, env_range) {
  m <- nrow(net)
  if (is.null(env)) {
    env <- net$distance / max(net$distance) * diff(env_range) + env_range[1]
  }
  if (!is.numeric(env) || length(env) != m || !all(is.finite(env))) {
    stop("`env` must be NULL or one finite number per site of `network` (", m, ")",
         call. = FALSE)
  }
  env
}

#' dispersal along the channel: an exponential kernel of the distances travelled downstream
#' (`down`) and upstream, each with its scale; rows sum to one and the diagonal is zero
#' @noRd
river_channel <- function(down, scale_down, scale_up) {
  w <- -down / scale_down - t(down) / scale_up
  diag(w) <- -Inf
  channel <- exp(w - apply(w, 1, max))
  channel / rowSums(channel)
}

#' the network table of river_network() from the downstream neighbour and reach length of each
#' site; checks that the sites form a single tree draining to one outlet
#' @noRd
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
  children <- split(seq_len(m)[-outlet], factor(downstream[-outlet], levels = seq_len(m)))
  ord <- river_upwards(children, outlet)
  distance <- numeric(m)
  for (s in ord[-1]) distance[s] <- distance[downstream[s]] + length[s]
  size <- river_size(children, ord, length)
  angle <- river_angle(children, ord, size$magnitude, outlet)
  data.frame(site = paste0("site", seq_len(m)), downstream = downstream, length = length,
             distance = distance, magnitude = size$magnitude, upstream_length = size$upstream,
             order = as.integer(size$order),
             headwater = lengths(children) == 0, x = distance * cos(angle),
             y = distance * sin(angle), stringsAsFactors = FALSE)
}

#' sites from the outlet upwards, every site after its downstream neighbour; stops unless the
#' `children` of each site form a single tree draining to the `outlet`
#' @noRd
river_upwards <- function(children, outlet) {
  ord <- outlet
  i <- 1
  while (i <= length(ord)) {
    ord <- c(ord, children[[ord[i]]])
    i <- i + 1
  }
  if (length(ord) != length(children)) {
    stop("`downstream` must form a tree draining to the outlet, without loops", call. = FALSE)
  }
  ord
}

#' size of the network above each site: the `magnitude` (number of headwaters), the `upstream`
#' length of the reaches, the site included, and the Strahler `order`
#' @noRd
river_size <- function(children, ord, reach) {
  n_children <- lengths(children)
  magnitude <- as.numeric(n_children == 0)
  upstream <- reach
  order <- as.integer(n_children == 0)
  for (s in rev(ord)) {
    ch <- children[[s]]
    if (!length(ch)) next
    magnitude[s] <- sum(magnitude[ch])
    upstream[s] <- upstream[s] + sum(upstream[ch])
    top <- max(order[ch])
    order[s] <- top + as.integer(sum(order[ch] == top) > 1)
  }
  list(magnitude = magnitude, upstream = upstream, order = order)
}

#' angular layout: the angle of each site, in the middle of its wedge, which is split among
#' the branches above it by their magnitude
#' @noRd
river_angle <- function(children, ord, magnitude, outlet) {
  lo <- hi <- numeric(length(children))
  lo[outlet] <- pi / 6
  hi[outlet] <- 5 * pi / 6
  for (s in ord) {
    ch <- children[[s]]
    if (!length(ch)) next
    cuts <- lo[s] + (hi[s] - lo[s]) * c(0, cumsum(magnitude[ch])) / sum(magnitude[ch])
    lo[ch] <- cuts[-length(cuts)]
    hi[ch] <- cuts[-1]
  }
  (lo + hi) / 2
}

#' distance travelled downstream from each site (rows) to each other site (columns): from the
#' site to the confluence where their paths meet
#' @noRd
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

#' Covariance of the environment along a river network
#'
#' Covariance among the sites of a river network as a sum of exponential components, the
#' variance-component approach of Peterson and Ver Hoef (2010) for stream networks:
#'
#' * `tail_up`: correlation only between flow-connected sites (one downstream of the other),
#'   \eqn{\sigma^2 \exp(-h / \alpha) \sqrt{A_u / A_d}}{sigma^2 exp(-h / alpha) sqrt(A_u / A_d)},
#'   where \eqn{h}{h} is the distance along the channel and \eqn{A_u}{A_u} and \eqn{A_d}{A_d}
#'   are the additive function at the upstream and downstream site (Ver Hoef et al. 2006). It
#'   suits variables carried by the flow, such as water chemistry or temperature, which mix at
#'   confluences in proportion to the flow of each branch.
#' * `tail_down`: correlation between all sites, \eqn{\sigma^2 \exp(-h / \alpha)}{sigma^2
#'   exp(-h / alpha)}, which for the exponential model is the same for flow-connected and
#'   flow-unconnected pairs at the same distance along the channel (Ver Hoef and Peterson
#'   2010). It suits processes that spread up and down the channel.
#' * `euclidean`: \eqn{\sigma^2 \exp(-d / \alpha)}{sigma^2 exp(-d / alpha)} with the
#'   straight-line distance \eqn{d}{d}, for drivers of the terrestrial landscape and the
#'   atmosphere, such as climate or land use, that ignore the channel.
#'
#' The additive function comes from weights of the branches: at each confluence the
#' moving average of the tail-up model splits among the branches upstream in proportion to
#' their weights, and the additive function of a site is the product of these shares down to
#' the outlet, which keeps the variance equal at all sites. With `weights = "magnitude"` (the
#' default) the shares follow the Shreve magnitude, with `"upstream_length"` the length of
#' channel upstream (a proxy for catchment area or discharge), and a numeric vector gives
#' any other positive measure per site, such as discharge.
#'
#' @param landscape A landscape from [landscape_river()].
#' @param tail_up,tail_down,euclidean `NULL` (the default) to leave the component out, or a
#'   numeric vector with the partial sill (variance) and the range \eqn{\alpha}{alpha}, named
#'   `sill` and `range` or in that order. Ranges are in the units of the reach lengths (tail-up
#'   and tail-down) or of the coordinates (Euclidean).
#' @param nugget Variance independent among sites: one number, 0 or more.
#' @param weights Weights of the branches for the tail-up component: `"magnitude"`,
#'   `"upstream_length"` or a positive numeric vector with one value per site.
#' @return A symmetric sites x sites covariance matrix.
#' @references Peterson, E. E., & Ver Hoef, J. M. (2010). A mixed-model moving-average approach
#'   to geostatistical modeling in stream networks. *Ecology*, 91, 644--651.
#'   \doi{10.1890/08-1668.1}
#'
#'   Ver Hoef, J. M., Peterson, E., & Theobald, D. (2006). Spatial statistical models that use
#'   flow and stream distance. *Environmental and Ecological Statistics*, 13, 449--464.
#'   \doi{10.1007/s10651-006-0022-8}
#'
#'   Ver Hoef, J. M., & Peterson, E. E. (2010). A moving average approach for spatial
#'   statistical models of stream networks. *Journal of the American Statistical
#'   Association*, 105, 6--18. \doi{10.1198/jasa.2009.ap08248}
#' @family landscapes
#' @examples
#' set.seed(1)
#' land <- landscape_river(river_network(30))
#' v <- river_covariance(land, tail_up = c(sill = 1, range = 5), nugget = 0.1)
#' # flow-unconnected sites are uncorrelated under the tail-up model
#' unconnected <- land$downstream > 0 & t(land$downstream) > 0
#' range(v[unconnected])
#' @export
river_covariance <- function(landscape, tail_up = NULL, tail_down = NULL, euclidean = NULL,
                             nugget = 0, weights = "magnitude") {
  check_river_landscape(landscape)
  check_number(nugget, "nugget", lower = 0)
  net <- landscape$network
  m <- nrow(net)
  tail_up <- covariance_component(tail_up, "tail_up")
  tail_down <- covariance_component(tail_down, "tail_down")
  euclidean <- covariance_component(euclidean, "euclidean")
  no_component <- is.null(tail_up) && is.null(tail_down) && is.null(euclidean)
  if (no_component && nugget == 0) {
    stop("give at least one component or a nugget", call. = FALSE)
  }
  out <- diag(nugget, m)
  h <- landscape$dist
  if (!is.null(tail_down)) out <- out + tail_down[1] * exp(-h / tail_down[2])
  if (!is.null(euclidean)) {
    d <- unname(as.matrix(stats::dist(landscape$coords)))
    out <- out + euclidean[1] * exp(-d / euclidean[2])
  }
  if (!is.null(tail_up)) {
    af <- river_additive(net, weights)
    down <- landscape$downstream
    connected <- down == 0 | t(down) == 0
    ratio <- outer(af, af, pmin) / outer(af, af, pmax)
    out <- out + tail_up[1] * exp(-h / tail_up[2]) * sqrt(ratio) * connected
  }
  out
}

#' a component of river_covariance() as its sill and range, in that order and without names,
#' checked; NULL when not given
#' @noRd
covariance_component <- function(x, name) {
  if (is.null(x)) return(NULL)
  not_pair <- !is.numeric(x) || length(x) != 2
  if (not_pair || !all(is.finite(x)) || any(x <= 0)) {
    stop("`", name, "` must be NULL or two positive numbers, the sill and the range",
         call. = FALSE)
  }
  if (!is.null(names(x)) && all(c("sill", "range") %in% names(x))) x <- x[c("sill", "range")]
  unname(x)
}

#' Environment along a river network
#'
#' Draws environments for the sites of a river network from a Gaussian field with the
#' covariance of [river_covariance()], for [simulate_metacommunity()]: a single value per site
#' or, with `n_steps`, a series that changes in time as a first-order autoregressive process
#' with the same spatial covariance at every step, as `env_time`.
#'
#' @inheritParams river_covariance
#' @param n_steps `NULL` (the default) for one value per site, or the number of steps of a
#'   series: one whole number, 1 or more.
#' @param mean Mean of the environment (one value, or one per site, such as a gradient from
#'   the outlet to the headwaters).
#' @param temporal_scale Steps over which the correlation in time falls by a factor of
#'   \eqn{e}{e}, one number, 0 or more, or `Inf` (the correlation between consecutive steps is
#'   \eqn{\exp(-1 /}{exp(-1 /}
#'   `temporal_scale`)); `Inf` keeps the first draw at every step and 0 draws every step
#'   anew.
#' @param ... Components of the covariance, passed to [river_covariance()] (`tail_up`,
#'   `tail_down`, `euclidean`, `nugget` and `weights`).
#' @return A numeric vector with one value per site, or with `n_steps` a matrix of sites x
#'   steps.
#' @family landscapes
#' @examples
#' set.seed(1)
#' land <- landscape_river(river_network(40))
#' # water temperature: carried by the flow, warmer towards the outlet
#' temp <- env_river(land, tail_up = c(4, 10), mean = 25 - land$network$distance / 2)
#' land$env <- temp
#' # a series that changes over about 50 steps, for env_time
#' series <- env_river(land, n_steps = 200, tail_up = c(4, 10), nugget = 0.5,
#'                     mean = 25 - land$network$distance / 2, temporal_scale = 50)
#' dim(series)
#' @export
env_river <- function(landscape, n_steps = NULL, mean = 0, temporal_scale = Inf, ...) {
  check_count(n_steps, "n_steps", null_ok = TRUE)
  # Inf keeps the first draw at every step
  check_number(temporal_scale, "temporal_scale", lower = 0, finite = FALSE)
  v <- river_covariance(landscape, ...)
  m <- nrow(v)
  if (!is.numeric(mean) || !length(mean) %in% c(1, m) || !all(is.finite(mean))) {
    stop("`mean` must be one finite number or one per site (", m, ")", call. = FALSE)
  }
  ev <- eigen(v, symmetric = TRUE)
  root <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), m)
  draw <- function() as.vector(root %*% stats::rnorm(m))
  if (is.null(n_steps)) return(mean + draw())
  rho <- exp(-1 / temporal_scale)
  out <- matrix(0, m, n_steps)
  out[, 1] <- draw()
  if (n_steps > 1) for (s in 2:n_steps) {
    out[, s] <- if (rho == 1) out[, s - 1] else rho * out[, s - 1] + sqrt(1 - rho^2) * draw()
  }
  out + mean
}

#' stops unless `landscape` comes from landscape_river()
#' @noRd
check_river_landscape <- function(landscape) {
  if (!inherits(landscape, "mc_landscape") || is.null(landscape$network)) {
    stop("`landscape` must come from landscape_river()", call. = FALSE)
  }
  invisible(landscape)
}

#' additive function of the tail-up model: 1 at the outlet, and at each confluence the share of
#' each branch upstream in proportion to its weight, multiplied down to the outlet
#' @noRd
river_additive <- function(net, weights) {
  m <- nrow(net)
  w <- if (is.character(weights)) {
    weights <- match.arg(weights, c("magnitude", "upstream_length"))
    net[[weights]]
  } else weights
  if (!is.numeric(w) || length(w) != m || !all(is.finite(w)) || any(w <= 0)) {
    stop("`weights` must be \"magnitude\", \"upstream_length\" or one positive value per site",
         call. = FALSE)
  }
  outlet <- which(is.na(net$downstream))
  children <- split(seq_len(m)[-outlet], factor(net$downstream[-outlet], levels = seq_len(m)))
  af <- numeric(m)
  af[outlet] <- 1
  ord <- outlet
  i <- 1
  while (i <= length(ord)) {
    s <- ord[i]
    ch <- children[[s]]
    if (length(ch)) af[ch] <- af[s] * w[ch] / sum(w[ch])
    ord <- c(ord, ch)
    i <- i + 1
  }
  af
}

#' Seasonal floods along a river network
#'
#' Probabilities that each site of a river network loses all its populations at each step
#' because of floods, for `disturbance` in [simulate_metacommunity()]. Each period (a year)
#' has a flood season of `flood_steps` steps, starting at step `start` of the period, and a
#' flood happens in a given year with probability `year_prob`, so that years differ. During
#' a flood, the probability of washout at a site is proportional to its size relative to the
#' largest site raised to `size_effect`, scaled so that the most exposed site has probability
#' `probability`: positive values make floods harsher in large rivers, negative values in
#' headwaters (flashier floods in small, steep streams), and 0 the same everywhere.
#'
#' Floods here remove populations; to change the environment during floods instead, use
#' [env_pulse()] on the sites and steps of the flood.
#'
#' @param landscape A landscape from [landscape_river()].
#' @param n_steps Number of steps: one whole number, 1 or more.
#' @param period Steps per year: one whole number, 1 or more.
#' @param flood_steps Steps of the flood season in each year: one whole number, from 1 up to
#'   `period`.
#' @param start Step of the year at which the flood season begins: one whole number, from 1 up
#'   to `period`.
#' @param probability Probability of washout during a flood at the most exposed site: one
#'   number between 0 and 1.
#' @param size_effect Exponent of the relative size of each site: one finite number.
#' @param size Size of each site: `"upstream_length"` (the default, a proxy for catchment
#'   area and discharge), `"magnitude"` (Shreve magnitude) or a positive numeric vector with
#'   one value per site.
#' @param year_prob Probability that a flood happens in a year: one number between 0 and 1.
#' @return A matrix of probabilities with one row per site and one column per step, with an
#'   attribute `flood_years`, the years in which a flood happened.
#' @family landscapes
#' @examples
#' set.seed(1)
#' land <- landscape_river(river_network(30))
#' floods <- river_floods(land, n_steps = 200, period = 20, flood_steps = 2,
#'                        probability = 0.3, year_prob = 0.7)
#' attr(floods, "flood_years")
#' # the mainstem floods more than the headwaters
#' tapply(apply(floods, 1, max), land$network$headwater, mean)
#' traits <- data.frame(optimum = 50, sigma = 30, r_max = 3, row.names = "sp1")
#' sim <- simulate_metacommunity(land, traits, matrix(0.05, dimnames = list("sp1", "sp1")),
#'                               n_steps = 200, disturbance = floods)
#' @export
river_floods <- function(landscape, n_steps, period, flood_steps = 1, start = 1,
                         probability = 0.5, size_effect = 1, size = "upstream_length",
                         year_prob = 1) {
  check_river_landscape(landscape)
  net <- landscape$network
  m <- nrow(net)
  check_count(n_steps, "n_steps")
  check_count(period, "period")
  check_count(flood_steps, "flood_steps")
  check_count(start, "start")
  if (flood_steps > period || start > period) {
    stop("`flood_steps` and `start` must be at most `period`", call. = FALSE)
  }
  check_number(probability, "probability", 0, 1)
  check_number(size_effect, "size_effect")
  check_number(year_prob, "year_prob", 0, 1)
  w <- if (is.character(size)) {
    net[[match.arg(size, c("upstream_length", "magnitude"))]]
  } else size
  if (!is.numeric(w) || length(w) != m || !all(is.finite(w)) || any(w <= 0)) {
    stop("`size` must be \"upstream_length\", \"magnitude\" or one positive value per site",
         call. = FALSE)
  }
  rel <- (w / max(w))^size_effect
  p <- probability * rel / max(rel)
  years <- seq_len(ceiling(n_steps / period))
  flood <- years[stats::runif(length(years)) < year_prob]
  steps <- unlist(lapply(flood, function(y) (y - 1) * period + start - 1 + seq_len(flood_steps)))
  steps <- steps[steps <= n_steps]
  out <- matrix(0, m, n_steps)
  out[, steps] <- p
  attr(out, "flood_years") <- flood
  out
}
