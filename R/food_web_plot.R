#' Plot the dynamics of a food web
#'
#' Draws the trajectories from [food_web_dynamics()] in up to three stacked panels sharing
#' the time axis: abundances, mean traits and rates of trait change
#' (\eqn{|d\bar z / dt|}). Levels that are absent during the whole run are left out, and the
#' times of `introductions` are marked by dashed vertical lines.
#'
#' With `time`, the curves stop at that time and a point marks the current state, while the
#' axes keep the limits of the whole run, so that successive plots are frames of an animation.
#' [animate_food_web()] uses this to write a GIF.
#'
#' @param x Result of [food_web_dynamics()].
#' @param panels Panels to draw, in order: any of `"abundance"`, `"trait"` and `"rate"`.
#' @param time Optional time at which to stop the curves; `NULL` draws the whole run.
#' @param col Colours of the resource, the consumer and the predator.
#' @param log_abundance Draw abundances on a logarithmic axis (zeros are left out)?
#' @param legend Add a legend to the first panel?
#' @param lwd Width of the lines.
#' @param ... Other graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @seealso [animate_food_web()], [food_web_dynamics()]
#' @examples
#' x <- food_web_dynamics(seq(0, 4000, by = 20), init = c(R = 50, N = 0.1, P = 0.1),
#'                        traits = c(R = 10, N = 4, P = 6), introductions = c(P = 2000),
#'                        dt = 0.5)
#' plot(x)
#' # one frame: the state at time 1000, on the axes of the whole run
#' plot(x, panels = c("abundance", "trait"), time = 1000)
#' @export
plot.food_web <- function(x, panels = c("abundance", "trait", "rate"), time = NULL,
                          col = c(R = "#009E73", N = "#0072B2", P = "#D55E00"),
                          log_abundance = FALSE, legend = TRUE, lwd = 2, ...) {
  check_food_web(x)
  panels <- unique(match.arg(panels, several.ok = TRUE))
  if (length(col) != 3) stop("`col` must give three colours, for R, N and P", call. = FALSE)
  if (!is.null(time) && (!is.numeric(time) || length(time) != 1 || !is.finite(time))) {
    stop("`time` must be one number", call. = FALSE)
  }
  col <- stats::setNames(unname(col), c("R", "N", "P"))
  levels <- c("R", "N", "P")[vapply(c("R", "N", "P"), function(l) any(x[[l]] > 0), logical(1))]
  if (length(levels) == 0) stop("every level is absent during the run", call. = FALSE)
  shown <- if (is.null(time)) rep(TRUE, nrow(x)) else x$time <= time
  arrivals <- attr(x, "introductions")
  labels <- c(abundance = "Abundance", trait = "Mean trait", rate = "Rate of trait change")
  columns <- list(abundance = levels, trait = paste0("z_", levels),
                  rate = paste0("rate_", levels))

  old <- graphics::par(mfrow = c(length(panels), 1), mar = c(0.5, 4.5, 0.5, 1),
                       oma = c(4, 0, 1, 0))
  on.exit(graphics::par(old))
  for (k in seq_along(panels)) {
    panel <- panels[k]
    values <- as.matrix(x[columns[[panel]]])
    # an absent level has no trait or rate to show
    values[as.matrix(x[levels]) <= 0] <- NA
    log_axis <- panel == "abundance" && log_abundance
    if (log_axis) values[values <= 0] <- NA
    ylim <- range(values, na.rm = TRUE)
    graphics::plot(range(x$time), ylim, type = "n", xaxt = "n", xlab = "",
                   ylab = labels[[panel]], log = if (log_axis) "y" else "", las = 1, ...)
    if (k == length(panels)) {
      graphics::axis(1)
      graphics::mtext("Time", side = 1, line = 2.5, cex = 0.8)
    }
    if (length(arrivals) > 0) graphics::abline(v = arrivals, lty = 2, col = "grey60")
    for (j in seq_along(levels)) {
      graphics::lines(x$time[shown], values[shown, j], col = col[levels[j]], lwd = lwd)
      if (!is.null(time) && any(shown)) {
        last <- max(which(shown))
        graphics::points(x$time[last], values[last, j], pch = 19, col = col[levels[j]])
      }
    }
    if (k == 1 && legend) {
      names_shown <- c(R = "Resource", N = "Consumer", P = "Predator")[levels]
      graphics::legend("topright", legend = names_shown, col = col[levels], lwd = lwd,
                       bg = grDevices::adjustcolor("white", 0.85), box.col = NA,
                       cex = 0.8)
    }
  }
  if (!is.null(time)) {
    graphics::mtext(paste("Time", format(time)), side = 3, outer = TRUE, line = -0.2,
                    cex = 0.8, adj = 0.02)
  }
  invisible(x)
}

#' Animate the dynamics of a food web
#'
#' Writes a GIF in which the trajectories from [food_web_dynamics()] unfold over time, one
#' frame per time step shown, with the panels and axes of [plot.food_web()]. It needs the
#' package gifski.
#'
#' @param x Result of [food_web_dynamics()].
#' @param file Path of the GIF to write.
#' @param frames Number of frames, taken at evenly spaced rows of `x`; it cannot exceed the
#'   number of rows.
#' @param delay Seconds between frames.
#' @param width,height Size of the GIF in pixels.
#' @param ... Other arguments of [plot.food_web()], such as `panels` or `col`.
#' @return The path of the GIF, invisibly.
#' @seealso [plot.food_web()]
#' @examples
#' \donttest{
#' if (requireNamespace("gifski", quietly = TRUE)) {
#'   x <- food_web_dynamics(seq(0, 3000, by = 30), dt = 0.5)
#'   gif <- animate_food_web(x, tempfile(fileext = ".gif"), frames = 30)
#'   file.exists(gif)
#' }
#' }
#' @export
animate_food_web <- function(x, file = "food_web.gif", frames = 100, delay = 0.08,
                             width = 640, height = 640, ...) {
  check_food_web(x)
  if (!has_gifski()) {
    stop("`animate_food_web()` needs the package gifski: install.packages(\"gifski\")",
         call. = FALSE)
  }
  if (!is.numeric(frames) || length(frames) != 1 || frames < 2 || frames > nrow(x)) {
    stop("`frames` must be one number from 2 to the number of rows of `x` (", nrow(x), ")",
         call. = FALSE)
  }
  if (!is.numeric(delay) || length(delay) != 1 || delay <= 0) {
    stop("`delay` must be one positive number", call. = FALSE)
  }
  rows <- unique(round(seq(1, nrow(x), length.out = frames)))
  gifski::save_gif(for (r in rows) plot.food_web(x, time = x$time[r], ...),
                   gif_file = file, width = width, height = height, delay = delay,
                   progress = FALSE)
  invisible(file)
}

has_gifski <- function() requireNamespace("gifski", quietly = TRUE)

check_food_web <- function(x) {
  needed <- c("time", "R", "N", "P", "z_R", "z_N", "z_P", "rate_R", "rate_N", "rate_P")
  if (!inherits(x, "food_web") || !all(needed %in% names(x))) {
    stop("`x` must be the result of `food_web_dynamics()`, with all its columns",
         call. = FALSE)
  }
  invisible(x)
}
