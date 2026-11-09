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
#' @param colours Colours of the resource, the consumer and the predator.
#' @param log_abundance Draw abundances as their base-10 logarithm (zeros are left out)?
#' @param legend Draw the legend of the levels?
#' @param linewidth Width of the lines.
#' @param col,lwd,... Deprecated: `col` is now `colours` and `lwd` is `linewidth`; graphical
#'   arguments of base graphics are no longer used.
#' @return A ggplot object.
#' @seealso [animate_food_web()], [food_web_dynamics()]
#' @family coevolution in food webs
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' x <- food_web_dynamics(seq(0, 4000, by = 20), init = c(R = 50, N = 0.1, P = 0.1),
#'                        traits = c(R = 10, N = 4, P = 6), introductions = c(P = 2000),
#'                        dt = 0.5)
#' plot(x)
#' # one frame: the state at time 1000, on the axes of the whole run
#' plot(x, panels = c("abundance", "trait"), time = 1000)
#' @export
plot.food_web <- function(x, panels = c("abundance", "trait", "rate"), time = NULL,
                          colours = c(R = "#009E73", N = "#0072B2", P = "#D55E00"),
                          log_abundance = FALSE, legend = TRUE, linewidth = 0.8, ...,
                          col = NULL, lwd = NULL) {
  need_ggplot("plot() of food_web_dynamics()")
  colours <- renamed_arg(colours, col, "colours", "col")
  linewidth <- renamed_arg(linewidth, lwd, "linewidth", "lwd")
  legend <- legend_flag(legend)
  check_food_web(x)
  panels <- unique(match.arg(panels, several.ok = TRUE))
  if (length(colours) != 3) {
    stop("`colours` must give three colours, for R, N and P", call. = FALSE)
  }
  if (!is.null(time) && (!is.numeric(time) || length(time) != 1 || !is.finite(time))) {
    stop("`time` must be one number", call. = FALSE)
  }
  colours <- stats::setNames(unname(colours), c("Resource", "Consumer", "Predator"))
  levels <- c("R", "N", "P")[vapply(c("R", "N", "P"), function(l) any(x[[l]] > 0), logical(1))]
  if (length(levels) == 0) stop("every level is absent during the run", call. = FALSE)
  level_names <- c(R = "Resource", N = "Consumer", P = "Predator")
  titles <- c(abundance = if (log_abundance) "Abundance (log10)" else "Abundance",
              trait = "Mean trait", rate = "Rate of trait change")
  columns <- list(abundance = levels, trait = paste0("z_", levels),
                  rate = paste0("rate_", levels))
  long <- do.call(rbind, lapply(panels, function(panel) {
    values <- as.matrix(x[columns[[panel]]])
    # an absent level has no trait or rate to show
    values[as.matrix(x[levels]) <= 0] <- NA
    if (panel == "abundance" && log_abundance) {
      values[values <= 0] <- NA
      values <- log10(values)
    }
    data.frame(time = rep(x$time, length(levels)), value = as.vector(values),
               level = rep(unname(level_names[levels]), each = nrow(x)), panel = titles[[panel]])
  }))
  long$panel <- factor(long$panel, unname(titles[panels]))
  long$level <- factor(long$level, unname(level_names[levels]))
  shown <- if (is.null(time)) long else long[long$time <= time, , drop = FALSE]
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(shown, ggplot2::aes(.data$time, .data$value, colour = .data$level)) +
    # the whole run sets the limits, so that frames of an animation share their axes
    ggplot2::geom_blank(data = long)
  arrivals <- attr(x, "introductions")
  if (length(arrivals) > 0) {
    p <- p + ggplot2::geom_vline(xintercept = arrivals, linetype = "dashed", colour = "grey60")
  }
  p <- p + ggplot2::geom_line(linewidth = linewidth, na.rm = TRUE)
  if (!is.null(time) && nrow(shown) > 0) {
    now <- shown[shown$time == max(shown$time), , drop = FALSE]
    p <- p + ggplot2::geom_point(data = now, size = 2, na.rm = TRUE) +
      ggplot2::labs(title = paste("Time", format(time)))
  }
  p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 1, scales = "free_y",
                          strip.position = "left") +
    ggplot2::scale_colour_manual(values = colours[levels(long$level)], name = NULL,
                                 guide = if (legend) "legend" else "none") +
    ggplot2::labs(x = "Time", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.placement = "outside", legend.position = "top")
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
#' @param ... Other arguments of [plot.food_web()], such as `panels` or `colours`.
#' @return The path of the GIF, invisibly.
#' @seealso [plot.food_web()]
#' @family coevolution in food webs
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
  need_ggplot("animate_food_web()")
  gifski::save_gif(for (r in rows) print(plot.food_web(x, time = x$time[r], ...)),
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
