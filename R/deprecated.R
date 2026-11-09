# Functions and arguments renamed in 0.7.0. The old names keep working, with a warning, until
# 1.0.0.

#' Deprecated functions
#'
#' Renamed in 0.7.0 for consistency; the old names still work, with a warning, until 1.0.0.
#'
#' * `traits_ou()` is now [trait_ou()], like [trait_bm()], [trait_mk()] and [trait_rates()].
#' * `fit_damocles()` is now [damocles_fit()], like the other functions of DAMOCLES.
#'
#' @param ... Passed to the new function.
#' @return The result of the new function.
#' @name comdynamite-deprecated
#' @keywords internal
NULL

#' @rdname comdynamite-deprecated
#' @export
traits_ou <- function(...) {
  .Deprecated("trait_ou", package = "comdynamite")
  trait_ou(...)
}

#' @rdname comdynamite-deprecated
#' @export
fit_damocles <- function(...) {
  .Deprecated("damocles_fit", package = "comdynamite")
  damocles_fit(...)
}

# value of a renamed argument: the old one, with a warning, when it was given
renamed_arg <- function(new, old, new_name, old_name) {
  if (is.null(old)) return(new)
  warning("`", old_name, "` is deprecated; use `", new_name, "`", call. = FALSE)
  old
}
