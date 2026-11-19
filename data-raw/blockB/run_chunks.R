# Shared runner for the block B studies (dev/ROADMAP.md, block B): more replicates and sweeps
# of the central parameters of the three replications. A design (one row per run) is split
# into chunks of `chunk` rows; each chunk is run in parallel and saved to `out_dir` as
# chunk_XXXXX.rds, so that an interrupted run resumes by running the script again (finished
# chunks are skipped). `collect_chunks()` binds the chunks once all are done.
# Environment: CORES (default: all cores but two), CHUNK (rows per chunk), OUT (root folder of
# the outputs, default data-raw/blockB/out), MAX_ROWS (only the first rows of each design, for
# a quick test).

block_settings <- function() {
  list(cores = as.integer(Sys.getenv("CORES", max(1, parallel::detectCores() - 2))),
       chunk = as.integer(Sys.getenv("CHUNK", "200")),
       out = Sys.getenv("OUT", "data-raw/blockB/out"),
       max_rows = as.integer(Sys.getenv("MAX_ROWS", "0")))
}

run_chunks <- function(design, f, name, settings = block_settings()) {
  if (settings$max_rows > 0) design <- utils::head(design, settings$max_rows)
  dir <- file.path(settings$out, name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(design, file.path(dir, "design.rds"))
  starts <- seq(1, nrow(design), by = settings$chunk)
  t0 <- Sys.time()
  for (c in seq_along(starts)) {
    file <- file.path(dir, sprintf("chunk_%05d.rds", c))
    if (file.exists(file)) next
    rows <- starts[c]:min(starts[c] + settings$chunk - 1, nrow(design))
    out <- parallel::mclapply(rows, function(i) {
      res <- try(f(design[i, , drop = FALSE]), silent = TRUE)
      if (inherits(res, "try-error")) {
        stop("row ", i, " of ", name, " failed: ", res, call. = FALSE)
      }
      res
    }, mc.cores = settings$cores, mc.preschedule = TRUE)
    bad <- vapply(out, inherits, logical(1), "try-error")
    if (any(bad)) stop(out[[which(bad)[1]]])
    saveRDS(do.call(rbind, out), file)
    done <- c / length(starts)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("[%s] %s: chunk %d of %d (%.0f%%), %.1f min", format(Sys.time(), "%H:%M"),
                    name, c, length(starts), 100 * done, elapsed))
  }
  invisible(dir)
}

collect_chunks <- function(name, settings = block_settings()) {
  dir <- file.path(settings$out, name)
  files <- sort(list.files(dir, "^chunk_\\d+\\.rds$", full.names = TRUE))
  design <- readRDS(file.path(dir, "design.rds"))
  expected <- ceiling(nrow(design) / settings$chunk)
  if (length(files) < expected) {
    stop(name, ": ", length(files), " of ", expected, " chunks done", call. = FALSE)
  }
  do.call(rbind, lapply(files, readRDS))
}
