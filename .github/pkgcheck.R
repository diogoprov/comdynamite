# Run pkgcheck inside the image of ropensci-review-tools/pkgcheck-action. This follows the
# action's own R/check.R, but flattens the markdown with unlist(): checks_to_markdown() can
# return a list when goodpractice runs (it inserts the list of the CI check), and the action
# then stops with "can only write character objects" (pkgcheck-action issue #77).
library(pkgcheck)

out_dir <- fs::dir_create(".pkgcheck")
Sys.setenv(PKGCACHE_HTTP_VERSION = "2")
pak::lockfile_create("local::.", lockfile = ".pkgcheck/pkg.lock", dependencies = "all")
pak::lockfile_install(".pkgcheck/pkg.lock")
usethis::use_build_ignore(".pkgcheck")
pkgstats::ctags_install(sudo = TRUE)

check <- pkgcheck()
print(check)

md <- unlist(checks_to_markdown(check, render = FALSE), use.names = FALSE)
s_break <- grep("---", md, fixed = TRUE)[1]
writeLines(md[seq_len(s_break - 1)], file.path(out_dir, "summary.md"))
writeLines(md, file.path(out_dir, "full.md"))
if (fs::file_exists(check$info$network_file)) fs::file_copy(check$info$network_file, out_dir)

errors <- grep("^- :heavy_multiplication_x:", md[seq_len(s_break)], value = TRUE)
errors <- sub("^- :heavy_multiplication_x: ", "", errors)
for (e in errors) message("pkgcheck failure: ", e)
quit(status = as.integer(length(errors) > 0))
