# options and session state for one test, restored when the test ends
local_cite <- function(env = parent.frame()) {
  old <- options(comdynamite.cite = TRUE)
  state <- mget(ls(cite_state), envir = cite_state)
  rm(list = ls(cite_state), envir = cite_state)
  restore <- function() {
    options(old)
    rm(list = ls(cite_state), envir = cite_state)
    if (length(state)) list2env(state, envir = cite_state)
  }
  do.call(on.exit, list(substitute(f(), list(f = restore)), add = TRUE), envir = env)
}

test_that("a method function shows its references once per session", {
  local_cite()
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  expect_message(node_contrasts(tree), class = "comdynamite_citation")
  expect_message(node_contrasts(tree), regexp = NA)
})

test_that("the message names the function, its papers and the package", {
  local_cite()
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  message <- tryCatch(node_contrasts(tree), comdynamite_citation = conditionMessage)
  expect_match(message, "node_contrasts() implements", fixed = TRUE)
  expect_match(message, "Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010)", fixed = TRUE)
  expect_match(message, "https://doi.org/10.1111/j.1461-0248.2010.01523.x", fixed = TRUE)
  expect_match(message, "citation(\"comdynamite\")", fixed = TRUE)
})

test_that("options(comdynamite.cite = FALSE) and suppressMessages() silence it", {
  local_cite()
  options(comdynamite.cite = FALSE)
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  expect_silent(node_contrasts(tree))
  options(comdynamite.cite = TRUE)
  expect_silent(suppressMessages(node_contrasts(tree)))
  expect_message(node_contrasts(tree), regexp = NA)
})

test_that("a method function called by another one does not cite itself", {
  local_cite()
  outer <- function() cite_method("node_test")
  inner <- function() cite_method("node_contrasts")
  node_test <- function() {
    outer()
    inner()
  }
  seen <- new.env()
  seen$messages <- character()
  withCallingHandlers(node_test(), comdynamite_citation = function(m) {
    seen$messages <- c(seen$messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  expect_length(seen$messages, 0)
  expect_false(isTRUE(cite_state[["node_contrasts"]]))
})

test_that("method_references() lists the papers of each function", {
  refs <- method_references("node_test")
  expect_length(refs, 2)
  expect_match(refs[2], "Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012)", fixed = TRUE)
  expect_match(method_references("food_web_dynamics")[1], "Princeton University Press.",
               fixed = TRUE)
  all <- method_references()
  expect_named(all, c("fun", "reference"))
  expect_setequal(unique(all$fun), names(method_functions))
  expect_error(method_references("plot_rates"), "not a function that implements")
  expect_error(method_references(1), "`fun`")
})

test_that("every function in the registry exists and every paper is used", {
  expect_true(all(names(method_functions) %in% getNamespaceExports("comdynamite")))
  used <- unique(unlist(method_functions))
  expect_setequal(used, names(method_papers))
})

test_that("inst/CITATION has the package and every method paper", {
  path <- system.file("CITATION", package = "comdynamite")
  skip_if(path == "")
  entries <- utils::readCitationFile(path, meta = list(Version = "0.0.0", Encoding = "UTF-8"))
  dois <- vapply(unclass(entries), function(x) x$doi, character(1))
  expect_true("10.5281/zenodo.23993533" %in% dois)
  expect_setequal(setdiff(dois, "10.5281/zenodo.23993533"),
                  vapply(method_papers, `[[`, character(1), "doi"))
})

test_that("call_name() drops the namespace and handles anonymous functions", {
  expect_identical(call_name(quote(node_test(x))), "node_test")
  expect_identical(call_name(quote(comdynamite::node_test(x))), "node_test")
  expect_identical(call_name(quote(comdynamite:::cite_method(x))), "cite_method")
  expect_identical(call_name(quote((function(x) x)(1))), "")
})
