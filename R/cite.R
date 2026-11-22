# Credit to the original papers of the methods implemented in the package. Each function
# that implements a published method calls cite_method() first; the first call of that
# function in a session shows the original references, unless
# options(comdynamite.cite = FALSE). The references are the same as in each function's
# @references; test-cite.R checks that inst/CITATION has the same DOIs.

method_papers <- list(
  leibold2010 = list(
    authors = c("Leibold, M. A.", "Economo, E. P.", "Peres-Neto, P."), year = 2010,
    title = paste("Metacommunity phylogenetics: separating the roles of environmental",
                  "filters and historical biogeography"),
    journal = "Ecology Letters", volume = "13", pages = "1290--1299",
    doi = "10.1111/j.1461-0248.2010.01523.x"),
  peresneto2012 = list(
    authors = c("Peres-Neto, P. R.", "Leibold, M. A.", "Dray, S."), year = 2012,
    title = paste("Assessing the effects of spatial contingency and environmental",
                  "filtering on metacommunity phylogenetics"),
    journal = "Ecology", volume = "93", pages = "S14--S30", doi = "10.1890/11-0494.1"),
  peresneto2006 = list(
    authors = c("Peres-Neto, P. R.", "Legendre, P.", "Dray, S.", "Borcard, D."),
    year = 2006,
    title = "Variation partitioning of species data matrices: estimation and comparison of fractions",
    journal = "Ecology", volume = "87", pages = "2614--2625",
    doi = "10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2"),
  peresneto2015 = list(
    authors = c("Peres-Neto, P. R.", "Kembel, S. W."), year = 2015,
    title = paste("Phylogenetic gradient analysis: environmental drivers of phylogenetic",
                  "variation across ecological communities"),
    journal = "Plant Ecology", volume = "216", pages = "709--724",
    doi = "10.1007/s11258-014-0405-0"),
  pillar2010 = list(
    authors = c("Pillar, V. D.", "Duarte, L. D. S."), year = 2010,
    title = "A framework for metacommunity analysis of phylogenetic structure",
    journal = "Ecology Letters", volume = "13", pages = "587--596",
    doi = "10.1111/j.1461-0248.2010.01456.x"),
  duarte2016 = list(
    authors = c("Duarte, L. D. S.", "Debastiani, V. J.", "Freitas, A. V. L.", "Pillar, V. D."),
    year = 2016,
    title = paste("Dissecting phylogenetic fuzzy weighting: theory and application in",
                  "metacommunity phylogenetics"),
    journal = "Methods in Ecology and Evolution", volume = "7", pages = "937--946",
    doi = "10.1111/2041-210X.12547"),
  debastiani2014 = list(
    authors = c("Debastiani, V. J.", "Duarte, L. D. S."), year = 2014,
    title = paste("PCPS -- an R-package for exploring phylogenetic eigenvectors across",
                  "metacommunities"),
    journal = "Frontiers of Biogeography", volume = "6", pages = NULL,
    doi = "10.21425/F5FBG22943"),
  violle2012 = list(
    authors = c("Violle, C.", "Enquist, B. J.", "McGill, B. J.", "Jiang, L.", "Albert, C. H.",
                "Hulshof, C.", "Jung, V.", "Messier, J."), year = 2012,
    title = "The return of the variance: intraspecific variability in community ecology",
    journal = "Trends in Ecology & Evolution", volume = "27", pages = "244--252",
    doi = "10.1016/j.tree.2011.11.014"),
  taudiere2016 = list(
    authors = c("Taudiere, A.", "Violle, C."), year = 2016,
    title = paste("cati: an R package using functional traits to detect and quantify",
                  "multi-level community assembly processes"),
    journal = "Ecography", volume = "39", pages = "699--708", doi = "10.1111/ecog.01433"),
  leps2011 = list(
    authors = c("Lep\u0161, J.", "de Bello, F.", "\u0160milauer, P.", "Dole\u017eal, J."),
    year = 2011,
    title = paste("Community trait response to environment: disentangling species turnover",
                  "vs intraspecific trait variability effects"),
    journal = "Ecography", volume = "34", pages = "856--863",
    doi = "10.1111/j.1600-0587.2010.06904.x"),
  messier2010 = list(
    authors = c("Messier, J.", "McGill, B. J.", "Lechowicz, M. J."), year = 2010,
    title = "How do traits vary across ecological scales? A case for trait-based ecology",
    journal = "Ecology Letters", volume = "13", pages = "838--848",
    doi = "10.1111/j.1461-0248.2010.01476.x"),
  ackerly2007 = list(
    authors = c("Ackerly, D. D.", "Cornwell, W. K."), year = 2007,
    title = paste("A trait-based approach to community assembly: partitioning of species",
                  "trait values into within- and among-community components"),
    journal = "Ecology Letters", volume = "10", pages = "135--145",
    doi = "10.1111/j.1461-0248.2006.01006.x"),
  pigot2015 = list(
    authors = c("Pigot, A. L.", "Etienne, R. S."), year = 2015,
    title = "A new dynamic null model for phylogenetic community structure",
    journal = "Ecology Letters", volume = "18", pages = "153--163", doi = "10.1111/ele.12395"),
  ackermann2004 = list(
    authors = c("Ackermann, R. R.", "Cheverud, J. M."), year = 2004,
    title = "Detecting genetic drift versus selection in human evolution",
    journal = "Proceedings of the National Academy of Sciences", volume = "101",
    pages = "17946--17951", doi = "10.1073/pnas.0405919102"),
  jetz2012 = list(
    authors = c("Jetz, W.", "Thomas, G. H.", "Joy, J. B.", "Hartmann, K.", "Mooers, A. O."),
    year = 2012, title = "The global diversity of birds in space and time",
    journal = "Nature", volume = "491", pages = "444--448", doi = "10.1038/nature11631"),
  januario2025 = list(
    authors = c("Januario, M.", "Pinsky, M. L.", "Rabosky, D. L."), year = 2025,
    title = paste("The metapopulation bridge to macroevolutionary speciation rates: a",
                  "conceptual framework and empirical test"),
    journal = "Ecology Letters", volume = "28", pages = "e70021", doi = "10.1111/ele.70021"),
  raynal2019 = list(
    authors = c("Raynal, L.", "Marin, J.-M.", "Pudlo, P.", "Ribatet, M.", "Robert, C. P.",
                "Estoup, A."), year = 2019,
    title = "ABC random forests for Bayesian parameter inference",
    journal = "Bioinformatics", volume = "35", pages = "1720--1728",
    doi = "10.1093/bioinformatics/bty867"),
  pudlo2016 = list(
    authors = c("Pudlo, P.", "Marin, J.-M.", "Estoup, A.", "Cornuet, J.-M.", "Gautier, M.",
                "Robert, C. P."), year = 2016,
    title = "Reliable ABC model choice via random forests",
    journal = "Bioinformatics", volume = "32", pages = "859--866",
    doi = "10.1093/bioinformatics/btv684"),
  thompson2020 = list(
    authors = c("Thompson, P. L.", "Guzman, L. M.", "De Meester, L.", "Horv\u00e1th, Z.",
                "Ptacnik, R.", "Vanschoenwinkel, B.", "Viana, D. S.", "Chase, J. M."),
    year = 2020,
    title = paste("A process-based metacommunity framework linking local and regional scale",
                  "community ecology"),
    journal = "Ecology Letters", volume = "23", pages = "1314--1329", doi = "10.1111/ele.13568"),
  mcpeek2017 = list(
    authors = "McPeek, M. A.", year = 2017, title = "Evolutionary Community Ecology",
    publisher = "Princeton University Press",
    doi = "10.23943/princeton/9780691088778.001.0001"),
  mcpeek2017an = list(
    authors = "McPeek, M. A.", year = 2017,
    title = paste("The ecological dynamics of natural selection: traits and the coevolution",
                  "of community structure"),
    journal = "The American Naturalist", volume = "189", pages = "E91--E117",
    doi = "10.1086/691101"),
  goldberg2006 = list(
    authors = c("Goldberg, E. E.", "Lande, R."), year = 2006,
    title = "Ecological and reproductive character displacement on an environmental gradient",
    journal = "Evolution", volume = "60", pages = "1344--1357", doi = "10.1554/05-696.1"),
  adams2007 = list(
    authors = c("Adams, D. C.", "Collyer, M. L."), year = 2007,
    title = "Analysis of character divergence along environmental gradients and other covariates",
    journal = "Evolution", volume = "61", pages = "510--515",
    doi = "10.1111/j.1558-5646.2007.00063.x")
)

# functions that implement a published method, and the papers they implement
method_functions <- list(
  node_contrasts = "leibold2010",
  node_test = c("leibold2010", "peresneto2012"),
  phylo_vector = "peresneto2012",
  trait_partition = c("peresneto2012", "peresneto2006"),
  phylo_beta = "peresneto2015",
  phylo_beta_test = "peresneto2015",
  pcc_test = "peresneto2015",
  pcps_test = c("pillar2010", "duarte2016", "debastiani2014"),
  tstats = c("violle2012", "taudiere2016"),
  community_index = "taudiere2016",
  trait_flex = "leps2011",
  nested_variance = "messier2010",
  trait_gradient = "ackerly2007",
  trait_gradient_test = "ackerly2007",
  trait_gradient_phylo = "ackerly2007",
  simulate_damocles = "pigot2015",
  damocles_fit = "pigot2015",
  damocles_test = "pigot2015",
  drift_test = "ackermann2004",
  dr_rates = "jetz2012",
  persistence = "januario2025",
  abc_fit = "raynal2019",
  abc_recovery = "raynal2019",
  abc_choose = "pudlo2016",
  simulate_metacommunity = "thompson2020",
  food_web_dynamics = c("mcpeek2017", "mcpeek2017an"),
  population_shifts = c("goldberg2006", "adams2007")
)

# functions whose citation was already shown in this session
cite_state <- new.env(parent = emptyenv())

# one paper as text, in the style of the help pages
format_paper <- function(paper) {
  authors <- paper$authors
  if (length(authors) > 1) {
    authors <- paste0(paste(authors[-length(authors)], collapse = ", "), ", & ",
                      authors[length(authors)])
  }
  source <- if (is.null(paper$journal)) {
    paste0(paper$publisher, ".")
  } else {
    paste0(paste(c(paper$journal, paper$volume, gsub("--", "-", paper$pages, fixed = TRUE)),
                 collapse = ", "), ".")
  }
  paste0(authors, " (", paper$year, "). ", paper$title, ". ", source,
         " https://doi.org/", paper$doi)
}

#' Original references of the methods in comdynamite
#'
#' Many functions of comdynamite implement methods published by other authors. The first
#' time one of these functions is used in a session, it shows a message with the original
#' papers of its method, which should be cited together with the package
#' (`citation("comdynamite")`). `method_references()` gives the same references at any
#' time.
#'
#' The message is shown once per function and session. Set
#' `options(comdynamite.cite = FALSE)` to turn it off, for example in a script that runs
#' many sessions; `suppressMessages()` also hides it. A function called by another one
#' (such as [node_contrasts()] inside [node_test()]) does not show its own message.
#'
#' @param fun Name of a function of comdynamite, such as `"node_test"`, or `NULL` for all
#'   functions that implement a published method.
#' @return A character vector of references for one function. With `fun = NULL`, a data
#'   frame with columns `fun` and `reference`, one row per function and paper.
#' @seealso [citation()]
#' @examples
#' method_references("node_test")
#' head(method_references())
#' @export
method_references <- function(fun = NULL) {
  check_string(fun, "fun", null_ok = TRUE)
  if (is.null(fun)) {
    keys <- unlist(method_functions, use.names = FALSE)
    return(data.frame(
      fun = rep(names(method_functions), lengths(method_functions)),
      reference = vapply(method_papers[keys], format_paper, character(1)),
      row.names = NULL, stringsAsFactors = FALSE
    ))
  }
  if (!fun %in% names(method_functions)) {
    stop("`", fun, "` is not a function that implements a published method; ",
         "`method_references()` lists them", call. = FALSE)
  }
  unname(vapply(method_papers[method_functions[[fun]]], format_paper, character(1)))
}

# name of the function of a call, without the namespace
call_name <- function(call) {
  fn <- call[[1]]
  if (is.call(fn) && length(fn) == 3 && as.character(fn[[1]]) %in% c("::", ":::")) {
    fn <- fn[[3]]
  }
  if (is.name(fn)) as.character(fn) else ""
}

# shows the original references of `fun` the first time it is called in the session, unless
# options(comdynamite.cite = FALSE) or `fun` was called by another function that cites
cite_method <- function(fun) {
  if (!isTRUE(getOption("comdynamite.cite", TRUE))) return(invisible(NULL))
  if (isTRUE(cite_state[[fun]])) return(invisible(NULL))
  calls <- sys.calls()
  callers <- vapply(calls[seq_len(max(0, length(calls) - 2))], call_name, character(1))
  if (any(callers %in% names(method_functions))) return(invisible(NULL))
  assign(fun, TRUE, envir = cite_state)
  text <- paste0(
    fun, "() implements the method of:\n",
    paste0("  ", method_references(fun), collapse = "\n"),
    "\nPlease cite these papers together with comdynamite (citation(\"comdynamite\")).",
    "\nThis message is shown once per session; options(comdynamite.cite = FALSE) turns it off."
  )
  condition <- structure(class = c("comdynamite_citation", "message", "condition"),
                         list(message = paste0(text, "\n"), call = NULL))
  message(condition)
  invisible(NULL)
}
