# Run from the package root after reviewing the R source documentation.
# Removes only the clearly named bootstrap documentation supplied in this ZIP.
stopifnot(file.exists("DESCRIPTION"), dir.exists("R"))
if (!requireNamespace("roxygen2", quietly = TRUE)) stop("Install roxygen2 first")
namespace <- readLines("NAMESPACE", warn = FALSE)
if (length(namespace) && grepl("Bootstrap namespace", namespace[[1]], fixed = TRUE)) unlink("NAMESPACE")
rd <- "man/importABCatlas.Rd"
if (file.exists(rd) && grepl("Bootstrap help", readLines(rd, n = 1), fixed = TRUE)) unlink(rd)
roxygen2::roxygenise()
