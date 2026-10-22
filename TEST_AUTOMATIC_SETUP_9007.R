# Test that importABCatlas works without a manually selected Python path.
# Run in a fresh R session after installing 0.2.0.9007.

test_automatic_setup_9007 <- function(run_live_plan = FALSE,
                                      cache = "~/importABCatlas_9007_setup_test_cache") {
  if (!requireNamespace("importABCatlas", quietly = TRUE)) {
    stop("Install importABCatlas before running this test.", call. = FALSE)
  }
  if (utils::packageVersion("importABCatlas") < numeric_version("0.2.0.9007")) {
    stop("This test requires importABCatlas >= 0.2.0.9007. Restart R after installation.", call. = FALSE)
  }

  old_env <- Sys.getenv(
    c("IMPORTABCATLAS_PYTHON", "IMPORTABCATLAS_PYTHON_MODE", "RETICULATE_PYTHON"),
    unset = NA_character_
  )
  old_options <- options(
    importABCatlas.python = NULL,
    importABCatlas.python_mode = "managed",
    importABCatlas.auto_setup = TRUE
  )
  on.exit({
    options(old_options)
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
    }
  }, add = TRUE)

  # Deliberately provide invalid ambient paths. Managed mode must ignore them.
  Sys.setenv(
    IMPORTABCATLAS_PYTHON = "/this/path/must/not/be/used/python",
    RETICULATE_PYTHON = "/this/path/must/not/be/used/python"
  )
  Sys.unsetenv("IMPORTABCATLAS_PYTHON_MODE")

  importABCatlas::setup_importABCatlas()
  status <- importABCatlas::diagnose()
  if (!isTRUE(status$python$available)) {
    stop("Managed Python is not healthy: ", status$python$message, call. = FALSE)
  }
  if (identical(status$python$executable, "/this/path/must/not/be/used/python")) {
    stop("Managed mode incorrectly used an ambient Python path.", call. = FALSE)
  }

  example <- importABCatlas::example_data("matrix")
  stopifnot(identical(dim(example$expression), c(4L, 3L)))

  live <- NULL
  if (isTRUE(run_live_plan)) {
    d <- importABCatlas::load_data(
      dataset = "Zhuang4",
      download_base = cache,
      release = "20260711",
      preview_n = 100L
    )
    cells <- head(d$cell_metadata$cell_label, 3L)
    live <- importABCatlas::plan_fetch(
      data = d,
      cell_ids = cells,
      genes = NULL,
      data_type = "raw"
    )
    stopifnot(length(live$files) >= 1L)
  }

  result <- list(
    package = as.character(utils::packageVersion("importABCatlas")),
    python_mode = status$python$mode,
    python_available = status$python$available,
    python_version = status$python$version,
    synthetic_dimensions = dim(example$expression),
    live_plan = live
  )
  print(result)
  invisible(result)
}
