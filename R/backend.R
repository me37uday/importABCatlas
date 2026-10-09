# Internal subprocess bridge. No embedded Python session and no automatic installs.
.default_release <- "20260711"

# Resolve Python deterministically. Keep this deliberately simple: do not let
# reticulate, package load hooks, or a previously discovered interpreter mutate
# the choice. An explicit argument always wins.
.resolve_abc_python <- function(python = NULL) {
  pick <- function(x) {
    if (is.null(x) || length(x) == 0L || is.na(x[[1L]]) || !nzchar(x[[1L]])) return(NULL)
    path.expand(as.character(x[[1L]]))
  }

  path <- pick(python)
  if (is.null(path)) path <- pick(Sys.getenv("IMPORTABCATLAS_PYTHON", unset = ""))
  if (is.null(path)) path <- pick(getOption("importABCatlas.python", ""))
  if (is.null(path)) path <- pick(Sys.getenv("RETICULATE_PYTHON", unset = ""))
  if (is.null(path)) path <- pick(Sys.getenv("ABCATLAS_PYTHON", unset = "")) # legacy

  if (is.null(path)) {
    on_path <- unname(Sys.which(c("python3", "python")))
    on_path <- on_path[nzchar(on_path)]
    if (length(on_path)) path <- on_path[[1L]]
  }

  if (is.null(path)) {
    stop("Python not found. Run setup_environment(install=TRUE), set IMPORTABCATLAS_PYTHON, or set options(importABCatlas.python='/absolute/path/to/python').", call. = FALSE)
  }
  if (!file.exists(path)) stop(sprintf("Configured Python does not exist: %s", path), call. = FALSE)
  # Preserve virtualenv/venv symlink paths; resolving them bypasses environment site-packages.
  path
}

# Backwards-compatible internal name used by existing code/tests.
.abc_python <- function(python = NULL) .resolve_abc_python(python = python)

.abc_backend <- function(action, args = list(), reader = function(out, result) result,
                         python = NULL, show_progress = action %in% c("load", "fetch", "download_file")) {
  script <- system.file("python", "abc_backend.py", package = "importABCatlas")
  if (!nzchar(script)) stop("Backend missing: install the package or use devtools::load_all().", call. = FALSE)

  # Resolve exactly once and pass the resulting executable directly to system2.
  # This prevents a backend call from rediscovering a different Python.
  python_exe <- .resolve_abc_python(python = python)

  out <- tempfile("importABCatlas-")
  dir.create(out, recursive = TRUE)
  ok <- FALSE
  on.exit({ if (ok) unlink(out, recursive = TRUE) }, add = TRUE)
  request <- file.path(out, "request.json")
  jsonlite::write_json(list(action = action, args = args, out_dir = out), request,
                       auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  log <- file.path(out, "backend.log")
  # Allen uses tqdm-style progress output for downloads. In interactive download-capable
  # calls, inherit stdout/stderr so bytes, rate and ETA remain visible to the R user.
  # The backend still writes error.json, so structured failures are retained.
  if (isTRUE(show_progress) && interactive()) {
    status <- system2(python_exe, c(shQuote(script), shQuote(request)), stdout = "", stderr = "")
  } else {
    status <- system2(python_exe, c(shQuote(script), shQuote(request)), stdout = log, stderr = log)
  }
  if (!identical(as.integer(status), 0L)) {
    error_file <- file.path(out, "error.json")
    msg <- if (file.exists(error_file)) jsonlite::read_json(error_file)$message else paste(readLines(log, warn = FALSE), collapse = "\n")
    stop(sprintf("ABC backend failed using Python %s: %s\nDiagnostic files retained at: %s",
                 python_exe, msg, out), call. = FALSE)
  }
  result <- jsonlite::read_json(file.path(out, "result.json"), simplifyVector = FALSE)
  value <- reader(out, result)
  ok <- TRUE
  value
}

.read_frame <- function(path) {
  x <- utils::read.csv(path, colClasses = "character", check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = "", fileEncoding = "UTF-8")
  schema <- jsonlite::read_json(paste0(path, ".schema.json"), simplifyVector = TRUE)
  for (nm in intersect(names(x), names(schema))) {
    if (schema[[nm]] == "numeric") x[[nm]] <- as.numeric(x[[nm]])
    if (schema[[nm]] == "logical") x[[nm]] <- tolower(x[[nm]]) == "true"
  }
  x
}

.need_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is needed for this return_type. See docs/BEGINNER_GUIDE.md for installation, or choose return_type='matrix'.", package), call. = FALSE)
  }
}

.cache_args <- function(download_base, release, offline) {
  list(download_base = normalizePath(path.expand(download_base), mustWork = FALSE, winslash = "/"),
       release = as.character(release), offline = isTRUE(offline))
}
