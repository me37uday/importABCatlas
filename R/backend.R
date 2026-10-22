# Internal subprocess bridge. Python runs out-of-process; a private environment is prepared automatically on first use.
.default_release <- "20260711"

# Resolve Python deterministically. The default is a package-managed, isolated
# environment. Ambient Python settings are ignored unless the user deliberately
# enables advanced external mode. An explicit function argument always wins.
.resolve_abc_python <- function(python = NULL, auto_setup = TRUE) {
  pick <- function(x) {
    if (is.null(x) || length(x) == 0L || is.na(x[[1L]]) || !nzchar(x[[1L]])) return(NULL)
    path.expand(as.character(x[[1L]]))
  }

  explicit <- pick(python)
  if (!is.null(explicit)) {
    if (!file.exists(explicit)) {
      stop(sprintf("Configured Python does not exist: %s", explicit), call. = FALSE)
    }
    # Preserve virtualenv/managed-environment symlink paths.
    return(explicit)
  }

  mode <- Sys.getenv(
    "IMPORTABCATLAS_PYTHON_MODE",
    unset = getOption("importABCatlas.python_mode", "managed")
  )
  mode <- tolower(as.character(mode[[1L]]))
  if (!mode %in% c("managed", "external")) {
    stop("importABCatlas.python_mode must be 'managed' or 'external'.", call. = FALSE)
  }

  if (identical(mode, "external")) {
    path <- pick(Sys.getenv("IMPORTABCATLAS_PYTHON", unset = ""))
    if (is.null(path)) path <- pick(getOption("importABCatlas.python", ""))
    if (is.null(path)) {
      stop(
        "External Python mode is enabled but no executable is configured. ",
        "Either run setup_importABCatlas() to return to automatic managed mode, ",
        "or call setup_environment(python='/absolute/path/to/python').",
        call. = FALSE
      )
    }
    if (!file.exists(path)) stop(sprintf("Configured Python does not exist: %s", path), call. = FALSE)
    return(path)
  }

  status <- .iab_managed_python_status(probe = TRUE)
  if (!isTRUE(status$healthy)) {
    if (!isTRUE(auto_setup) || !isTRUE(getOption("importABCatlas.auto_setup", TRUE))) {
      stop(
        "The managed Python environment is not ready. Run setup_importABCatlas().",
        call. = FALSE
      )
    }
    setup_importABCatlas(quiet = FALSE)
    status <- .iab_managed_python_status(probe = TRUE)
  }
  if (!isTRUE(status$healthy)) {
    stop(
      "The managed Python environment is not healthy. Run ",
      "setup_importABCatlas(repair=TRUE).",
      call. = FALSE
    )
  }
  options(importABCatlas.python = path.expand(status$python),
          importABCatlas.python_mode = "managed")
  path.expand(status$python)
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
    status <- .iab_system2_python(
      python_exe,
      c(shQuote(script), shQuote(request)),
      stdout = "",
      stderr = ""
    )
  } else {
    status <- .iab_system2_python(
      python_exe,
      c(shQuote(script), shQuote(request)),
      stdout = log,
      stderr = log
    )
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
  list(
    download_base = normalizePath(path.expand(download_base), mustWork = FALSE, winslash = "/"),
    release = as.character(release),
    offline = isTRUE(offline),
    check_disk_space = isTRUE(getOption("importABCatlas.check_disk_space", TRUE)),
    disk_reserve_gb = as.numeric(getOption("importABCatlas.disk_reserve_gb", 1)),
    disk_safety_factor = as.numeric(getOption("importABCatlas.disk_safety_factor", 1.05)),
    allow_unknown_size = isTRUE(getOption("importABCatlas.allow_unknown_size", FALSE)),
    network_retries = as.integer(getOption("importABCatlas.network_retries", 5L)),
    network_backoff_seconds = as.numeric(getOption("importABCatlas.network_backoff_seconds", 2)),
    network_connect_timeout_seconds = as.integer(getOption("importABCatlas.network_connect_timeout_seconds", 30L)),
    network_read_timeout_seconds = as.integer(getOption("importABCatlas.network_read_timeout_seconds", 300L)),
    download_chunk_mb = as.numeric(getOption("importABCatlas.download_chunk_mb", 8))
  )
}
