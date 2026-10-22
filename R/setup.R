# Managed Python environment -------------------------------------------------
#
# importABCatlas runs its Python backend in a separate process. By default the
# package provisions a private, reproducible Python 3.11 environment with
# reticulate/uv. Users do not need to locate Python, create a virtualenv, or set
# environment variables. The environment is created on first use and reused.

.iab_python_requirements_id <- function() {
  "py311-abc-da54664-20261020-v1"
}

.iab_python_config_dir <- function() {
  path.expand(tools::R_user_dir("importABCatlas", "config"))
}

.iab_python_config_file <- function() {
  file.path(.iab_python_config_dir(), "managed-python.json")
}

.iab_python_log_file <- function() {
  file.path(.iab_python_config_dir(), "python-setup.log")
}

.iab_python_subprocess_env <- function() {
  # Prevent global Python settings from contaminating the isolated backend.
  # An empty PYTHONHOME/PYTHONPATH disables inherited values without requiring
  # users to edit shell profiles or .Renviron files.
  c("PYTHONHOME=", "PYTHONPATH=", "PYTHONNOUSERSITE=1")
}

.iab_system2_python <- function(command, args = character(),
                                stdout = "", stderr = "") {
  if (.Platform$OS.type != "windows") {
    return(system2(
      command,
      args,
      stdout = stdout,
      stderr = stderr,
      env = .iab_python_subprocess_env()
    ))
  }

  vars <- c("PYTHONHOME", "PYTHONPATH", "PYTHONNOUSERSITE")
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit({
    Sys.unsetenv(vars)
    restore <- !is.na(old)
    if (any(restore)) {
      do.call(
        Sys.setenv,
        as.list(stats::setNames(old[restore], vars[restore]))
      )
    }
  }, add = TRUE)

  Sys.unsetenv(c("PYTHONHOME", "PYTHONPATH"))
  Sys.setenv(PYTHONNOUSERSITE = "1")

  system2(
    command,
    args,
    stdout = stdout,
    stderr = stderr
  )
}

.iab_python_probe <- function(python) {
  if (is.null(python) || !length(python) || is.na(python[[1L]]) || !nzchar(as.character(python[[1L]]))) {
    return(list(ok = FALSE, python = NA_character_, message = "Python executable is missing."))
  }
  python <- path.expand(as.character(python[[1L]]))
  if (!file.exists(python)) {
    return(list(ok = FALSE, python = python, message = "Python executable is missing."))
  }

  code <- paste(
    "import importlib.util, json, sys",
    "mods=['numpy','pandas','scipy','h5py','anndata','boto3','botocore','tqdm','abc_atlas_access']",
    "missing=[m for m in mods if importlib.util.find_spec(m) is None]",
    "ok=(sys.version_info >= (3,10) and not missing)",
    "print(json.dumps({'ok':ok,'version':sys.version.split()[0],'missing':missing,'executable':sys.executable}))",
    "raise SystemExit(0 if ok else 2)",
    sep = "; "
  )

  out <- tryCatch(
    .iab_system2_python(
      python,
      c("-c", shQuote(code)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )

  status <- attr(out, "status") %||% 0L
  parsed <- tryCatch(
    jsonlite::fromJSON(tail(out, 1L), simplifyVector = TRUE),
    error = function(e) NULL
  )

  if (identical(as.integer(status), 0L) && is.list(parsed) && isTRUE(parsed$ok)) {
    return(list(
      ok = TRUE,
      python = python,
      version = as.character(parsed$version),
      missing = character(),
      message = "ready"
    ))
  }

  list(
    ok = FALSE,
    python = python,
    version = if (is.list(parsed)) as.character(parsed$version %||% NA_character_) else NA_character_,
    missing = if (is.list(parsed)) as.character(parsed$missing %||% character()) else character(),
    message = paste(out, collapse = "\n")
  )
}

.iab_read_python_config <- function() {
  path <- .iab_python_config_file()
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

.iab_managed_python_status <- function(probe = TRUE) {
  cfg <- .iab_read_python_config()
  if (is.null(cfg)) {
    return(list(
      configured = FALSE,
      healthy = FALSE,
      python = NA_character_,
      requirements_id = .iab_python_requirements_id(),
      message = "The managed Python environment has not been prepared yet."
    ))
  }

  python <- as.character(cfg$python %||% NA_character_)
  current <- identical(as.character(cfg$requirements_id %||% ""),
                       .iab_python_requirements_id())
  if (!isTRUE(probe)) {
    return(list(
      configured = TRUE,
      healthy = current && !is.na(python) && nzchar(python) && file.exists(path.expand(python)),
      python = python,
      requirements_id = as.character(cfg$requirements_id %||% ""),
      message = if (current) "configured" else "The managed environment belongs to an older package requirement set."
    ))
  }

  checked <- .iab_python_probe(python)
  list(
    configured = TRUE,
    healthy = current && isTRUE(checked$ok),
    python = python,
    version = checked$version %||% NA_character_,
    missing = checked$missing %||% character(),
    requirements_id = as.character(cfg$requirements_id %||% ""),
    message = if (!current) {
      "The managed environment belongs to an older package requirement set."
    } else {
      checked$message %||% "not ready"
    }
  )
}

.iab_acquire_setup_lock <- function(lock_dir, timeout_seconds = 7200) {
  dir.create(dirname(lock_dir), recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE)) return(invisible(TRUE))

    info <- suppressWarnings(file.info(lock_dir))
    if (nrow(info) && !is.na(info$mtime[[1L]]) &&
        as.numeric(difftime(Sys.time(), info$mtime[[1L]], units = "secs")) > timeout_seconds) {
      unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }

    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > timeout_seconds) {
      stop(
        "Timed out waiting for another importABCatlas Python setup process. ",
        "Close other R sessions and run setup_importABCatlas(repair=TRUE).",
        call. = FALSE
      )
    }
    Sys.sleep(1)
  }
}

#' Prepare Python automatically
#'
#' Creates a private Python 3.11 environment for `importABCatlas` using
#' `reticulate` and `uv`. No Python path, Conda environment, Homebrew Python, or
#' administrator access is required. The environment is reused in future R
#' sessions and is automatically recreated if its cache is removed.
#'
#' This function is optional: the first function that needs Python will call it
#' automatically. Running it once after installation makes the initial setup
#' explicit and lets the user see progress before downloading atlas data. It
#' also restores managed mode in the current session if an older development
#' version left package-specific external-Python settings behind.
#'
#' @param repair Re-run managed-environment resolution even when a stored
#'   environment appears healthy.
#' @param quiet Suppress setup messages.
#' @return The managed Python executable path, invisibly.
#' @export
setup_importABCatlas <- function(repair = FALSE, quiet = FALSE) {
  # Calling this function is an explicit request to use the package-managed
  # environment. Override any stale package-specific external-mode setting in
  # the current R session; ambient RETICULATE_PYTHON is intentionally ignored.
  Sys.setenv(IMPORTABCATLAS_PYTHON_MODE = "managed")
  options(importABCatlas.python_mode = "managed")
  # Reference the public API here so R package dependency checks can verify the
  # required reticulate version without initializing Python in this process.
  if (!is.function(reticulate::py_require)) {
    stop("reticulate >= 1.41.0 is required for automatic Python setup.", call. = FALSE)
  }
  if (!is.logical(repair) || length(repair) != 1L || is.na(repair)) {
    stop("repair must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE.", call. = FALSE)
  }

  status <- .iab_managed_python_status(probe = TRUE)
  if (!isTRUE(repair) && isTRUE(status$healthy)) {
    options(importABCatlas.python = path.expand(status$python),
            importABCatlas.python_mode = "managed")
    if (!quiet) message("importABCatlas is ready. Python is managed automatically.")
    return(invisible(path.expand(status$python)))
  }

  config_dir <- .iab_python_config_dir()
  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  lock_dir <- file.path(config_dir, "setup.lock")
  .iab_acquire_setup_lock(lock_dir)
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Another process may have completed setup while this one waited.
  status <- .iab_managed_python_status(probe = TRUE)
  if (!isTRUE(repair) && isTRUE(status$healthy)) {
    options(importABCatlas.python = path.expand(status$python),
            importABCatlas.python_mode = "managed")
    if (!quiet) message("importABCatlas is ready. Python is managed automatically.")
    return(invisible(path.expand(status$python)))
  }

  helper <- system.file("python", "provision_environment.R", package = "importABCatlas")
  requirements <- system.file("python", "managed-requirements.txt", package = "importABCatlas")
  if (!nzchar(helper) || !file.exists(helper) || !nzchar(requirements) || !file.exists(requirements)) {
    stop("Managed Python setup files are missing. Reinstall importABCatlas.", call. = FALSE)
  }

  config_file <- .iab_python_config_file()
  log_file <- .iab_python_log_file()
  if (isTRUE(repair) && file.exists(config_file)) unlink(config_file, force = TRUE)
  if (file.exists(log_file)) unlink(log_file, force = TRUE)

  if (!quiet) {
    message(
      "Preparing a private Python environment for importABCatlas.\n",
      "This is a one-time setup and may take a few minutes.\n",
      "Setup log: ", log_file
    )
  }

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  args <- c(
    "--vanilla",
    shQuote(helper),
    shQuote(config_file),
    shQuote(requirements),
    shQuote(.iab_python_requirements_id()),
    shQuote(as.character(utils::packageVersion("importABCatlas")))
  )
  result <- tryCatch(
    system2(rscript, args, stdout = log_file, stderr = log_file),
    error = function(e) structure(1L, message = conditionMessage(e))
  )
  exit_status <- suppressWarnings(as.integer(result))
  if (length(exit_status) != 1L || is.na(exit_status)) exit_status <- 1L

  status <- .iab_managed_python_status(probe = TRUE)
  if (exit_status != 0L || !isTRUE(status$healthy)) {
    detail <- if (file.exists(log_file)) {
      paste(tail(readLines(log_file, warn = FALSE), 80L), collapse = "\n")
    } else {
      attr(result, "message") %||% "No setup log was produced."
    }
    stop(
      "Automatic Python setup did not complete.\n",
      "Setup log: ", log_file, "\n\n",
      detail, "\n\n",
      "Check your internet connection, then run setup_importABCatlas(repair=TRUE).",
      call. = FALSE
    )
  }

  options(importABCatlas.python = path.expand(status$python),
          importABCatlas.python_mode = "managed")
  if (!quiet) {
    message("Setup complete. No Python path needs to be set in future sessions.")
  }
  invisible(path.expand(status$python))
}

#' Configure Python for importABCatlas
#'
#' `setup_environment()` is retained for backwards compatibility. With its
#' defaults it delegates to [setup_importABCatlas()] and creates the managed
#' environment automatically. Supplying `python` opts into advanced external-
#' Python mode for the current R session.
#'
#' @param envname Deprecated; retained for compatibility.
#' @param python Optional external Python executable. Most users should leave
#'   this `NULL`.
#' @param install Whether setup should be performed.
#' @param upstream_ref Deprecated; the managed environment uses the SDK commit
#'   pinned by the package.
#' @return The selected Python executable path, invisibly.
#' @export
setup_environment <- function(envname = NULL, python = NULL, install = TRUE,
                              upstream_ref = "da54664e8a9c8e2e2ba573d1d8f590b114a149c4") {
  if (!is.null(python)) {
    path <- path.expand(as.character(python[[1L]]))
    if (!file.exists(path)) stop("The supplied Python executable does not exist: ", path, call. = FALSE)
    checked <- .iab_python_probe(path)
    if (!isTRUE(checked$ok)) {
      stop(
        "The supplied Python does not contain the required dependencies. ",
        "Use setup_importABCatlas() for automatic setup.\n",
        checked$message,
        call. = FALSE
      )
    }
    Sys.setenv(
      IMPORTABCATLAS_PYTHON_MODE = "external",
      IMPORTABCATLAS_PYTHON = path
    )
    options(importABCatlas.python = path, importABCatlas.python_mode = "external")
    return(invisible(path))
  }

  if (!isTRUE(install)) {
    status <- .iab_managed_python_status(probe = TRUE)
    if (!isTRUE(status$healthy)) {
      stop("The managed environment is not ready. Run setup_importABCatlas().", call. = FALSE)
    }
    return(invisible(path.expand(status$python)))
  }
  setup_importABCatlas(repair = FALSE)
}

#' Report the R and Python software environment
#'
#' This function never installs anything. It reports whether the package-managed
#' Python environment is ready and, when available, probes its dependencies.
#'
#' @return A list of R, package, Python, and session information.
#' @export
diagnose <- function() {
  mode <- Sys.getenv(
    "IMPORTABCATLAS_PYTHON_MODE",
    unset = as.character(getOption("importABCatlas.python_mode", "managed"))
  )
  mode <- tolower(as.character(mode[[1L]]))
  if (identical(mode, "external")) {
    fallback <- as.character(getOption("importABCatlas.python", "") %||% "")
    py <- Sys.getenv("IMPORTABCATLAS_PYTHON", unset = fallback)
    checked <- if (nzchar(py)) .iab_python_probe(py) else list(ok = FALSE, message = "No external Python is configured.")
    managed <- FALSE
  } else {
    state <- .iab_managed_python_status(probe = TRUE)
    py <- state$python %||% NA_character_
    checked <- if (isTRUE(state$healthy)) .iab_python_probe(py) else list(ok = FALSE, message = state$message)
    managed <- TRUE
  }

  has_py <- length(py) && !is.na(py[[1L]]) && nzchar(as.character(py[[1L]]))
  py_info <- list(
    mode = if (managed) "managed" else "external",
    executable = if (has_py) as.character(py[[1L]]) else NA_character_,
    available = isTRUE(checked$ok),
    version = checked$version %||% NA_character_,
    missing = checked$missing %||% character(),
    message = checked$message %||% "unknown",
    setup_command = if (!isTRUE(checked$ok)) "setup_importABCatlas()" else NA_character_,
    setup_log = .iab_python_log_file()
  )

  list(
    R = R.version.string,
    platform = R.version$platform,
    package = as.character(utils::packageVersion("importABCatlas")),
    python = py_info,
    session = utils::capture.output(utils::sessionInfo())
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
