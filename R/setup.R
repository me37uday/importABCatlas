#' Configure an explicit Python virtual environment
#'
#' Does not install anything unless install is TRUE. The upstream Git reference
#' should be replaced by a tested commit for a paper release; main is not a lock.
#' @param envname Directory for a dedicated virtual environment.
#' @param python Python executable used to create it; Python >= 3.10 is required.
#' @param install Whether to create the environment and install dependencies.
#' @param upstream_ref Official abc_atlas_access Git ref or commit.
#' @return Absolute Python executable path, invisibly.
#' @export
setup_environment <- function(envname = "~/.virtualenvs/importABCatlas", python = NULL,
                              install = FALSE, upstream_ref = "main") {
  envname <- normalizePath(path.expand(envname), mustWork = FALSE, winslash = "/")
  exe <- file.path(envname, if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python")
  run <- function(executable, args) {
    status <- system2(executable, vapply(args, shQuote, character(1)))
    if (status != 0L) stop("Environment command failed; inspect the output above.", call. = FALSE)
  }
  if (install) {
    if (is.null(python)) python <- .abc_python()
    ver <- system2(python, c("-c", shQuote("import sys; assert sys.version_info >= (3,10), 'Python >=3.10 required'")))
    if (ver != 0L) stop("Python >= 3.10 is required.", call. = FALSE)
    if (!file.exists(exe)) run(python, c("-m", "venv", envname))
    if (!grepl("^[A-Za-z0-9._/-]+$", upstream_ref)) stop("Invalid upstream_ref.", call. = FALSE)
    run(exe, c("-m", "pip", "install", "numpy>=1.24", "pandas>=2.0", "scipy>=1.10", "h5py>=3.8",
               paste0("git+https://github.com/AllenInstitute/abc_atlas_access.git@", upstream_ref)))
  }
  if (!file.exists(exe)) stop("Environment not found. Run setup_environment(install=TRUE) with internet access.", call. = FALSE)
  options(importABCatlas.python = path.expand(exe))
  message("Python selected for this R session: ", getOption("importABCatlas.python"))
  invisible(getOption("importABCatlas.python"))
}

#' Report the R and Python software environment
#'
#' This diagnostic deliberately does not invoke the ABC backend, so it remains
#' useful when optional Python packages such as h5py are missing.
#' @return A list of versions, interpreter details and dependency probes.
#' @export
diagnose <- function() {
  py <- tryCatch(.abc_python(), error = function(e) NA_character_)
  py_info <- list(executable = py, available = !is.na(py), version = NA_character_,
                  dependencies = character(), dependency_status = NA_integer_)
  if (!is.na(py)) {
    py_info$version <- tryCatch(system2(py, c("-c", shQuote("import sys; print(sys.version.replace('\\n',' '))")),
                                      stdout = TRUE, stderr = TRUE),
                                error = function(e) conditionMessage(e))
    probe <- paste(
      "import importlib.util",
      "mods=['numpy','pandas','scipy','h5py','abc_atlas_access']",
      "print('\\n'.join(m + '=' + ('OK' if importlib.util.find_spec(m) else 'MISSING') for m in mods))",
      sep = ";")
    dep <- tryCatch(system2(py, c("-c", shQuote(probe)), stdout = TRUE, stderr = TRUE),
                    error = function(e) structure(conditionMessage(e), status = 1L))
    py_info$dependencies <- unname(dep)
    py_info$dependency_status <- attr(dep, "status") %||% 0L
  }
  list(R = R.version.string, platform = R.version$platform,
       package = as.character(utils::packageVersion("importABCatlas")),
       python = py_info, session = utils::capture.output(utils::sessionInfo()))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
