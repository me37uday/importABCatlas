args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Expected config_file, requirements_file, requirements_id, package_version.")
}
config_file <- path.expand(args[[1L]])
requirements_file <- path.expand(args[[2L]])
requirements_id <- args[[3L]]
package_version <- args[[4L]]

dir.create(dirname(config_file), recursive = TRUE, showWarnings = FALSE)

requirements <- trimws(readLines(requirements_file, warn = FALSE))
requirements <- requirements[nzchar(requirements) & !startsWith(requirements, "#")]
if (!length(requirements)) stop("Managed Python requirements are empty.")

# This runs in a dedicated Rscript process, so it cannot conflict with a Python
# interpreter already initialized by another package in the user's R session.
# Remove global Python path settings that can otherwise leak into a virtual
# environment and produce hard-to-diagnose import errors.
Sys.unsetenv(c("PYTHONHOME", "PYTHONPATH"))
Sys.setenv(PYTHONNOUSERSITE = "1")
Sys.setenv(
  RETICULATE_PYTHON = "managed",
  RETICULATE_USE_MANAGED_VENV = "yes"
)
# Keep the managed environment available for reproducible command-line backend
# calls. If a user clears the reticulate cache, importABCatlas recreates it.
options(reticulate.max_cache_age = as.difftime(3650, units = "days"))

reticulate::py_require(
  packages = requirements,
  python_version = "3.11",
  exclude_newer = "2026-10-20",
  action = "set"
)
config <- reticulate::py_config()
python <- path.expand(config$python)
if (!file.exists(python)) stop("reticulate returned a missing Python executable: ", python)

probe <- paste(
  "import importlib.util, json, sys",
  "mods=['numpy','pandas','scipy','h5py','anndata','boto3','botocore','tqdm','abc_atlas_access']",
  "missing=[m for m in mods if importlib.util.find_spec(m) is None]",
  "ok=(sys.version_info >= (3,10) and not missing)",
  "print(json.dumps({'ok':ok,'version':sys.version.split()[0],'missing':missing,'executable':sys.executable}))",
  "raise SystemExit(0 if ok else 2)",
  sep = "; "
)
out <- system2(
  python, c("-c", shQuote(probe)),
  stdout = TRUE, stderr = TRUE,
  env = c("PYTHONHOME=", "PYTHONPATH=", "PYTHONNOUSERSITE=1")
)
status <- attr(out, "status")
if (!is.null(status) && as.integer(status) != 0L) {
  stop("Managed Python dependency probe failed:\n", paste(out, collapse = "\n"))
}
parsed <- jsonlite::fromJSON(tail(out, 1L), simplifyVector = TRUE)
if (!isTRUE(parsed$ok)) stop("Managed Python is incomplete: ", paste(parsed$missing, collapse = ", "))

payload <- list(
  python = python,
  python_version = as.character(parsed$version),
  requirements_id = requirements_id,
  package_version = package_version,
  manager = "reticulate-uv",
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
tmp <- paste0(config_file, ".tmp")
jsonlite::write_json(payload, tmp, auto_unbox = TRUE, pretty = TRUE)
if (file.exists(config_file)) unlink(config_file, force = TRUE)
if (!file.rename(tmp, config_file)) {
  unlink(tmp, force = TRUE)
  stop("Could not save managed Python configuration: ", config_file)
}
cat(python, "\n", sep = "")
