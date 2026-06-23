#' Set up the Python environment for importABCatlas
#'
#' Declares the Python packages the loaders depend on via
#' \code{reticulate::py_require()} and initialises Python so any provisioning
#' errors surface early. reticulate uses \code{uv} to fetch a suitable prebuilt
#' Python interpreter (no system Python or compilation required) and resolve the
#' declared dependencies into an isolated environment.
#'
#' Calling this is optional: the dependencies are also declared automatically
#' when the package is loaded (see \code{.onLoad}), and Python is initialised
#' lazily on first use. Use it to provision and check the environment up front.
#'
#' @return Invisibly, the reticulate Python configuration.
#' @export
setup_environment <- function() {
  requireNamespace("reticulate")
  reticulate::py_require(.abc_python_packages())
  cfg <- reticulate::py_config()   # triggers provisioning / initialisation
  message("Python ready: ", cfg$python)
  invisible(cfg)
}

# Declarative Python dependencies, in PEP 508 form for reticulate/uv.
#
# abc_atlas_access is pinned to a specific commit for reproducibility (it has no
# stable PyPI release and otherwise tracks a moving branch). The scientific
# stack is left unpinned so uv can resolve versions compatible with whatever
# prebuilt Python it selects -- keeping setup portable across machines. numpy,
# scipy and anndata are listed explicitly because fetch_data()'s helper imports
# them directly (not only transitively through abc_atlas_access).
.abc_python_packages <- function() {
  c(
    "abc_atlas_access @ git+https://github.com/alleninstitute/abc_atlas_access@da54664e8a9c8e2e2ba573d1d8f590b114a149c4",
    "anndata",
    "numpy",
    "pandas",
    "scipy"
  )
}
