.onLoad <- function(libname, pkgname) {
  # Declare the Python dependencies as soon as the package is loaded, before
  # Python is initialised, so reticulate (via uv) can provision a suitable
  # prebuilt interpreter and resolve them on first use. Wrapped in try() so a
  # package load never fails on this (e.g. if Python is already initialised).
  if (requireNamespace("reticulate", quietly = TRUE)) {
    try(reticulate::py_require(.abc_python_packages()), silent = TRUE)
  }
}
