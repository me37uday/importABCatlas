#' Resolve the ABC Atlas download/cache directory
#'
#' Determines where ABC Atlas data is downloaded and cached. The location is
#' resolved with the following precedence (first non-empty wins):
#' \enumerate{
#'   \item the explicit \code{path} argument;
#'   \item \code{getOption("importABCatlas.download_base")};
#'   \item the \code{ABC_ATLAS_CACHE} environment variable;
#'   \item a per-user cache directory, \code{tools::R_user_dir("importABCatlas", "cache")}.
#' }
#' The resolved directory is created if it does not yet exist and returned as an
#' absolute path. This lets users point downloads at a custom folder without
#' depending on the working directory.
#'
#' @param path Optional explicit directory. If \code{NULL} or empty, the other
#'   sources above are consulted in order.
#' @return An absolute path to the (existing) download directory.
#' @examples
#' \dontrun{
#' # Use a custom folder for this call only
#' load_data("WHB", download_base = "/data/abc_atlas")
#'
#' # Or set it once for the whole session
#' options(importABCatlas.download_base = "/data/abc_atlas")
#'
#' # Or from outside R, before launching:  export ABC_ATLAS_CACHE=/data/abc_atlas
#'
#' abc_download_dir()  # see where data will go
#' }
#' @export
abc_download_dir <- function(path = NULL) {
  is_set <- function(x) !is.null(x) && length(x) == 1L && !is.na(x) && nzchar(x)

  dir <- NULL
  if (is_set(path)) {
    dir <- path
  } else if (is_set(getOption("importABCatlas.download_base"))) {
    dir <- getOption("importABCatlas.download_base")
  } else {
    env <- Sys.getenv("ABC_ATLAS_CACHE", unset = "")
    if (is_set(env)) {
      dir <- env
    } else {
      dir <- tools::R_user_dir("importABCatlas", which = "cache")
    }
  }

  dir <- path.expand(dir)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(dir, mustWork = FALSE)
}

#' Connect to the ABC Atlas project cache
#'
#' Creates an \code{AbcProjectCache} pointing at \code{download_base}. With
#' \code{offline = FALSE} (default) it connects to S3 and downloads any missing
#' files; with \code{offline = TRUE} it uses only data already present on disk
#' (no network), via the package's local-cache constructor. The Python
#' environment must already be initialised (see \code{\link{setup_environment}}).
#'
#' @param download_base Download/cache directory; resolved by
#'   \code{\link{abc_download_dir}}.
#' @param offline If \code{TRUE}, use only locally cached data and never contact
#'   S3. Requires a manifest to have been downloaded previously.
#' @return A reticulate handle to an \code{AbcProjectCache} instance.
#' @export
abc_cache_connect <- function(download_base = NULL, offline = FALSE) {
  requireNamespace("reticulate")
  download_base <- abc_download_dir(download_base)

  # Validate offline preconditions before importing Python, so an offline call
  # with no local data fails fast without bootstrapping an interpreter.
  if (isTRUE(offline) && !.abc_manifest_exists(download_base)) {
    stop("offline = TRUE but no local manifest was found under '", download_base,
         "'. Run once with offline = FALSE to download the manifest and data.")
  }

  py_path <- reticulate::import("pathlib")$Path(download_base)
  AbcProjectCache <- reticulate::import(
    "abc_atlas_access.abc_atlas_cache.abc_project_cache")$AbcProjectCache

  if (isTRUE(offline)) {
    return(AbcProjectCache$from_local_cache(py_path))
  }
  AbcProjectCache$from_s3_cache(py_path)
}

# TRUE if a release manifest has already been downloaded into download_base.
.abc_manifest_exists <- function(download_base) {
  releases <- file.path(download_base, "releases")
  dir.exists(releases) &&
    length(list.files(releases, pattern = "manifest\\.json$", recursive = TRUE)) > 0
}
