# =============================================================================
# ABC Atlas cache: release manifests, file downloads and metadata tables.
#
# Re-implements the parts of abc_atlas_access's AbcProjectCache that the package
# uses, with the same on-disk layout so a download directory can be shared with
# the Python package:
#   <download_base>/releases/<version>/manifest.json   release manifests
#   <download_base>/<relative_path>                    data files, as in the manifest
#   <download_base>/_manifest_last_used.txt            manifest pinned for this cache
#   <download_base>/_downloaded_data.json              files fully downloaded and verified
# =============================================================================

# Public S3 bucket holding every ABC Atlas release.
.abc_bucket_url <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

# Manifest versions the loaders understand (same window as abc_atlas_access).
.abc_manifest_compatibility <- c("20230101", "20300101")

# Strings pandas.read_csv treats as missing by default; reused so metadata tables
# have the same missing values as with the Python package.
.abc_na_strings <- c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN",
  "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN",
  "None", "n/a", "nan", "null")

#' Connect to the ABC Atlas project cache
#'
#' Selects a release manifest and returns a handle used to locate, download and
#' read ABC Atlas files under \code{download_base}. With \code{offline = FALSE}
#' (default) the available releases are listed from S3 and missing files are
#' downloaded (and md5-verified) on first use; with \code{offline = TRUE} only
#' manifests and files already on disk are used, and the network is never
#' contacted.
#'
#' The manifest is chosen as in the Python \code{abc_atlas_access} package: the
#' one last used with this directory if any (with a warning when a newer release
#' exists), otherwise the latest release. The directory layout is also the same,
#' so a download directory can be shared with the Python package.
#'
#' @param download_base Download/cache directory; resolved by
#'   \code{\link{abc_download_dir}}.
#' @param offline If \code{TRUE}, use only locally cached data and never contact
#'   S3. Requires a manifest to have been downloaded previously.
#' @return An \code{abc_cache} object: a list with the \code{download_base},
#'   \code{manifest_name}, \code{version} and \code{offline} fields, plus the
#'   functions \code{get_file_path(directory, file_name)} (returns the local path
#'   of a file, downloading it if needed) and
#'   \code{get_metadata_dataframe(directory, file_name, col_classes = NULL, ...)}
#'   (reads a metadata CSV into a data.frame; \code{...} goes to
#'   \code{data.table::fread}).
#' @examples
#' \dontrun{
#' cache <- abc_cache_connect()
#' donors <- cache$get_metadata_dataframe("ASAP-PMDBS-10X", "donor")
#' }
#' @export
abc_cache_connect <- function(download_base = NULL, offline = FALSE) {
  download_base <- abc_download_dir(download_base)
  offline <- isTRUE(offline)

  # Validate offline preconditions up front so an offline call with no local
  # data fails fast with an actionable message.
  if (offline && !.abc_manifest_exists(download_base)) {
    stop("offline = TRUE but no local manifest was found under '", download_base,
         "'. Run once with offline = FALSE to download the manifest and data.")
  }

  available <- if (offline) .abc_local_manifests(download_base) else .abc_remote_manifests()
  if (length(available) == 0L) stop("No ABC Atlas release manifests found.")
  manifest_name <- .abc_choose_manifest(download_base, available)

  manifest_path <- file.path(download_base, manifest_name)
  if (!file.exists(manifest_path)) {
    .abc_download(paste0(.abc_bucket_url, "/", manifest_name), manifest_path)
  }
  manifest <- jsonlite::read_json(manifest_path)
  .abc_check_manifest_version(manifest$version)

  # Pin this manifest for later sessions. Written without a trailing newline, as
  # abc_atlas_access compares the file content verbatim.
  tryCatch(
    cat(manifest_name, file = file.path(download_base, "_manifest_last_used.txt")),
    error = function(e) warning("Could not record the manifest in use: ",
                                conditionMessage(e), call. = FALSE)
  )

  cache <- list(
    download_base = download_base,
    manifest_name = manifest_name,
    version = manifest$version,
    offline = offline,
    manifest = manifest
  )
  cache$get_file_path <- function(directory, file_name, force_download = FALSE,
                                  skip_hash_check = FALSE) {
    .abc_get_file_path(cache, directory, file_name, force_download, skip_hash_check)
  }
  cache$get_metadata_dataframe <- function(directory, file_name, col_classes = NULL,
                                           ...) {
    .abc_read_csv(cache$get_file_path(directory, file_name), col_classes, ...)
  }
  structure(cache, class = "abc_cache")
}

#' @export
print.abc_cache <- function(x, ...) {
  cat("<abc_cache> ", x$manifest_name, if (x$offline) " (offline)", "\n",
      "  download_base: ", x$download_base, "\n",
      "  directories:   ", length(x$manifest$file_listing), "\n", sep = "")
  invisible(x)
}

# ---- Manifests --------------------------------------------------------------

# TRUE if a release manifest has already been downloaded into download_base.
.abc_manifest_exists <- function(download_base) {
  length(.abc_local_manifests(download_base)) > 0L
}

# Downloaded manifests as sorted "releases/<version>/manifest.json" names.
.abc_local_manifests <- function(download_base) {
  found <- Sys.glob(file.path(download_base, "releases", "*", "manifest.json"))
  if (length(found) == 0L) return(character())
  sort(paste0("releases/", basename(dirname(found)), "/manifest.json"),
       method = "radix")
}

# Manifests published on S3 (every .json key under releases/), sorted.
.abc_remote_manifests <- function() {
  keys <- character()
  token <- NULL
  repeat {
    url <- paste0(.abc_bucket_url, "/?list-type=2&prefix=releases/")
    if (!is.null(token)) {
      url <- paste0(url, "&continuation-token=", curl::curl_escape(token))
    }
    resp <- curl::curl_fetch_memory(url)
    if (resp$status_code != 200L) {
      stop("Could not list ABC Atlas releases (HTTP ", resp$status_code, ").")
    }
    doc <- xml2::xml_ns_strip(xml2::read_xml(resp$content))
    keys <- c(keys, xml2::xml_text(xml2::xml_find_all(doc, "//Contents/Key")))
    if (xml2::xml_text(xml2::xml_find_first(doc, "//IsTruncated")) != "true") break
    token <- xml2::xml_text(xml2::xml_find_first(doc, "//NextContinuationToken"))
  }
  sort(keys[endsWith(keys, ".json")], method = "radix")
}

# Pick the manifest to load: the one last used with this cache if still valid,
# otherwise the latest available. Warns as abc_atlas_access does when that is
# not the newest release, or when it is newer than what was used before.
.abc_choose_manifest <- function(download_base, available) {
  latest <- utils::tail(available, 1L)
  last_used_file <- file.path(download_base, "_manifest_last_used.txt")

  if (file.exists(last_used_file)) {
    to_load <- paste(readLines(last_used_file, warn = FALSE), collapse = "\n")
    if (!to_load %in% available) {
      warning("The manifest recorded as last used for this cache -- ", to_load,
              " -- is not a valid manifest. Loading the latest version -- ",
              latest, " -- instead.", call. = FALSE)
      return(latest)
    }
    if (to_load != latest) {
      warning("You are loading ", to_load, ". A more up to date version of the ",
              "dataset -- ", latest, " -- exists. To use it, delete '",
              last_used_file, "' or point download_base to a new directory.",
              call. = FALSE)
    }
    return(to_load)
  }

  downloaded <- .abc_local_manifests(download_base)
  if (length(downloaded) > 0L && utils::tail(downloaded, 1L) != latest) {
    warning("You are loading ", latest, ", which is newer than the most recent ",
            "manifest previously used here (", utils::tail(downloaded, 1L), "). ",
            "Some data files may have changed between the two releases and ",
            "will be downloaded again (existing files are not overwritten).",
            call. = FALSE)
  }
  latest
}

.abc_check_manifest_version <- function(version) {
  rng <- .abc_manifest_compatibility
  if (is.null(version) || version < rng[1] || version >= rng[2]) {
    stop("The manifest has version ", version, " but this package supports ",
         "versions on or after ", rng[1], " and before ", rng[2], ".")
  }
}

# ---- Files ------------------------------------------------------------------

# Resolve a file in the manifest. `file_name` is either a file name or
# "<name>/<kind>" for entries with several kinds (e.g. "WHB-10Xv3-Neurons/raw").
.abc_file_attributes <- function(cache, directory, file_name) {
  parts <- strsplit(file_name, "/", fixed = TRUE)[[1]]
  name <- parts[1]
  kind <- if (length(parts) > 1L) parts[2] else NULL

  dir_data <- cache$manifest$file_listing[[directory]]
  if (is.null(dir_data)) {
    stop("Directory '", directory, "' not found in manifest ", cache$manifest_name, ".")
  }

  entry <- NULL
  for (sub_dir in names(dir_data)) {
    files_data <- dir_data[[sub_dir]][[name]]
    if (is.null(files_data)) next
    if (!is.null(files_data[["files"]])) {
      entry <- files_data[["files"]]
    } else if (!is.null(kind) && !is.null(files_data[[kind]])) {
      entry <- files_data[[kind]][["files"]]
    } else if (is.null(kind)) {
      stop("File '", name, "' has several kinds in directory '", directory,
           "'. Specify the file name as one of: ",
           paste0(name, "/", names(files_data), collapse = ", "))
    }
  }
  if (is.null(entry)) {
    stop("File '", name, "' not found in directory '", directory, "'.")
  }

  f <- entry[[1]]
  list(url = f$url, version = f$version, size = f$size,
       relative_path = f$relative_path, file_type = names(entry)[1],
       file_hash = f$file_hash,
       local_path = file.path(cache$download_base, f$relative_path))
}

# Local path of a manifest file, downloading it first if needed (online only).
# A file on disk counts as present only if it is recorded as fully downloaded;
# otherwise it is fetched again, as abc_atlas_access does.
.abc_get_file_path <- function(cache, directory, file_name, force_download = FALSE,
                               skip_hash_check = FALSE) {
  attrs <- .abc_file_attributes(cache, directory, file_name)
  path <- attrs$local_path
  if (dir.exists(path)) stop(path, " exists but is not a file.")

  if (cache$offline) {
    if (!file.exists(path)) {
      stop("'", attrs$relative_path, "' is not available locally (offline = TRUE). ",
           "Run once with offline = FALSE to download it.")
    }
    return(path)
  }

  base <- cache$download_base
  if (file.exists(path) && !.abc_is_recorded(base, path)) force_download <- TRUE
  if (force_download) {
    .abc_unrecord(base, path)
    unlink(path)
  }
  if (!file.exists(path)) {
    message("Downloading ", attrs$relative_path, " (", .abc_format_size(attrs$size), ")")
    .abc_download(attrs$url, path, if (!skip_hash_check) attrs$file_hash)
    .abc_record(base, path)
  }
  path
}

# Download `url` to `dest` through a ".part" file that is resumed across retries
# and only moved into place once complete (and matching `md5`, if given).
.abc_download <- function(url, dest, md5 = NULL, max_tries = 10L) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  part <- paste0(dest, ".part")
  for (i in seq_len(max_tries)) {
    res <- curl::multi_download(url, part, resume = TRUE, progress = interactive())
    ok <- isTRUE(res$success)
    if (ok && !is.null(md5) && unname(tools::md5sum(part)) != md5) {
      message("Checksum mismatch for ", basename(dest), "; downloading again.")
      unlink(part)
      ok <- FALSE
    }
    if (ok) {
      file.rename(part, dest)
      return(invisible(dest))
    }
  }
  stop("Could not download ", url, " in ", max_tries, " attempts",
       if (!is.null(res$error) && !is.na(res$error)) paste0(": ", res$error), ".")
}

.abc_format_size <- function(bytes) {
  format(structure(as.numeric(bytes), class = "object_size"), units = "auto",
         standard = "SI")
}

# ---- Record of completed downloads (_downloaded_data.json) -------------------
# A JSON array of absolute paths, shared with abc_atlas_access.

.abc_record_path <- function(base) file.path(base, "_downloaded_data.json")

.abc_recorded <- function(base) {
  p <- .abc_record_path(base)
  if (!file.exists(p)) return(character())
  as.character(unlist(jsonlite::read_json(p)))
}

.abc_write_record <- function(base, paths) {
  jsonlite::write_json(unique(paths), .abc_record_path(base), pretty = TRUE)
}

.abc_is_recorded <- function(base, path) {
  normalizePath(path) %in% .abc_recorded(base)
}

.abc_record <- function(base, path) {
  .abc_write_record(base, c(.abc_recorded(base), normalizePath(path)))
}

.abc_unrecord <- function(base, path) {
  rec <- .abc_recorded(base)
  abs <- normalizePath(path, mustWork = FALSE)
  if (abs %in% rec) .abc_write_record(base, setdiff(rec, abs))
}

# ---- Metadata tables --------------------------------------------------------

# Read a metadata CSV into a data.frame. `col_classes` is a named vector of
# column classes (e.g. c(cell_label = "character")); other types are inferred.
# Integers too large for R's integer type are read as doubles; surrounding
# whitespace is kept, as pandas does.
.abc_read_csv <- function(path, col_classes = NULL, ...) {
  data.table::fread(path, colClasses = col_classes, na.strings = .abc_na_strings,
                    integer64 = "double", strip.white = FALSE, data.table = FALSE,
                    showProgress = FALSE, ...)
}
