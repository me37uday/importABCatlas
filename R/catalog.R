#' List curated dataset routes or live manifest directories
#' @param live Query the upstream manifest instead of the bundled registry.
#' @param download_base Cache directory, outside the Git repository.
#' @param release Pinned YYYYMMDD manifest release, or latest for exploration.
#' @param offline Use existing local files only.
#' @return A data frame describing each expression route, its modality, published representations, curated default, and release-scoped validation status. Registry membership alone is not a live validation certificate.
#' @export
list_datasets <- function(live = FALSE, download_base = "~/abc_atlas_cache",
                          release = .default_release, offline = FALSE) {
  if (live) {
    x <- list_files(download_base = download_base, release = release, offline = offline)
    return(unique(x["directory"]))
  }
  path <- system.file("python", "registry.json", package = "importABCatlas")
  registry <- jsonlite::read_json(path)
  do.call(rbind, lapply(names(registry), function(id) {
    s <- registry[[id]]
    available <- unlist(s$available_data_types %||% c("raw", "log2"), use.names = FALSE)
    validation_release <- s$validation_release %||% NA_character_
    requested_release <- sub("^releases/([0-9]+)/manifest\\.json$", "\\1", as.character(release))
    validation <- s$validation_status %||% "NOT_LIVE_VALIDATED"
    if (identical(tolower(requested_release), "latest")) {
      validation <- "DYNAMIC_RELEASE_NOT_VALIDATED"
    } else if (!is.na(validation_release) && !identical(requested_release, validation_release)) {
      validation <- sprintf("NOT_VALIDATED_FOR_RELEASE_%s", requested_release)
    }
    data.frame(
      dataset = id,
      directory = s$directory,
      gene_directory = s$gene_directory,
      taxonomy = if (is.null(s$taxonomy)) NA_character_ else s$taxonomy,
      modality = s$modality %||% if (isTRUE(s$spatial)) "spatial_transcriptomics" else "transcriptomics",
      spatial = isTRUE(s$spatial),
      imputed = isTRUE(s$imputed),
      available_data_types = paste(available, collapse = ","),
      default_data_type = s$default_data_type %||% available[[1L]],
      expected_measured_genes = as.integer(s$expected_measured_genes %||% NA_integer_),
      validation_release = validation_release,
      validation = validation,
      stringsAsFactors = FALSE
    )
  }))
}

#' Discover all metadata, expression, image and MapMyCells files in a manifest
#' @param directory Exact manifest directory or NULL for all directories.
#' @inheritParams list_datasets
#' @return A data frame with directory, kind and file_name plus manifest attribute.
#' @export
list_files <- function(directory = NULL, download_base = "~/abc_atlas_cache",
                       release = .default_release, offline = FALSE) {
  .abc_backend("catalog", c(.cache_args(download_base, release, offline), list(directory = directory)),
    function(out, result) {
      rows <- result$files
      x <- if (!length(rows)) data.frame(directory = character(), kind = character(), file_name = character()) else
        do.call(rbind, lapply(rows, function(z) as.data.frame(z, stringsAsFactors = FALSE)))
      attr(x, "manifest") <- result$manifest
      x
    })
}

#' Read an arbitrary manifest metadata table without forcing cell-level joins
#' @param directory Exact manifest directory.
#' @param file_name Exact metadata file name from list_files.
#' @param preview_n Optional number of first rows to read. Download is still whole-file.
#' @inheritParams list_datasets
#' @return A data frame with exact identifiers and a manifest attribute.
#' @export
get_metadata <- function(directory, file_name, download_base = "~/abc_atlas_cache",
                         release = .default_release, offline = FALSE, preview_n = NULL) {
  .abc_backend("metadata", c(.cache_args(download_base, release, offline),
                            list(directory = directory, file_name = file_name, preview_n = preview_n)),
    function(out, result) {
      x <- .read_frame(file.path(out, "table.csv")); attr(x, "manifest") <- result$manifest; x
    })
}

#' Download any manifest asset, including images and geometry files
#' @param allow_downloads Explicit consent to download a possibly very large file.
#' @inheritParams get_metadata
#' @return A local path with manifest attribute. It is not converted to an RNA assay.
#' @export
download_atlas_file <- function(directory, file_name, download_base = "~/abc_atlas_cache",
                                release = .default_release, offline = FALSE, allow_downloads = FALSE) {
  if (!offline && !allow_downloads) stop("Set allow_downloads=TRUE after checking size and disk space, or use offline=TRUE.", call. = FALSE)
  .abc_backend("download_file", c(.cache_args(download_base, release, offline),
                                  list(directory = directory, file_name = file_name)),
    function(out, result) structure(result$path, manifest = result$manifest))
}
