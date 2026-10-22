.prepare_fetch <- function(data, metadata, gene_data, dataset, download_base, release, overrides) {
  if (!is.null(data)) {
    if (!inherits(data, "abc_data")) stop("data must be returned by load_data().", call. = FALSE)
    metadata <- data$cell_metadata; gene_data <- data$gene_data
    dataset <- data$dataset; download_base <- data$download_base; release <- data$release
    overrides <- data$overrides
  }
  if (!is.data.frame(metadata) || !is.data.frame(gene_data)) stop("Provide data=load_data(...) or both metadata and gene_data.", call. = FALSE)
  if (!is.character(metadata$cell_label) || !is.character(gene_data$gene_identifier)) stop("cell_label and gene_identifier must be character vectors.", call. = FALSE)
  if (is.null(dataset)) dataset <- attr(metadata, "abc_dataset")
  if (is.null(dataset)) stop("Provide dataset explicitly; it cannot safely be inferred from gene names.", call. = FALSE)
  known <- attr(metadata, "abc_manifest")
  if (is.null(release)) release <- if (is.null(known)) .default_release else known
  if (!is.null(known) && !identical(sub("^releases/([0-9]+)/manifest.json$", "\\1", release),
                                  sub("^releases/([0-9]+)/manifest.json$", "\\1", known)))
    stop("Metadata and fetch manifest differ. Reload metadata for the intended release.", call. = FALSE)
  if (!length(overrides)) overrides <- attr(metadata, "abc_overrides")
  default_data_type <- if (!is.null(data)) {
    data$provenance$spec$default_data_type %||% "raw"
  } else {
    attr(metadata, "abc_default_data_type") %||% "raw"
  }
  available_data_types <- if (!is.null(data)) {
    unlist(data$provenance$spec$available_data_types %||% c("raw", "log2"),
           use.names = FALSE)
  } else {
    unlist(attr(metadata, "abc_available_data_types") %||% c("raw", "log2"),
           use.names = FALSE)
  }
  list(metadata = metadata, gene_data = gene_data, dataset = dataset,
       download_base = download_base, release = release, overrides = overrides,
       default_data_type = default_data_type,
       available_data_types = available_data_types)
}

.resolve_data_type <- function(data_type, prepared) {
  available <- as.character(prepared$available_data_types %||% c("raw", "log2"))
  if (!length(available)) {
    stop("This dataset route declares no accessible expression representation.",
         call. = FALSE)
  }
  if (is.null(data_type) || !length(data_type) || !nzchar(as.character(data_type[[1L]]))) {
    data_type <- prepared$default_data_type %||% available[[1L]]
  }
  value <- as.character(data_type[[1L]])
  if (!value %in% c("raw", "log2")) {
    stop("data_type must be one of: 'raw', 'log2'.", call. = FALSE)
  }
  if (!value %in% available) {
    stop(sprintf(
      "%s does not publish data_type='%s'. Available data_type values: %s; default: %s.",
      prepared$dataset, value, paste(sprintf("'%s'", available), collapse = ", "),
      prepared$default_data_type %||% available[[1L]]
    ), call. = FALSE)
  }
  value
}

.fetch_call <- function(action, p, args, reader) {
  dir <- tempfile("abc-tables-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  mp <- file.path(dir, "cells.csv"); gp <- file.path(dir, "genes.csv")
  utils::write.csv(p$metadata, mp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  utils::write.csv(p$gene_data, gp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  .abc_backend(action, c(.cache_args(p$download_base, p$release, args$offline),
                        list(dataset = p$dataset, overrides = p$overrides,
                             metadata_path = mp, gene_path = gp), args[names(args) != "offline"]), reader)
}

#' Preview expression-file requirements without downloading expression
#' @param data An abc_data object; preferred over separate metadata arguments.
#' @param metadata Cell metadata, for compatibility with the earlier API.
#' @param gene_data Gene metadata, for compatibility with the earlier API.
#' @param genes Exact gene identifiers or unambiguous symbols. `NULL` reconciles the authoritative H5AD feature axis, keeps at most one metadata row per measured feature in every selected source, and records unmatched/collapsed alias rows in provenance. Explicit unavailable or duplicate-feature requests remain errors and are never replaced by zeros or aggregated.
#' @param data_type Published representation for the selected route (`raw` or `log2`). NULL uses the route's declared default. `MERFISH_imputed` and `HMBA_PatchSeq` are log2-only in the validated release; unavailable representations fail before any expression download. Raw is not automatically proof of UMI counts.
#' @param file_map Optional named list mapping feature_matrix_label to directory/file_name lists.
#' @inheritParams load_data
#' @return An abc_fetch_plan with exact selected-file sizes when available, cache status, additional download bytes, disk-space status, and selected output dimensions. Directory totals are never reported as file sizes.
#' @export
plan_fetch <- function(data = NULL, metadata = NULL, gene_data = NULL, genes = NULL,
                       dataset = NULL, download_base = "~/abc_atlas_cache", release = NULL,
                       offline = FALSE, filters = list(), cell_ids = NULL, n_cells = NULL,
                       seed = 1L, data_type = NULL, file_map = NULL, overrides = list()) {
  p <- .prepare_fetch(data, metadata, gene_data, dataset, download_base, release, overrides)
  data_type <- .resolve_data_type(data_type, p)
  ans <- .fetch_call("plan", p, list(offline = offline, filters = filters, cell_ids = cell_ids,
              n_cells = n_cells, seed = seed, genes = genes, data_type = data_type, file_map = file_map),
              function(out, result) result)
  class(ans) <- c("abc_fetch_plan", class(ans))
  ans
}

#' Read a sparse expression subset and construct an analysis object
#'
#' Complete source H5AD files may be downloaded even for a very small output.
#' Subsets are read from X, not a guessed layer; cell and gene matching is exact.
#' H5AD feature axes may use stable identifiers or symbols. The importer performs
#' deterministic alias reconciliation against AnnData `var` metadata and records
#' the match methods in provenance.
#' Matrix return type is a list containing expression, metadata, features and provenance.
#' Imputed or transformed data is never placed in a counts assay.
#' @param assay_name Seurat assay name.
#' @param return_type seurat, sce, matrix, or spatial (SpatialExperiment).
#' @param allow_downloads Explicit consent to expression downloads.
#' @param max_output_mb Dense-equivalent output-size guard, not a total RAM guarantee.
#' @param obs_column Metadata column whose values exactly match the H5AD obs index.
#' @param counts_semantics NULL infers integer nonnegative raw values; FALSE prevents a counts label.
#' @param verify_reference Compare each selected row to an independent AnnData read (small live tests).
#' @param coordinate_columns Explicit anatomical coordinate columns for spatial output.
#' @inheritParams plan_fetch
#' @return A Seurat, SingleCellExperiment, SpatialExperiment, or sparse-matrix list.
#' @export
fetch_data <- function(download_base = "~/abc_atlas_cache", metadata = NULL, filters = list(),
                        gene_data = NULL, genes = NULL, assay_name = "RNA",
                        return_type = c("seurat", "sce", "matrix", "spatial"), data = NULL,
                        dataset = NULL, release = NULL, offline = FALSE, allow_downloads = FALSE,
                        data_type = NULL, max_output_mb = 1024, cell_ids = NULL, n_cells = NULL,
                        seed = 1L, file_map = NULL, obs_column = NULL, counts_semantics = NULL,
                        coordinate_columns = NULL, overrides = list(), verify_reference = FALSE) {
  return_type <- match.arg(return_type)
  p <- .prepare_fetch(data, metadata, gene_data, dataset, download_base, release, overrides)
  data_type <- .resolve_data_type(data_type, p)
  .fetch_call("fetch", p,
    list(offline = offline, filters = filters, genes = genes, data_type = data_type,
         allow_downloads = allow_downloads, max_output_mb = max_output_mb, cell_ids = cell_ids,
         n_cells = n_cells, seed = seed, file_map = file_map, obs_column = obs_column,
         counts_semantics = counts_semantics, verify_reference = verify_reference),
    function(out, result) {
      result$dataset <- p$dataset
      if (!is.null(data)) result$metadata_provenance <- data$provenance
      .materialize(out, result, return_type, assay_name, coordinate_columns,
                   spatial = if (is.null(data)) FALSE else isTRUE(data$provenance$spec$spatial))
    })
}

.materialize <- function(out, provenance, return_type, assay_name, coordinate_columns = NULL, spatial = FALSE) {
  cm <- .read_frame(file.path(out, "cell_metadata.csv"))
  gd <- .read_frame(file.path(out, "gene_data.csv"))
  m <- methods::as(Matrix::readMM(file.path(out, "matrix.mtx")), "CsparseMatrix")
  stopifnot(nrow(m) == nrow(gd), ncol(m) == nrow(cm), !anyDuplicated(cm$cell_label), !anyDuplicated(gd$gene_identifier))
  rownames(cm) <- cm$cell_label; rownames(gd) <- gd$gene_identifier
  dimnames(m) <- list(gd$gene_identifier, cm$cell_label)
  kind <- provenance$assay_kind
  if (return_type == "matrix") return(list(expression = m, cell_metadata = cm, gene_data = gd, provenance = provenance))
  if (return_type == "seurat") {
    .need_package("SeuratObject")
    if (any(grepl("[_|]", rownames(m)))) stop("Seurat would alter these feature IDs; choose return_type='sce' or 'matrix' to preserve them exactly.", call. = FALSE)
    if (kind == "counts") {
      object <- SeuratObject::CreateSeuratObject(counts = m, assay = assay_name, meta.data = cm,
                                                min.cells = 0, min.features = 0)
    } else {
      if (utils::packageVersion("SeuratObject") < "5.0.0") stop("Transformed assays require SeuratObject >=5.", call. = FALSE)
      if (kind %in% c("log2", "imputed_log2"))
        warning("This is upstream log2-scale data, not necessarily Seurat's usual natural-log normalization. Do not apply count-based normalization or invert it with expm1 without checking the source transformation.", call. = FALSE)
      assay <- SeuratObject::CreateAssay5Object(data = m)
      object <- SeuratObject::CreateSeuratObject(counts = assay, assay = assay_name, meta.data = cm)
    }
    object[[assay_name]] <- SeuratObject::AddMetaData(object[[assay_name]], metadata = gd)
    object@misc$importABCatlas <- provenance
    return(object)
  }
  .need_package("SingleCellExperiment"); .need_package("S4Vectors")
  assays <- stats::setNames(list(m), kind)
  if (return_type == "sce") return(SingleCellExperiment::SingleCellExperiment(
    assays = assays, colData = S4Vectors::DataFrame(cm), rowData = S4Vectors::DataFrame(gd),
    metadata = list(importABCatlas = provenance)))
  .need_package("SpatialExperiment")
  if (is.null(coordinate_columns)) {
    if (!spatial) stop("No declared spatial modality. Supply verified anatomical coordinate_columns explicitly; UMAP is not spatial location.", call. = FALSE)
    coordinate_columns <- if (all(c("x_slab_mm", "y_slab_mm") %in% names(cm))) c("x_slab_mm", "y_slab_mm") else c("x", "y")
  }
  if (length(coordinate_columns) < 2L || !all(coordinate_columns %in% names(cm)))
    stop("Anatomical coordinates absent. Inspect metadata or supply coordinate_columns.", call. = FALSE)
  if (!all(vapply(cm[coordinate_columns], is.numeric, logical(1)))) stop("Spatial coordinates must be numeric.", call. = FALSE)
  xy <- as.matrix(cm[coordinate_columns])
  if (any(!is.finite(xy))) stop("Missing spatial coordinates: filter deliberately before constructing SpatialExperiment.", call. = FALSE)
  provenance$coordinate_columns <- coordinate_columns
  SpatialExperiment::SpatialExperiment(assays = assays, colData = S4Vectors::DataFrame(cm),
    rowData = S4Vectors::DataFrame(gd), spatialCoords = xy, metadata = list(importABCatlas = provenance))
}

#' Build an analysis object from a tiny bundled synthetic H5AD fixture
#' @inheritParams fetch_data
#' @return An object in the requested format; never biological benchmark evidence.
#' @export
example_data <- function(return_type = c("matrix", "seurat", "sce", "spatial")) {
  return_type <- match.arg(return_type)
  .abc_backend("demo", reader = function(out, result)
    .materialize(out, result, return_type, "RNA", coordinate_columns = c("x", "y"), spatial = TRUE))
}


#' @export
print.abc_fetch_plan <- function(x, ...) {
  fmt_bytes <- function(n) {
    if (is.null(n) || length(n) == 0L || is.na(n)) return("unknown")
    n <- as.numeric(n)
    units <- c("B", "KiB", "MiB", "GiB", "TiB")
    i <- 1L
    value <- n
    while (value >= 1024 && i < length(units)) {
      value <- value / 1024
      i <- i + 1L
    }
    sprintf("%.2f %s (%s bytes)", value, units[[i]],
            format(n, scientific = FALSE, big.mark = ",", trim = TRUE))
  }

  cat("importABCatlas fetch plan

")
  cat("Expression representation:", x$data_type %||% "unknown", "
")
  cat("Requested output:", x$n_genes, "genes x", x$n_cells, "cells
")
  if (!is.null(x$gene_count_basis)) {
    cat("Gene-count basis:", x$gene_count_basis, "
")
  }
  cat("Expression files that would be used:
")
  fs <- x$files
  if (length(fs)) {
    for (f in fs) {
      cat(" - ", f$directory, "/", f$file_name, " [", f$n_cells, " cells]
", sep = "")
      if (!is.null(f$relative_path)) cat("   S3/cache object: ", f$relative_path, "
", sep = "")
      cat("   exact selected-file size: ", fmt_bytes(f$file_bytes),
          " | source: ", f$size_source %||% "unknown",
          " | cached: ", if (is.null(f$cached)) "unknown" else as.character(f$cached),
          " | resumable partial: ", fmt_bytes(f$partial_bytes),
          " | additional download: ", fmt_bytes(f$download_required_bytes), "
", sep = "")
      if (!is.null(f$size_error) && nzchar(f$size_error))
        cat("   size lookup note: ", f$size_error, "
", sep = "")
    }
  }
  cat("Total selected-file bytes: ", fmt_bytes(x$total_file_bytes), "
", sep = "")
  cat("Total additional download: ", fmt_bytes(x$total_download_bytes), "
", sep = "")
  if (!is.null(x$disk)) {
    cat("Free space in cache filesystem: ", fmt_bytes(x$disk$free_bytes), "
", sep = "")
    cat("Required free space including safety margin/reserve: ",
        fmt_bytes(x$disk$required_free_bytes), "
", sep = "")
    cat("Enough space: ",
        if (is.null(x$disk$enough_space)) "unknown" else as.character(x$disk$enough_space), "
", sep = "")
  }
  if (!is.null(x$warning) && nzchar(x$warning)) cat("
", x$warning, "
", sep = "")
  invisible(x)
}
