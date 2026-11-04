# Resolve either an abc_data object or the legacy metadata/gene_data pair.
.abc_fetch_inputs <- function(data = NULL, metadata = NULL, gene_data = NULL,
                              download_base = NULL, release = NULL) {
  if (!is.null(data)) {
    if (!inherits(data, "abc_data")) {
      stop(
        "data must be an abc_data object returned by load_data().",
        call. = FALSE
      )
    }

    if (!is.null(metadata) || !is.null(gene_data)) {
      stop(
        "Supply either data or metadata/gene_data, not both.",
        call. = FALSE
      )
    }

    metadata <- data$cell_metadata
    gene_data <- data$gene_data

    if (is.null(download_base)) {
      download_base <- data$download_base
    }

    if (is.null(release)) {
      release <- data$release
    }
  }

  if (is.null(metadata) || is.null(gene_data)) {
    stop(
      "Supply an abc_data object via data, or both metadata and gene_data.",
      call. = FALSE
    )
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(gene_data)) {
    stop("gene_data must be a data.frame.", call. = FALSE)
  }

  list(
    metadata = metadata,
    gene_data = gene_data,
    download_base = download_base,
    release = release
  )
}



.abc_require_output_package <- function(package, return_type) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "return_type='", return_type,
      "' requires package '", package,
      "'. Install it to use this output format.",
      call. = FALSE
    )
  }
}

.abc_spatial_coordinates <- function(metadata, coordinate_columns = NULL) {
  if (!is.null(coordinate_columns)) {
    if (!is.character(coordinate_columns) ||
        length(coordinate_columns) != 2L ||
        anyNA(coordinate_columns) ||
        any(!nzchar(coordinate_columns))) {
      stop(
        "coordinate_columns must contain exactly two non-missing column names.",
        call. = FALSE
      )
    }

    missing <- setdiff(coordinate_columns, names(metadata))

    if (length(missing)) {
      stop(
        "Spatial coordinate column",
        if (length(missing) == 1L) " " else "s ",
        paste(missing, collapse = ", "),
        " not found in cell metadata.",
        call. = FALSE
      )
    }

    chosen <- coordinate_columns
  } else {
    candidates <- list(
      c("x_slab_mm", "y_slab_mm"),
      c("x_rotated_um", "y_rotated_um"),
      c("x_experiment", "y_experiment"),
      c("x", "y")
    )

    present <- vapply(
      candidates,
      function(x) all(x %in% names(metadata)),
      logical(1)
    )

    if (!any(present)) {
      stop(
        "No supported spatial coordinate pair was found. ",
        "Provide coordinate_columns explicitly.",
        call. = FALSE
      )
    }

    chosen <- candidates[[which(present)[[1L]]]]
  }

  coordinates <- as.matrix(
    metadata[, chosen, drop = FALSE]
  )

  storage.mode(coordinates) <- "double"

  if (anyNA(coordinates) || any(!is.finite(coordinates))) {
    stop(
      "Spatial coordinates must be finite and non-missing.",
      call. = FALSE
    )
  }

  rownames(coordinates) <- as.character(metadata$cell_label)
  colnames(coordinates) <- chosen

  coordinates
}

.abc_output_provenance <- function(metadata, gene_data, data_type,
                                   release = NULL) {
  list(
    dataset = attr(metadata, "abc_dataset"),
    release = release,
    data_type = data_type,
    modality = attr(metadata, "abc_modality"),
    spatial = isTRUE(attr(metadata, "abc_spatial")),
    expression_directories = attr(
      metadata,
      "abc_expression_directories"
    ),
    n_cells = nrow(metadata),
    n_genes = nrow(gene_data)
  )
}

.abc_muffle_seurat_layer_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (conditionMessage(w) %in% c(
        "No layers found matching search pattern provided",
        "Cannot find the layer(s) specified"
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

.abc_materialize <- function(expression, cell_metadata, gene_data,
                             data_type, return_type = "seurat",
                             assay_name = "RNA",
                             coordinate_columns = NULL,
                             provenance = NULL,
                             spatial = FALSE) {
  return_type <- match.arg(
    return_type,
    c("seurat", "sce", "matrix", "spatial")
  )

  if (!inherits(expression, "dgCMatrix")) {
    expression <- methods::as(
      methods::as(expression, "generalMatrix"),
      "CsparseMatrix"
    )
    expression <- methods::as(
      expression,
      "dgCMatrix"
    )
  }

  if (ncol(expression) != nrow(cell_metadata)) {
    stop(
      "Expression columns do not match the selected cell metadata.",
      call. = FALSE
    )
  }

  if (nrow(expression) != nrow(gene_data)) {
    stop(
      "Expression rows do not match the selected gene metadata.",
      call. = FALSE
    )
  }

  if (!identical(
    as.character(colnames(expression)),
    as.character(cell_metadata$cell_label)
  )) {
    stop(
      "Expression columns are not aligned to cell_metadata$cell_label.",
      call. = FALSE
    )
  }

  gene_ids <- as.character(gene_data$gene_identifier)

  if (anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids)) {
    stop(
      "gene_data$gene_identifier must be a non-missing unique key.",
      call. = FALSE
    )
  }

  if (!identical(
    as.character(rownames(expression)),
    gene_ids
  )) {
    stop(
      "Expression rows are not aligned to gene_data$gene_identifier.",
      call. = FALSE
    )
  }

  rownames(cell_metadata) <- as.character(
    cell_metadata$cell_label
  )
  rownames(gene_data) <- gene_ids

  if (is.null(provenance)) {
    provenance <- list(
      data_type = data_type,
      modality = attr(
        cell_metadata,
        "abc_modality"
      ),
      spatial = isTRUE(
        attr(cell_metadata, "abc_spatial")
      ),
      n_cells = nrow(cell_metadata),
      n_genes = nrow(gene_data)
    )
  }

  if (return_type == "matrix") {
    return(
      structure(
        list(
          expression = expression,
          cell_metadata = cell_metadata,
          gene_data = gene_data,
          data_type = data_type,
          provenance = provenance
        ),
        class = c("abc_matrix", "list")
      )
    )
  }

  if (return_type == "seurat") {
    .abc_require_output_package(
      "SeuratObject",
      "seurat"
    )

    if (any(grepl("[_|]", rownames(expression)))) {
      stop(
        "Seurat would alter these feature identifiers; ",
        "choose return_type='sce' or 'matrix' to preserve them exactly.",
        call. = FALSE
      )
    }

    if (data_type == "raw") {
      object <- SeuratObject::CreateSeuratObject(
        counts = expression,
        assay = assay_name,
        meta.data = cell_metadata,
        min.cells = 0,
        min.features = 0
      )
    } else {
      if (utils::packageVersion("SeuratObject") < "5.0.0") {
        stop(
          "Transformed assays require SeuratObject >= 5.0.0.",
          call. = FALSE
        )
      }

      warning(
        "This is upstream log2-scale data, not necessarily Seurat's usual ",
        "natural-log normalization. Do not apply count-based normalization or ",
        "invert it with expm1 without checking the source transformation.",
        call. = FALSE
      )

      assay <- SeuratObject::CreateAssay5Object(
        data = expression
      )

      object <- withCallingHandlers(
        SeuratObject::CreateSeuratObject(
          counts = assay,
          assay = assay_name,
          meta.data = cell_metadata
        ),
        warning = function(w) {
          if (conditionMessage(w) %in% c(
            "No layers found matching search pattern provided",
            "Cannot find the layer(s) specified"
          )) {
            invokeRestart("muffleWarning")
          }
        }
      )
    }

    object <- .abc_muffle_seurat_layer_warnings({
      object[[assay_name]] <- SeuratObject::AddMetaData(
        object[[assay_name]],
        metadata = gene_data
      )
      object
    })

    object@misc$importABCatlas <- provenance

    return(object)
  }

  assay_name_bioc <- if (data_type == "raw") {
    "counts"
  } else {
    "log2"
  }

  assays <- stats::setNames(
    list(expression),
    assay_name_bioc
  )

  row_data <- gene_data
  col_data <- cell_metadata

  if (return_type == "sce") {
    .abc_require_output_package(
      "SingleCellExperiment",
      "sce"
    )
    .abc_require_output_package(
      "S4Vectors",
      "sce"
    )

    return(
      SingleCellExperiment::SingleCellExperiment(
        assays = assays,
        rowData = S4Vectors::DataFrame(row_data),
        colData = S4Vectors::DataFrame(col_data),
        metadata = list(
          importABCatlas = provenance
        )
      )
    )
  }

  .abc_require_output_package(
    "SpatialExperiment",
    "spatial"
  )
  .abc_require_output_package(
    "S4Vectors",
    "spatial"
  )

  if (is.null(coordinate_columns) && !isTRUE(spatial)) {
    stop(
      "No declared spatial modality. Supply verified coordinate_columns ",
      "explicitly; dimensionality-reduction coordinates are not anatomical ",
      "spatial locations.",
      call. = FALSE
    )
  }

  coordinates <- .abc_spatial_coordinates(
    cell_metadata,
    coordinate_columns = coordinate_columns
  )

  provenance$coordinate_columns <- colnames(
    coordinates
  )

  SpatialExperiment::SpatialExperiment(
    assays = assays,
    rowData = S4Vectors::DataFrame(row_data),
    colData = S4Vectors::DataFrame(col_data),
    spatialCoords = coordinates,
    metadata = list(
      importABCatlas = provenance
    )
  )
}

#' Fetch selected ABC Atlas expression into an R analysis representation
#'
#' @param data Optional \code{abc_data} object returned by \code{load_data()}.
#'   This is the preferred interface. When supplied, its metadata, gene data,
#'   cache directory, and release are reused automatically.
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#'   If `NULL` (default), it is resolved by [abc_download_dir()] (option,
#'   `ABC_ATLAS_CACHE` env var, or a per-user cache directory).
#' @param metadata A data.frame containing ABC atlas cell metadata (e.g., the
#'   `cell_metadata` element returned by [load_data()]).
#' @param filters A named list of column = value pairs used to filter the metadata.
#' @param gene_data A data.frame with 'gene_identifier' and 'gene_symbol'
#'   columns, typically the `gene_data` element returned by [load_data()].
#' @param genes A character vector of gene names to include. If NULL, all genes in gene_data$gene_symbol will be used.
#' @param assay_name Assay name used for Seurat output (default: "RNA").
#' @param return_type Output format: \code{"matrix"} for a sparse expression
#'   bundle, \code{"seurat"} for a Seurat object, \code{"sce"} for a
#'   SingleCellExperiment, or \code{"spatial"} for a SpatialExperiment.
#'   Optional object formats require their corresponding packages.
#' @param coordinate_columns Optional length-two character vector naming the
#'   spatial x/y columns. Used only for \code{return_type = "spatial"}.
#' @param data_type Published expression representation (`"raw"` or `"log2"`).
#'   If `NULL`, the dataset route's declared default is used. Raw expression
#'   retains count semantics. Upstream log2 expression remains transformed
#'   data and is never presented as counts.
#' @param offline If `TRUE`, use only locally cached expression data and never
#'   contact S3 (see [abc_cache_connect()]).
#' @param release Optional ABC Atlas release version in \code{YYYYMMDD} form.
#'   Defaults to the release validated against this package. Advanced users may
#'   supply another release explicitly.
#' @return The requested expression representation: a sparse matrix bundle,
#'   Seurat object, SingleCellExperiment, or SpatialExperiment.
#' @export
fetch_data <- function(data = NULL,
                       metadata = NULL,
                       filters = list(),
                       gene_data = NULL,
                       genes = NULL,
                       assay_name = "RNA",
                       return_type = c("seurat", "sce", "matrix", "spatial"),
                       coordinate_columns = NULL,
                       data_type = NULL,
                       offline = FALSE,
                       release = NULL,
                       download_base = NULL) {
  # A historical positional first argument was download_base. Preserve that
  # form while allowing an abc_data object to be the preferred first argument.
  if (!is.null(data) &&
      !inherits(data, "abc_data") &&
      is.character(data) &&
      length(data) == 1L &&
      is.null(download_base)) {
    download_base <- data
    data <- NULL
  }

  inputs <- .abc_fetch_inputs(
    data = data,
    metadata = metadata,
    gene_data = gene_data,
    download_base = download_base,
    release = release
  )

  metadata <- inputs$metadata
  gene_data <- inputs$gene_data
  download_base <- abc_download_dir(inputs$download_base)

  if (is.null(inputs$release)) {
    release <- .abc_default_release
  } else {
    release <- inputs$release
  }

  abc_cache <- abc_cache_connect(
    download_base,
    offline = offline,
    release = release
  )

  # Apply metadata filters and resolve selected genes
  prepared <- .abc_prepare_fetch(metadata, filters, gene_data, genes)
  data_type <- .abc_resolve_data_type(data_type, prepared$metadata)
  filtered_meta <- prepared$metadata
  selected_gene_data <- prepared$gene_data

  # Fetch gene expression data using the full ordered gene universe while
  # retaining the exact resolved rows and request order.
  expression <- .abc_expression_matrix(
    abc_cache = abc_cache,
    cells = filtered_meta,
    gene_data = gene_data,
    gene_sel = selected_gene_data$.abc_gene_index,
    data_type = data_type
  )

  return_type <- match.arg(return_type)

  selected_gene_data$.abc_gene_index <- NULL

  provenance <- .abc_output_provenance(
    filtered_meta,
    selected_gene_data,
    data_type = data_type,
    release = release
  )

  .abc_materialize(
    expression = expression,
    cell_metadata = filtered_meta,
    gene_data = selected_gene_data,
    data_type = data_type,
    return_type = return_type,
    assay_name = assay_name,
    coordinate_columns = coordinate_columns,
    provenance = provenance,
    spatial = isTRUE(
      attr(filtered_meta, "abc_spatial")
    )
  )
}

#' Plan an ABC Atlas expression fetch without downloading expression files
#'
#' Applies the same cell filtering, gene selection and expression-file routing
#' as [fetch_data()], but does not download or read expression matrices.
#'
#' @param data Optional \code{abc_data} object returned by \code{load_data()}.
#'   This is the preferred interface. When supplied, its metadata, gene data,
#'   cache directory, and release are reused automatically.
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#' @param metadata Cell metadata, typically from [load_data()].
#' @param filters Named list of column = value filters.
#' @param gene_data Gene metadata, typically from [load_data()].
#' @param genes Character vector of gene symbols. If `NULL`, all genes are used.
#' @param data_type Published expression representation (`"raw"` or `"log2"`).
#'   If `NULL`, the dataset route's declared default is used.
#' @param offline If `TRUE`, use only the locally cached manifest.
#' @param release ABC Atlas release. Defaults to the package-validated release.
#' @return An `abc_fetch_plan` object describing the requested output and
#'   expression files required to produce it.
#' @export
plan_fetch <- function(data = NULL, metadata = NULL, filters = list(),
                       gene_data = NULL, genes = NULL, data_type = NULL,
                       offline = FALSE, release = NULL,
                       download_base = NULL) {
  # Backward compatibility for the historical positional download_base.
  if (!is.null(data) &&
      !inherits(data, "abc_data") &&
      is.character(data) &&
      length(data) == 1L &&
      is.null(download_base)) {
    download_base <- data
    data <- NULL
  }

  inputs <- .abc_fetch_inputs(
    data = data,
    metadata = metadata,
    gene_data = gene_data,
    download_base = download_base,
    release = release
  )

  metadata <- inputs$metadata
  gene_data <- inputs$gene_data
  download_base <- abc_download_dir(inputs$download_base)

  if (is.null(inputs$release)) {
    release <- .abc_default_release
  } else {
    release <- inputs$release
  }

  abc_cache <- abc_cache_connect(
    download_base,
    offline = offline,
    release = release
  )

  prepared <- .abc_prepare_fetch(metadata, filters, gene_data, genes)
  data_type <- .abc_resolve_data_type(data_type, prepared$metadata)
  files <- .abc_expression_files(abc_cache, prepared$metadata,
                                 data_type = data_type)

  files$additional_download <- ifelse(
    files$cached, 0,
    pmax(files$size - files$partial_bytes, 0)
  )

  structure(
    list(
      data_type = data_type,
      n_cells = nrow(prepared$metadata),
      n_genes = nrow(prepared$gene_data),
      files = files,
      total_size = sum(files$size),
      total_partial = sum(files$partial_bytes),
      total_download = sum(files$additional_download),
      release = abc_cache$version
    ),
    class = "abc_fetch_plan"
  )
}

#' @export
print.abc_fetch_plan <- function(x, ...) {
  cat(
    "<abc_fetch_plan>\n",
    "  release:        ", x$release, "\n",
    "  representation: ", x$data_type, "\n",
    "  requested:      ", x$n_genes, " genes x ", x$n_cells, " cells\n",
    sep = ""
  )

  cat("\nExpression files:\n")
  for (i in seq_len(nrow(x$files))) {
    f <- x$files[i, ]
    cat(
      "  ", f$dataset_label, "/", f$feature_matrix_label, "/", f$data_type,
      " [", f$n_cells, " cells]\n",
      "    size:                ", .abc_format_size(f$size), "\n",
      "    cached:              ", if (f$cached) "yes" else "no", "\n",
      "    resumable partial:   ", .abc_format_size(f$partial_bytes), "\n",
      "    additional download: ", .abc_format_size(f$additional_download), "\n",
      sep = ""
    )
  }

  cat(
    "\nTotals:\n",
    "  selected-file size: ", .abc_format_size(x$total_size), "\n",
    "  resumable partial:  ", .abc_format_size(x$total_partial), "\n",
    "  additional download:", " ", .abc_format_size(x$total_download), "\n",
    "\nNo expression files were downloaded.\n",
    sep = ""
  )

  invisible(x)
}

# Resolve an expression representation from dataset metadata.
.abc_resolve_data_type <- function(data_type, metadata) {
  default <- attr(metadata, "abc_default_data_type")
  available <- attr(metadata, "abc_available_data_types")

  if (is.null(default) || !length(default)) {
    default <- "raw"
  }
  if (is.null(available) || !length(available)) {
    available <- c("raw", "log2")
  }

  default <- as.character(default[[1L]])
  available <- as.character(available)

  if (is.null(data_type) || !length(data_type) ||
      !nzchar(as.character(data_type[[1L]]))) {
    data_type <- default
  }

  if (length(data_type) != 1L || is.na(data_type)) {
    stop("data_type must be one of: 'raw', 'log2'.", call. = FALSE)
  }

  data_type <- as.character(data_type)

  if (!data_type %in% c("raw", "log2")) {
    stop("data_type must be one of: 'raw', 'log2'.", call. = FALSE)
  }

  if (!data_type %in% available) {
    dataset <- attr(metadata, "abc_dataset")
    if (is.null(dataset) || !length(dataset)) {
      dataset <- "This dataset"
    }

    stop(
      dataset, " does not publish data_type='", data_type,
      "'. Available data_type values: ",
      paste(sprintf("'%s'", available), collapse = ", "),
      "; default: ", default, ".",
      call. = FALSE
    )
  }

  data_type
}

# Resolve requested genes against stable identifiers first, then source
# identifiers when available, and finally symbols. Each request must identify
# exactly one gene; absent genes are never silently replaced by zero rows.
.abc_resolve_genes <- function(gene_data, genes = NULL) {
  if (!"gene_identifier" %in% names(gene_data)) {
    stop("'gene_identifier' column not found in gene_data.", call. = FALSE)
  }

  ids <- as.character(gene_data$gene_identifier)

  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(
      "'gene_identifier' must be a non-missing unique key in gene_data.",
      call. = FALSE
    )
  }

  if (is.null(genes)) {
    out <- gene_data
    out$.abc_gene_index <- seq_len(nrow(out))
    rownames(out) <- NULL
    return(out)
  }

  query <- as.character(genes)

  if (!length(query) || anyNA(query) || any(!nzchar(query)) ||
      anyDuplicated(query)) {
    stop(
      "genes must be a nonempty vector without duplicate requests.",
      call. = FALSE
    )
  }

  source_ids <- if ("source_gene_identifier" %in% names(gene_data)) {
    as.character(gene_data$source_gene_identifier)
  } else {
    NULL
  }

  symbols <- if ("gene_symbol" %in% names(gene_data)) {
    as.character(gene_data$gene_symbol)
  } else {
    NULL
  }

  selected <- integer(length(query))

  for (i in seq_along(query)) {
    q <- query[[i]]

    matches <- which(ids == q)

    if (!length(matches) && !is.null(source_ids)) {
      matches <- which(!is.na(source_ids) & source_ids == q)
    }

    if (!length(matches) && !is.null(symbols)) {
      matches <- which(!is.na(symbols) & symbols == q)
    }

    if (!length(matches)) {
      stop(
        "Gene '", q, "' was not found.",
        call. = FALSE
      )
    }

    if (length(matches) > 1L) {
      choices <- utils::head(ids[matches], 8L)

      stop(
        "Gene '", q,
        "' is ambiguous. Please use one of these gene identifiers: ",
        paste(choices, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    selected[[i]] <- matches[[1L]]
  }

  out <- gene_data[selected, , drop = FALSE]

  if (anyDuplicated(as.character(out$gene_identifier))) {
    stop(
      "Multiple requests resolve to the same gene identifier.",
      call. = FALSE
    )
  }

  out$.abc_gene_index <- selected
  rownames(out) <- NULL
  out
}

# Select cells consistently for metadata loading, planning and fetching.
.abc_select_cells <- function(metadata, filters = list(), cell_ids = NULL,
                              n_cells = NULL, seed = 1L) {
  out <- metadata

  if (!is.list(filters) || is.null(names(filters)) ||
      any(!nzchar(names(filters)))) {
    if (length(filters)) {
      stop(
        "filters must be a named list of metadata column = allowed values.",
        call. = FALSE
      )
    }
  }

  for (filter_col in names(filters)) {
    filter_val <- filters[[filter_col]]

    if (!filter_col %in% names(out)) {
      stop(
        "Column '", filter_col,
        "' not found in metadata.",
        call. = FALSE
      )
    }

    out <- out[out[[filter_col]] %in% filter_val, , drop = FALSE]
  }

  if (!is.null(cell_ids)) {
    if (!is.character(cell_ids)) {
      stop("cell_ids must be a character vector.", call. = FALSE)
    }

    if (!length(cell_ids) || anyNA(cell_ids) ||
        any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
      stop(
        "cell_ids must be non-missing, nonempty and contain no duplicates.",
        call. = FALSE
      )
    }

    if (!"cell_label" %in% names(out)) {
      stop(
        "Cell metadata is missing required column 'cell_label'.",
        call. = FALSE
      )
    }

    labels <- as.character(out$cell_label)

    if (anyNA(labels) || anyDuplicated(labels)) {
      stop(
        "cell_label must be a non-missing unique key before cell_ids selection.",
        call. = FALSE
      )
    }

    pos <- match(cell_ids, labels)

    if (anyNA(pos)) {
      missing <- cell_ids[is.na(pos)]

      stop(
        length(missing),
        " requested cell",
        if (length(missing) == 1L) "" else "s",
        " absent after filtering: ",
        paste(utils::head(missing, 3L), collapse = ", "),
        if (length(missing) > 3L) ", ..." else "",
        call. = FALSE
      )
    }

    out <- out[pos, , drop = FALSE]
  }

  if (!is.null(n_cells)) {
    if (length(n_cells) != 1L || is.na(n_cells) ||
        n_cells < 0 || n_cells != as.integer(n_cells)) {
      stop(
        "n_cells must be a single non-negative integer.",
        call. = FALSE
      )
    }

    n_cells <- as.integer(n_cells)

    if (n_cells > nrow(out)) {
      stop(
        "n_cells=", n_cells,
        " exceeds the ", nrow(out),
        " cells available after filtering.",
        call. = FALSE
      )
    }

    if (n_cells < nrow(out)) {
      if (length(seed) != 1L || is.na(seed)) {
        stop("seed must be a single non-missing value.", call. = FALSE)
      }

      old_seed_exists <- exists(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )

      if (old_seed_exists) {
        old_seed <- get(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      }

      on.exit(
        {
          if (old_seed_exists) {
            assign(
              ".Random.seed",
              old_seed,
              envir = .GlobalEnv
            )
          } else if (exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
          )) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        },
        add = TRUE
      )

      set.seed(seed)

      keep <- sample.int(
        nrow(out),
        size = n_cells,
        replace = FALSE
      )

      out <- out[keep, , drop = FALSE]
    }
  }

  if (!nrow(out)) {
    stop(
      "No cells left after selection. Check filters, cell_ids or n_cells.",
      call. = FALSE
    )
  }

  rownames(out) <- as.character(out$cell_label)
  out
}

# Apply cell filters and resolve the requested genes for a fetch.
.abc_prepare_fetch <- function(metadata, filters, gene_data, genes = NULL) {
  filtered_meta <- .abc_select_cells(
    metadata,
    filters = filters
  )

  resolved_genes <- .abc_resolve_genes(
    gene_data,
    genes
  )

  list(
    metadata = filtered_meta,
    gene_data = resolved_genes
  )
}

# Resolve the expression files needed for `cells` without downloading them.
.abc_expression_files <- function(abc_cache, cells, data_type = "raw") {
  if (!"feature_matrix_label" %in% names(cells)) {
    stop(
      "Cell metadata is missing required column: feature_matrix_label.",
      call. = FALSE
    )
  }

  matrix_labels <- as.character(cells$feature_matrix_label)
  if (anyNA(matrix_labels) || any(!nzchar(matrix_labels))) {
    stop(
      "Selected cells contain a missing feature_matrix_label.",
      call. = FALSE
    )
  }

  # Some logical routes reuse cell metadata from one published dataset while
  # expression is published under a differently named feature matrix. Preserve
  # the published metadata column and translate labels only for internal
  # expression routing.
  matrix_map <- attr(cells, "abc_expression_matrix_map")

  if (!is.null(matrix_map) && length(matrix_map)) {
    map_names <- names(matrix_map)
    matrix_map <- stats::setNames(
      as.character(matrix_map),
      map_names
    )

    if (
      is.null(map_names) ||
      anyNA(map_names) ||
      any(!nzchar(map_names)) ||
      anyDuplicated(map_names) ||
      anyNA(matrix_map) ||
      any(!nzchar(matrix_map))
    ) {
      stop(
        "Invalid ABC expression-matrix mapping provenance.",
        call. = FALSE
      )
    }

    mapped <- unname(matrix_map[matrix_labels])
    use_map <- !is.na(mapped)
    matrix_labels[use_map] <- mapped[use_map]
  }

  expression_dirs <- attr(cells, "abc_expression_directories")

  if (is.null(expression_dirs) || !length(expression_dirs)) {
    # Backward compatibility for metadata not produced by load_data().
    # Historically dataset_label directly named the expression directory.
    if (!"dataset_label" %in% names(cells)) {
      stop(
        "Cell metadata has no ABC expression-route provenance and is missing ",
        "dataset_label. Use metadata returned by load_data(), or provide a ",
        "dataset_label column for legacy routing.",
        call. = FALSE
      )
    }

    dataset_labels <- as.character(cells$dataset_label)
    if (anyNA(dataset_labels) || any(!nzchar(dataset_labels))) {
      stop("Selected cells contain a missing dataset_label.", call. = FALSE)
    }

    pairs <- unique(data.frame(
      directory = dataset_labels,
      feature_matrix_label = matrix_labels,
      stringsAsFactors = FALSE
    ))
  } else {
    expression_dirs <- unique(as.character(expression_dirs))

    inventory <- do.call(
      rbind,
      lapply(
        expression_dirs,
        .abc_expression_inventory,
        cache = abc_cache
      )
    )

    inventory <- inventory[
      inventory$data_type == data_type,
      ,
      drop = FALSE
    ]

    labels <- unique(matrix_labels)
    pairs <- vector("list", length(labels))

    for (i in seq_along(labels)) {
      label <- labels[[i]]
      matches <- inventory[
        inventory$feature_matrix_label == label,
        c("directory", "feature_matrix_label"),
        drop = FALSE
      ]

      if (nrow(matches) == 0L) {
        stop(
          "Cannot route feature_matrix_label='", label,
          "' for data_type='", data_type,
          "' within this dataset's expression directories: ",
          paste(expression_dirs, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      if (nrow(matches) > 1L) {
        choices <- paste(
          paste0(matches$directory, "/", matches$feature_matrix_label),
          collapse = ", "
        )
        stop(
          "Ambiguous expression routing for feature_matrix_label='", label,
          "' and data_type='", data_type, "': ",
          choices, ".",
          call. = FALSE
        )
      }

      pairs[[i]] <- matches
    }

    pairs <- do.call(rbind, pairs)
  }

  pairs <- pairs[
    order(pairs$directory, pairs$feature_matrix_label, method = "radix"),
    ,
    drop = FALSE
  ]

  out <- vector("list", nrow(pairs))

  for (k in seq_len(nrow(pairs))) {
    directory <- as.character(pairs$directory[k])
    matrix_file <- as.character(pairs$feature_matrix_label[k])

    attrs <- .abc_file_attributes(
      abc_cache,
      directory,
      paste0(matrix_file, "/", data_type)
    )
    part <- paste0(attrs$local_path, ".part")

    out[[k]] <- data.frame(
      dataset_label = directory,
      feature_matrix_label = matrix_file,
      data_type = data_type,
      n_cells = sum(matrix_labels == matrix_file),
      size = as.numeric(attrs$size),
      relative_path = attrs$relative_path,
      local_path = attrs$local_path,
      cached = file.exists(attrs$local_path),
      partial_bytes = if (file.exists(part)) {
        as.numeric(file.info(part)$size)
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

# Sparse genes x cells expression matrix (dgCMatrix) for `cells` (rownames are
# cell labels) and the genes of `gene_data` whose symbol is in `selected_genes`.
#
# Cells are grouped by the expression file that holds them (dataset_label +
# feature_matrix_label); each file is downloaded if needed and read as an
# on-disk sparse matrix, from which only the requested cells and genes are
# realized, `block_size` cells at a time. Nothing is ever densified. Rows follow
# the resolved gene request order and are named by (made unique) gene symbol;
# columns
# follow the order of `cells`. Cells not found in any file are all-zero columns.
.abc_expression_matrix <- function(abc_cache, cells, gene_data, gene_sel,
                                   data_type = "raw", block_size = 50000L) {
  # Initialize HDF5Array before the first H5SparseMatrix() construction.
  # Lazy namespace initialization during H5SparseMatrix() can trigger an
  # R protection-stack imbalance warning with this dependency stack.
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop(
      "Package 'HDF5Array' is required to read expression matrices.",
      call. = FALSE
    )
  }
  gene_sel <- as.integer(gene_sel)
  gene_filtered <- gene_data[gene_sel, , drop = FALSE]
  cell_ids <- rownames(cells)

  files <- .abc_expression_files(abc_cache, cells, data_type = data_type)

  rows <- cols <- vals <- list()
  for (k in seq_len(nrow(files))) {
    matrix_file <- files$feature_matrix_label[k]
    message("loading file: ", matrix_file)

    path <- abc_cache$get_file_path(directory = files$dataset_label[k],
                                    file_name = paste0(matrix_file, "/", data_type))
    gene_pos <- .abc_gene_positions(
      gene_data,
      gene_sel,
      .h5ad_index(path, "var"),
      basename(path)
    )
    file_pos <- match(cell_ids, .h5ad_index(path, "obs"))
    wanted <- which(!is.na(file_pos))
    wanted <- wanted[order(file_pos[wanted])]  # read cells in on-disk order

    X <- HDF5Array::H5SparseMatrix(path, "/X")  # genes x cells, on disk
    for (block in split(wanted, ceiling(seq_along(wanted) / block_size))) {
      m <- methods::as(X[gene_pos, file_pos[block], drop = FALSE], "dgCMatrix")
      nz <- Matrix::summary(m)
      rows[[length(rows) + 1L]] <- nz$i
      cols[[length(cols) + 1L]] <- block[nz$j]
      vals[[length(vals) + 1L]] <- nz$x
    }
  }

  Matrix::sparseMatrix(
    i = unlist(rows, use.names = FALSE) %||% integer(),
    j = unlist(cols, use.names = FALSE) %||% integer(),
    x = unlist(vals, use.names = FALSE) %||% numeric(),
    dims = c(length(gene_sel), length(cell_ids)),
    dimnames = list(
      as.character(gene_filtered$gene_identifier),
      cell_ids
    )
  )
}

# Positions in an expression file's var index for selected gene metadata rows.
# Canonical gene_identifier is preferred. For normalized duplicate source IDs,
# occurrence metadata permits one-to-one alignment only when the H5AD contains
# exactly the expected number of source-ID occurrences. No aggregation or
# arbitrary duplicate selection is allowed.
#
# Some ABC Atlas files store malformed identifiers for a few genes: in
# Zeng-Aging-Mouse-10Xv3, genes with a duplicated symbol have a printed pandas
# Series embedding the real identifier. The existing positional fallback for
# that published-data defect is retained after ordinary identifier matching.
.abc_gene_positions <- function(gene_data, sel, var, file) {
  sel <- as.integer(sel)

  ids <- as.character(
    gene_data$gene_identifier %||% rownames(gene_data)
  )

  pos <- match(ids[sel], var)
  unresolved <- which(is.na(pos))

  if (!length(unresolved)) {
    return(pos)
  }

  has_source <- all(
    c(
      "source_gene_identifier",
      "source_gene_occurrence",
      "source_gene_group_size"
    ) %in% names(gene_data)
  )

  if (has_source) {
    for (j in unresolved) {
      row <- sel[[j]]

      source <- as.character(
        gene_data$source_gene_identifier[[row]]
      )
      occurrence <- as.integer(
        gene_data$source_gene_occurrence[[row]]
      )
      group_size <- as.integer(
        gene_data$source_gene_group_size[[row]]
      )

      if (
        is.na(source) ||
        !nzchar(source) ||
        is.na(occurrence) ||
        occurrence < 1L ||
        is.na(group_size) ||
        group_size < 1L
      ) {
        next
      }

      hits <- which(var == source)

      if (length(hits) == 1L && group_size == 1L) {
        pos[[j]] <- hits[[1L]]
        next
      }

      if (
        group_size > 1L &&
        length(hits) == group_size &&
        occurrence <= length(hits)
      ) {
        pos[[j]] <- hits[[occurrence]]
      }
    }
  }

  if (!anyNA(pos)) {
    return(pos)
  }

  # Preserve the existing AgingMouse compatibility fallback.
  if (length(var) == length(ids)) {
    bad <- which(var != ids)

    if (
      length(bad) &&
      all(mapply(
        grepl,
        ids[bad],
        var[bad],
        MoreArgs = list(fixed = TRUE)
      ))
    ) {
      warning(
        length(bad),
        " gene identifiers in ",
        file,
        " are malformed; genes were aligned by position.",
        call. = FALSE
      )
      return(sel)
    }
  }

  stop(
    sum(is.na(pos)),
    " of the selected genes are missing from ",
    file,
    "; is gene_data from the same dataset?",
    call. = FALSE
  )
}

.h5ad_index <- function(path, group) {
  index_col <- rhdf5::h5readAttributes(path, paste0("/", group))[["_index"]]
  if (is.null(index_col)) {
    stop("Cannot find the ", group, " index in ", basename(path),
         " (unsupported h5ad layout).")
  }
  as.character(rhdf5::h5read(path, paste0("/", group, "/", index_col)))
}
