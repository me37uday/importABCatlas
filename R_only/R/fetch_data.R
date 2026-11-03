#' Fetch data from ABC Atlas based on filters and return a Seurat object
#'
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#'   If `NULL` (default), it is resolved by [abc_download_dir()] (option,
#'   `ABC_ATLAS_CACHE` env var, or a per-user cache directory).
#' @param metadata A data.frame containing ABC atlas cell metadata (e.g., the
#'   `cell_metadata` element returned by [load_data()]).
#' @param filters A named list of column = value pairs used to filter the metadata.
#' @param gene_data A data.frame with 'gene_identifier' and 'gene_symbol'
#'   columns, typically the `gene_data` element returned by [load_data()].
#' @param genes A character vector of gene names to include. If NULL, all genes in gene_data$gene_symbol will be used.
#' @param assay_name The name to use for the Seurat assay (default: "RNA").
#' @param offline If `TRUE`, use only locally cached expression data and never
#'   contact S3 (see [abc_cache_connect()]).
#' @param release Optional ABC Atlas release version in \code{YYYYMMDD} form.
#'   Defaults to the release validated against this package. Advanced users may
#'   supply another release explicitly.
#' @return A Seurat object containing the filtered cells and gene expression data.
#' @export
fetch_data <- function(download_base = NULL,
                       metadata,
                       filters = list(),
                       gene_data,
                       genes = NULL,
                       assay_name = "RNA",
                       offline = FALSE,
                       release = .abc_default_release) {
  # Resolve the download directory and connect to the ABC Atlas cache
  download_base <- abc_download_dir(download_base)
  abc_cache <- abc_cache_connect(
   download_base,
   offline = offline,
   release = release
  )

  # Apply metadata filters and resolve selected genes
  prepared <- .abc_prepare_fetch(metadata, filters, gene_data, genes)
  filtered_meta <- prepared$metadata
  genes <- prepared$genes

  # Fetch gene expression data as a sparse genes x cells matrix
  counts <- .abc_expression_matrix(
    abc_cache = abc_cache,
    cells = filtered_meta,
    gene_data = gene_data,
    selected_genes = genes,
    data_type = "raw"
  )

  # Create Seurat object (counts stays sparse end to end)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay_name,
    meta.data = filtered_meta
    )

  return(seurat_obj)
}

#' Plan an ABC Atlas expression fetch without downloading expression files
#'
#' Applies the same cell filtering, gene selection and expression-file routing
#' as [fetch_data()], but does not download or read expression matrices.
#'
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#' @param metadata Cell metadata, typically from [load_data()].
#' @param filters Named list of column = value filters.
#' @param gene_data Gene metadata, typically from [load_data()].
#' @param genes Character vector of gene symbols. If `NULL`, all genes are used.
#' @param data_type Expression representation to plan. Currently `"raw"`.
#' @param offline If `TRUE`, use only the locally cached manifest.
#' @param release ABC Atlas release. Defaults to the package-validated release.
#' @return An `abc_fetch_plan` object describing the requested output and
#'   expression files required to produce it.
#' @export
plan_fetch <- function(download_base = NULL, metadata, filters = list(),
                       gene_data, genes = NULL, data_type = "raw",
                       offline = FALSE, release = .abc_default_release) {
  download_base <- abc_download_dir(download_base)
  abc_cache <- abc_cache_connect(download_base, offline = offline,
                                 release = release)

  prepared <- .abc_prepare_fetch(metadata, filters, gene_data, genes)
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
      n_genes = sum(gene_data$gene_symbol %in% prepared$genes),
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

# Apply cell filters and resolve the requested genes for a fetch.
.abc_prepare_fetch <- function(metadata, filters, gene_data, genes = NULL) {
  filtered_meta <- metadata
  for (filter_col in names(filters)) {
    filter_val <- filters[[filter_col]]
    if (!filter_col %in% colnames(filtered_meta)) {
      stop(paste("Column", filter_col, "not found in metadata."))
    }
    filtered_meta <- filtered_meta[filtered_meta[[filter_col]] %in% filter_val, ]
  }

  if (nrow(filtered_meta) == 0) {
    stop("No cells left after filtering. Check your filter values.")
  }

  rownames(filtered_meta) <- filtered_meta$cell_label

  if (is.null(genes)) {
    if (!"gene_symbol" %in% colnames(gene_data)) {
      stop("'gene_symbol' column not found in gene_data.")
    }
    genes <- gene_data$gene_symbol
  }

  list(metadata = filtered_meta, genes = genes)
}

# Resolve the expression files needed for `cells` without downloading them.
.abc_expression_files <- function(abc_cache, cells, data_type = "raw") {
  required <- c("dataset_label", "feature_matrix_label")
  missing <- setdiff(required, names(cells))
  if (length(missing)) {
    stop("Cell metadata is missing required column(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  files <- unique(cells[, required, drop = FALSE])
  files <- files[order(files$dataset_label, files$feature_matrix_label,
                       method = "radix"), , drop = FALSE]

  out <- vector("list", nrow(files))
  for (k in seq_len(nrow(files))) {
    directory <- as.character(files$dataset_label[k])
    matrix_file <- as.character(files$feature_matrix_label[k])
    attrs <- .abc_file_attributes(abc_cache, directory,
                                  paste0(matrix_file, "/", data_type))
    part <- paste0(attrs$local_path, ".part")

    out[[k]] <- data.frame(
      dataset_label = directory,
      feature_matrix_label = matrix_file,
      data_type = data_type,
      n_cells = sum(cells$dataset_label == directory &
                    cells$feature_matrix_label == matrix_file, na.rm = TRUE),
      size = as.numeric(attrs$size),
      relative_path = attrs$relative_path,
      local_path = attrs$local_path,
      cached = file.exists(attrs$local_path),
      partial_bytes = if (file.exists(part)) as.numeric(file.info(part)$size) else 0,
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
# the order of `gene_data` and are named by (made unique) gene symbol; columns
# follow the order of `cells`. Cells not found in any file are all-zero columns.
.abc_expression_matrix <- function(abc_cache, cells, gene_data, selected_genes,
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
  gene_sel <- which(gene_data$gene_symbol %in% selected_genes)
  gene_filtered <- gene_data[gene_sel, , drop = FALSE]
  all_gene_ids <- as.character(gene_data$gene_identifier %||% rownames(gene_data))
  cell_ids <- rownames(cells)

  files <- .abc_expression_files(abc_cache, cells, data_type = data_type)

  rows <- cols <- vals <- list()
  for (k in seq_len(nrow(files))) {
    matrix_file <- files$feature_matrix_label[k]
    message("loading file: ", matrix_file)

    path <- abc_cache$get_file_path(directory = files$dataset_label[k],
                                    file_name = paste0(matrix_file, "/", data_type))
    gene_pos <- .abc_gene_positions(all_gene_ids, gene_sel,
                                    .h5ad_index(path, "var"), basename(path))
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
    dimnames = list(make.unique(as.character(gene_filtered$gene_symbol)), cell_ids)
  )
}

# Positions in an expression file's var index of the genes `ids[sel]`, where
# `ids` are all gene identifiers of gene_data. Genes are matched by identifier.
# Some ABC Atlas files store malformed identifiers for a few genes: in
# Zeng-Aging-Mouse-10Xv3, genes with a duplicated symbol have a printed pandas
# Series (embedding the real identifier) instead. When every identifier in the
# file either equals or embeds the gene_data identifier at the same position,
# genes are aligned by position instead, as abc_atlas_access does.
.abc_gene_positions <- function(ids, sel, var, file) {
  pos <- match(ids[sel], var)
  if (!anyNA(pos)) return(pos)
  if (length(var) == length(ids)) {
    bad <- which(var != ids)
    if (all(mapply(grepl, ids[bad], var[bad], fixed = TRUE))) {
      warning(length(bad), " gene identifiers in ", file, " are malformed; ",
              "genes were aligned by position.", call. = FALSE)
      return(sel)
    }
  }
  stop(sum(is.na(pos)), " of the selected genes are missing from ", file,
       "; is gene_data from the same dataset?")
}

# Row names of an h5ad's obs or var dataframe (the column named by `_index`).
.h5ad_index <- function(path, group) {
  index_col <- rhdf5::h5readAttributes(path, paste0("/", group))[["_index"]]
  if (is.null(index_col)) {
    stop("Cannot find the ", group, " index in ", basename(path),
         " (unsupported h5ad layout).")
  }
  as.character(rhdf5::h5read(path, paste0("/", group, "/", index_col)))
}
