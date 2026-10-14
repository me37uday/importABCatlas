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
#' @return A Seurat object containing the filtered cells and gene expression data.
#' @export
fetch_data <- function(download_base = NULL,
                       metadata,
                       filters = list(),
                       gene_data,
                       genes = NULL,
                       assay_name = "RNA",
                       offline = FALSE) {
  # Resolve the download directory and connect to the ABC Atlas cache
  download_base <- abc_download_dir(download_base)
  abc_cache <- abc_cache_connect(download_base, offline = offline)

  # Apply metadata filters
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

  # Cell IDs become the Seurat cell names
  rownames(filtered_meta) <- filtered_meta$cell_label

  # Resolve selected genes
  if (is.null(genes)) {
    if (!"gene_symbol" %in% colnames(gene_data)) {
      stop("'gene_symbol' column not found in gene_data.")
    }
    genes <- gene_data$gene_symbol
  }

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
  gene_sel <- which(gene_data$gene_symbol %in% selected_genes)
  gene_filtered <- gene_data[gene_sel, , drop = FALSE]
  all_gene_ids <- as.character(gene_data$gene_identifier %||% rownames(gene_data))
  cell_ids <- rownames(cells)

  files <- unique(cells[, c("dataset_label", "feature_matrix_label")])
  files <- files[order(files$dataset_label, files$feature_matrix_label,
                       method = "radix"), , drop = FALSE]

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
