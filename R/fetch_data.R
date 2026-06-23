#' Fetch data from ABC Atlas based on filters and return a Seurat object
#'
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#'   If `NULL` (default), it is resolved by [abc_download_dir()] (option,
#'   `ABC_ATLAS_CACHE` env var, or a per-user cache directory).
#' @param metadata A data.frame containing ABC atlas cell metadata (e.g., from `get_cell_metadata()`).
#' @param filters A named list of column = value pairs used to filter the metadata.
#' @param gene_data A data.frame with a 'gene_symbol' column, typically from `load_data()`.
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
  # Ensure required packages are available
  requireNamespace("reticulate")
  requireNamespace("Seurat")
  requireNamespace("Matrix")

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

  # Name cell IDs so reticulate exposes them as the pandas index downstream
  rownames(filtered_meta) <- filtered_meta$cell_label

  # Resolve selected genes
  if (is.null(genes)) {
    if (!"gene_symbol" %in% colnames(gene_data)) {
      stop("'gene_symbol' column not found in gene_data.")
    }
    genes <- gene_data$gene_symbol
  }

  # Fetch gene expression data as a sparse genes x cells matrix. The helper
  # mirrors abc_atlas_access::get_gene_data but never densifies: it slices each
  # h5ad chunk while it is still a scipy sparse matrix and accumulates COO
  # triplets, so the full cells x genes matrix is never materialized as dense.
  # reticulate auto-converts the returned scipy CSC matrix to a dgCMatrix.
  helper <- reticulate::py_run_string(.sparse_gene_data_py, convert = TRUE)
  res <- helper$get_gene_data_sparse(
    abc_atlas_cache = abc_cache,
    all_cells = filtered_meta,
    all_genes = gene_data,
    selected_genes = genes,
    data_type = "raw"
  )

  counts <- res$matrix                              # dgCMatrix, genes x cells
  rownames(counts) <- make.unique(as.character(res$genes))
  colnames(counts) <- as.character(res$cells)

  # Create Seurat object (counts stays sparse end to end)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = assay_name,
    meta.data = filtered_meta
    )

  return(seurat_obj)
}

# Python helper used by fetch_data(). Kept as an R string so it travels with the
# package regardless of how it is loaded (devtools::load_all or installed) and
# is independent of the package name.
.sparse_gene_data_py <- "
import numpy as np
import scipy.sparse as sp
import anndata


def get_gene_data_sparse(abc_atlas_cache, all_cells, all_genes, selected_genes,
                         data_type='raw', chunk_size=8192):
    '''Sparse counterpart of abc_atlas_access.get_gene_data.

    Returns a dict with a scipy CSC matrix of shape (n_genes, n_cells) plus the
    matching gene_symbol and cell_label labels. Nothing is ever densified into a
    full cells x genes array.
    '''
    selected_genes = np.asarray(selected_genes, dtype=object)
    gene_mask = np.isin(all_genes.gene_symbol.to_numpy(), selected_genes)
    gene_filtered = all_genes[gene_mask]

    cell_index = all_cells.index
    n_cells = len(cell_index)
    n_genes = int(gene_mask.sum())

    # Group cells by the expression file that holds them, as upstream does.
    matrices = all_cells.groupby(
        ['dataset_label', 'feature_matrix_label']
    )[[all_cells.columns[0]]].count()
    matrices.columns = ['cell_count']

    data_parts = []
    gene_rows = []   # row position = gene index within gene_filtered
    cell_cols = []   # col position = cell index within all_cells

    for matrix_index in matrices.index:
        directory = matrix_index[0]
        matrix_file = matrix_index[1]
        print('loading file:', matrix_file)

        file_path = abc_atlas_cache.get_file_path(
            directory=directory,
            file_name=f'{matrix_file}/{data_type}'
        )

        expression_data = anndata.read_h5ad(file_path, backed='r')
        obs = expression_data.obs
        for chunk, min_idx, max_idx in expression_data.chunked_X(chunk_size=chunk_size):
            cell_indexes = obs.index[min_idx:max_idx]
            cell_mask = cell_indexes.isin(cell_index)
            if not cell_mask.any():
                continue
            subcell_indexes = cell_indexes[cell_mask]

            # Keep the chunk sparse: slice rows (cells) then cols (genes).
            if not sp.issparse(chunk):
                chunk = sp.csr_matrix(chunk)
            else:
                chunk = chunk.tocsr()
            sub = sp.coo_matrix(chunk[np.asarray(cell_mask), :][:, gene_mask])

            # Map this chunk's local cell rows to global cell columns.
            global_cols = cell_index.get_indexer(subcell_indexes)

            data_parts.append(sub.data)
            gene_rows.append(sub.col)               # gene index, already aligned
            cell_cols.append(global_cols[sub.row])  # global cell position

        expression_data.file.close()
        del expression_data

    if data_parts:
        data = np.concatenate(data_parts)
        rows = np.concatenate(gene_rows)
        cols = np.concatenate(cell_cols)
    else:
        data = np.array([], dtype=np.float32)
        rows = np.array([], dtype=np.int64)
        cols = np.array([], dtype=np.int64)

    out = sp.coo_matrix((data, (rows, cols)), shape=(n_genes, n_cells)).tocsc()

    return {
        'matrix': out,
        'genes': gene_filtered.gene_symbol.to_numpy().astype(str),
        'cells': cell_index.to_numpy().astype(str),
    }
"
