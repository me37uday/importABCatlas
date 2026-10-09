#' Load and validate ABC Atlas cell and gene metadata
#'
#' Metadata files may be large and are downloaded in full on first use. Supplying
#' exact cell_ids streams CSV filtering; n_cells samples only after annotation.
#' preview_n is exploratory, not a representative biological sample.
#' @param dataset Registry ID from list_datasets or an exact manifest directory.
#' @param filters Named list of exact allowed values in metadata columns.
#' @param cell_ids Character vector of exact cell labels, preserving requested order.
#' @param n_cells Optional fixed sample size after filtering.
#' @param seed Reproducible sampling seed.
#' @param annotations Add validated donor, taxonomy and spatial annotations.
#' @param overrides Named list overriding registry fields for a verified new schema.
#' @inheritParams get_metadata
#' @return An abc_data list with cell_metadata, gene_data, unique_values and provenance.
#' @export
load_data <- function(dataset, download_base = "~/abc_atlas_cache", release = .default_release,
                      offline = FALSE, filters = list(), cell_ids = NULL, n_cells = NULL,
                      seed = 1L, preview_n = NULL, annotations = TRUE, overrides = list()) {
  if (!is.null(cell_ids) && !is.character(cell_ids)) stop("cell_ids must be character, not numeric (large IDs lose precision).", call. = FALSE)
  args <- c(.cache_args(download_base, release, offline),
            list(dataset = dataset, filters = filters, cell_ids = cell_ids, n_cells = n_cells,
                 seed = seed, preview_n = preview_n, annotations = annotations, overrides = overrides))
  .abc_backend("load", args, function(out, result) {
    cm <- .read_frame(file.path(out, "cell_metadata.csv"))
    gd <- .read_frame(file.path(out, "gene_data.csv"))
    attr(cm, "abc_dataset") <- dataset
    attr(cm, "abc_manifest") <- result$manifest
    attr(cm, "abc_overrides") <- overrides
    # Avoid duplicating millions of cell labels in unique_values.
    categories <- names(cm)[vapply(cm, function(x) is.character(x) && length(unique(x)) <= 200L, logical(1))]
    structure(list(cell_metadata = cm, gene_data = gd,
                   unique_values = lapply(cm[categories], unique), provenance = result,
                   dataset = dataset, download_base = args$download_base, release = result$manifest,
                   overrides = overrides), class = "abc_data")
  })
}

#' @export
print.abc_data <- function(x, ...) {
  cat(sprintf("ABC Atlas metadata: %s\n%d cells; %d candidate genes\nManifest: %s\nExpression is not loaded. Use plan_fetch(data=x).\n",
              x$dataset, nrow(x$cell_metadata), nrow(x$gene_data), x$release))
  invisible(x)
}
