# =============================================================================
# Unified metadata loader.
#
# Datasets differ only in declarative details (which directory holds the cell
# metadata, which auxiliary tables to join and on what keys, whether there is a
# taxonomy/region step, and the columns of interest). Those details live in the
# .abc_datasets registry below; .build_dataset() turns a spec into the assembled
# result, and load_data() is the single user-facing entry point.
# =============================================================================

# ---- Dataset registry -------------------------------------------------------
# Each entry is a declarative spec. Fields:
#   cell             list(dir, file = "cell_metadata", dtype = <named list>)
#   taxonomy         taxonomy directory (triggers the cluster reshape + join)
#   roi              list(dir, file) for the region-of-interest colour join
#   joins            ordered list of list(dir, file, key, suffix) left-joins
#   gene             list(dir, file = "gene")
#   cols_of_interest columns to summarise into `unique_values`
.abc_datasets <- list(
  WHB = list(
    cell = list(dir = "WHB-10Xv3", dtype = list(cell_label = "str")),
    taxonomy = "WHB-taxonomy",
    roi = list(dir = "WHB-10Xv3", file = "region_of_interest_structure_map"),
    gene = list(dir = "WHB-10Xv3"),
    cols_of_interest = c("feature_matrix_label", "brain_section_label",
      "region_of_interest_label", "anatomical_division_label", "donor_sex",
      "subcluster", "cluster", "supercluster", "neurotransmitter")
  ),
  WMB = list(
    cell = list(dir = "WMB-10X", dtype = list(cell_label = "str")),
    taxonomy = "WMB-taxonomy",
    gene = list(dir = "WMB-10X"),
    cols_of_interest = c("feature_matrix_label", "library_method",
      "region_of_interest_acronym", "donor_sex", "class", "subclass",
      "supertype", "cluster", "neurotransmitter")
  ),
  AgingMouse = list(
    cell = list(dir = "Zeng-Aging-Mouse-10Xv3",
                dtype = list(cell_label = "str", wmb_cluster_alias = "Int64")),
    joins = list(
      list(dir = "Zeng-Aging-Mouse-WMB-taxonomy",
           file = "cell_cluster_mapping_annotations",
           key = "cell_label", suffix = c("", "_cl_map")),
      list(dir = "Zeng-Aging-Mouse-10Xv3", file = "cell_annotation_colors",
           key = "cell_label", suffix = c("", "_cl_colors")),
      list(dir = "Zeng-Aging-Mouse-10Xv3", file = "cluster",
           key = "cluster_alias", suffix = c("", "_cl_info"))
    ),
    # Gene annotation is shared with the human WHB reference, as in the original.
    gene = list(dir = "WHB-10Xv3"),
    cols_of_interest = c("anatomical_division_label", "donor_sex", "donor_age",
      "cluster_name_cl_info", "neurotransmitter_combined_label")
  ),
  HMBA = list(
    cell = list(dir = "HMBA-10xMultiome-BG-Aligned", dtype = list(cell_label = "str")),
    joins = list(
      list(dir = "HMBA-10xMultiome-BG-Aligned", file = "donor", key = "donor_label"),
      list(dir = "HMBA-10xMultiome-BG-Aligned", file = "library",
           key = "library_label", suffix = c("", "_library_table"))
    ),
    gene = list(dir = "HMBA-10xMultiome-BG-Aligned"),
    cols_of_interest = c("species_scientific_name", "species_common_name",
      "donor_sex", "donor_age", "region_of_interest_name",
      "anatomical_division_label")
  ),
  PMDBS = list(
    cell = list(dir = "ASAP-PMDBS-10X", dtype = list(cell_label = "str")),
    joins = list(
      list(dir = "ASAP-PMDBS-10X", file = "sample", key = "sample_label"),
      list(dir = "ASAP-PMDBS-10X", file = "donor", key = "donor_label")
    ),
    # NOTE: gene annotation points at WHB-10Xv3, mirroring the original code.
    gene = list(dir = "WHB-10Xv3"),
    cols_of_interest = c("source_dataset_label", "region_of_interest_label",
      "donor_race", "donor_sex", "primary_diagnosis", "age_at_death",
      "apoe4_status", "braak_stage", "cerad_score", "cognitive_status")
  )
)

#' Load assembled cell metadata for an ABC Atlas dataset
#'
#' Single entry point for every supported dataset. The dataset is selected by
#' name; what differs between datasets (directories, joins, taxonomy, columns of
#' interest) is described declaratively in an internal registry rather than in
#' per-dataset functions.
#'
#' @param dataset Dataset name. One of \code{"WHB"}, \code{"WMB"},
#'   \code{"AgingMouse"}, \code{"HMBA"}, \code{"PMDBS"} (partial matching
#'   allowed).
#' @param download_base Directory where ABC Atlas data is downloaded and cached.
#'   If \code{NULL} (default), it is resolved by \code{\link{abc_download_dir}}.
#' @param offline If \code{TRUE}, use only locally cached data and never contact
#'   S3 (see \code{\link{abc_cache_connect}}).
#' @return A list with \code{cell_metadata}, \code{gene_data}, and
#'   \code{unique_values} (the unique entries of each column of interest, used to
#'   build filters for \code{\link{fetch_data}}).
#' @examples
#' \dontrun{
#' WHB <- load_data("WHB")
#' names(WHB$unique_values)
#' }
#' @export
load_data <- function(dataset, download_base = NULL, offline = FALSE) {
  datasets <- names(.abc_datasets)
  if (missing(dataset)) {
    stop("Specify a dataset, one of: ", paste(datasets, collapse = ", "))
  }
  dataset <- match.arg(dataset, datasets)
  spec <- .abc_datasets[[dataset]]

  requireNamespace("reticulate")
  requireNamespace("dplyr")

  download_base <- abc_download_dir(download_base)
  # Python deps are declared at load time (.onLoad) and provisioned lazily on
  # first reticulate use inside abc_cache_connect(); no explicit setup needed.
  abc_cache <- abc_cache_connect(download_base, offline = offline)

  .build_dataset(spec, abc_cache, dataset)
}

# ---- Generic builder --------------------------------------------------------

# Turn a dataset spec into the assembled list(cell_metadata, gene_data,
# unique_values). All dataset-specific behaviour comes from `spec`.
.build_dataset <- function(spec, abc_cache, dataset = "") {
  get_df <- function(dir, file, dtype = NULL) {
    args <- list(directory = dir, file_name = file)
    if (!is.null(dtype)) args$dtype <- do.call(reticulate::dict, dtype)
    do.call(abc_cache$get_metadata_dataframe, args)
  }

  cell <- get_df(spec$cell$dir, spec$cell$file %||% "cell_metadata", spec$cell$dtype)
  message(dataset, ": ", nrow(cell), " cells")
  cell_extended <- cell

  # Taxonomy reshape + join on cluster_alias (WHB / WMB).
  if (!is.null(spec$taxonomy)) {
    tax <- .build_taxonomy(abc_cache, spec$taxonomy)
    cell_extended$cluster_alias <- as.character(cell_extended$cluster_alias)
    cell_extended <- dplyr::left_join(cell_extended, tax$cluster_details,
                                      by = "cluster_alias")
    cell_extended <- dplyr::left_join(cell_extended, tax$cluster_colors,
                                      by = "cluster_alias", suffix = c("", "_color"))
  }

  # Region-of-interest colour join (WHB).
  if (!is.null(spec$roi)) {
    roi <- get_df(spec$roi$dir, spec$roi$file)
    roi$region_of_interest_label <- make.unique(as.character(roi$region_of_interest_label))
    names(roi)[names(roi) == "color_hex_triplet"] <- "region_of_interest_color"
    roi <- roi[, c("region_of_interest_label", "region_of_interest_color")]
    cell_extended <- dplyr::left_join(cell_extended, roi,
                                      by = "region_of_interest_label")
  }

  # Simple left-joins, applied in order onto the accumulating table.
  for (j in spec$joins %||% list()) {
    aux <- get_df(j$dir, j$file)
    cell_extended <- dplyr::left_join(cell_extended, aux, by = j$key,
                                      suffix = j$suffix %||% c("", ""))
  }

  gene <- get_df(spec$gene$dir, spec$gene$file %||% "gene")
  rownames(gene) <- gene$gene_identifier
  message(dataset, ": ", nrow(gene), " genes")

  cols <- intersect(spec$cols_of_interest, colnames(cell_extended))
  unique_values <- lapply(cell_extended[cols], unique)

  list(cell_metadata = cell_extended, gene_data = gene, unique_values = unique_values)
}

# Build the wide cluster_details / cluster_colors lookup tables from a taxonomy
# directory's membership + term-set tables, keyed by cluster_alias. Shared by
# every dataset that declares a `taxonomy` (currently WHB and WMB).
.build_taxonomy <- function(abc_cache, taxonomy_dir) {
  membership <- abc_cache$get_metadata_dataframe(
    directory = taxonomy_dir, file_name = "cluster_to_cluster_annotation_membership")
  term_sets <- abc_cache$get_metadata_dataframe(
    directory = taxonomy_dir, file_name = "cluster_annotation_term_set")

  widen <- function(value_col) {
    agg <- aggregate(
      stats::as.formula(paste(value_col,
                              "~ cluster_alias + cluster_annotation_term_set_name")),
      data = membership, FUN = function(x) x[1])
    w <- stats::reshape(agg, idvar = "cluster_alias",
                        timevar = "cluster_annotation_term_set_name",
                        direction = "wide")
    colnames(w) <- gsub(paste0(value_col, "."), "", colnames(w), fixed = TRUE)
    # Keep the real cluster_alias column plus the term-set columns, in order.
    w <- w[, c("cluster_alias", term_sets$name), drop = FALSE]
    w$cluster_alias <- as.character(w$cluster_alias)
    w
  }

  cluster_details <- widen("cluster_annotation_term_name")
  cluster_details[is.na(cluster_details)] <- "Other"
  cluster_colors <- widen("color_hex_triplet")

  list(cluster_details = cluster_details, cluster_colors = cluster_colors)
}
