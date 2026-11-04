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
#   cell             list(dir, file = "cell_metadata", col_classes = <named chr>)
#   taxonomy         taxonomy directory (triggers the cluster reshape + join)
#   roi              list(dir, file) for the region-of-interest colour join
#   joins            ordered list of list(dir, file, key, suffix) left-joins
#   gene             list(dir, file = "gene")
#   cols_of_interest columns to summarise into `unique_values`
.abc_datasets <- list(
  WHB = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "WHB-10Xv3",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(dir = "WHB-10Xv3", col_classes = c(cell_label = "character")),
    taxonomy = "WHB-taxonomy",
    roi = list(dir = "WHB-10Xv3", file = "region_of_interest_structure_map"),
    gene = list(dir = "WHB-10Xv3"),
    cols_of_interest = c("feature_matrix_label", "brain_section_label",
      "region_of_interest_label", "anatomical_division_label", "donor_sex",
      "subcluster", "cluster", "supercluster", "neurotransmitter")
  ),
  WMB = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = c("WMB-10Xv2", "WMB-10Xv3",
                               "WMB-10XMulti"),
    modality = "single_cell_transcriptomics",
    spatial = FALSE,
    cell = list(dir = "WMB-10X", col_classes = c(cell_label = "character")),
    taxonomy = "WMB-taxonomy",
    gene = list(dir = "WMB-10X"),
    cols_of_interest = c("feature_matrix_label", "library_method",
      "region_of_interest_acronym", "donor_sex", "class", "subclass",
      "supertype", "cluster", "neurotransmitter")
  ),
  ConsensusMouseAIBS = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = c(
      "WMB-10Xv2",
      "WMB-10Xv3",
      "WMB-10XMulti"
    ),
    modality = "single_cell_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "Consensus-WMB-AIBS-10X",
      col_classes = c(cell_label = "character")
    ),
    taxonomy = "Consensus-WMB-integrated-taxonomy",
    gene = list(dir = "WMB-10X"),
    cols_of_interest = character()
  ),

  ConsensusMouseMacosko = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "Consensus-WMB-Macosko-10X",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "Consensus-WMB-Macosko-10X",
      col_classes = c(cell_label = "character")
    ),
    taxonomy = "Consensus-WMB-integrated-taxonomy",
    gene = list(dir = "Consensus-WMB-Macosko-10X"),
    cols_of_interest = character()
  ),

  DevelopingMouse = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "Developing-Mouse-Vis-Cortex-10X",
    modality = "single_cell_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "Developing-Mouse-Vis-Cortex-10X",
      col_classes = c(cell_label = "character")
    ),
    taxonomy = "Developing-Mouse-Vis-Cortex-taxonomy",
    gene = list(dir = "Developing-Mouse-Vis-Cortex-10X"),
    cols_of_interest = c("synchronized_age")
  ),

  AgingMouse = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "Zeng-Aging-Mouse-10Xv3",
    modality = "single_cell_transcriptomics",
    spatial = FALSE,
    cell = list(dir = "Zeng-Aging-Mouse-10Xv3",
                col_classes = c(cell_label = "character",
                                wmb_cluster_alias = "integer")),
    # Match the validated Python-reticulate annotation path: native Aging
    # annotations define cluster identity and age bias, while the separate
    # mapping directory contributes the optional WMB reference mapping.
    joins = list(
      list(
        dir = "Zeng-Aging-Mouse-10Xv3",
        file = "cell_cluster_annotations",
        key = "cell_label"
      ),
      list(
        dir = "Zeng-Aging-Mouse-WMB-taxonomy",
        file = "cell_cross_mapping_annotations",
        key = "cell_label"
      )
    ),
    # No gene table of its own: the expression matrix uses the WMB 10X gene
    # set (same 32285 genes, same order).
    gene = list(dir = "WMB-10X"),
    cols_of_interest = c(
      "anatomical_division_label",
      "donor_sex",
      "donor_age",
      "cluster_name",
      "cluster_age_bias",
      "neurotransmitter_combined_label",
      "wmb_class_name",
      "wmb_subclass_name",
      "wmb_supertype_name",
      "wmb_cluster_name"
    )
  ),
  HMBA = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "HMBA-10xMultiome-BG-Aligned",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(dir = "HMBA-10xMultiome-BG-Aligned",
                col_classes = c(cell_label = "character")),
    joins = list(
      list(dir = "HMBA-10xMultiome-BG-Aligned", file = "donor", key = "donor_label"),
      list(
        dir = "HMBA-10xMultiome-BG-Aligned",
        file = "library",
        key = "library_label",
        column_suffix = "_library_table"
      )
    ),
    gene = list(dir = "HMBA-10xMultiome-BG-Aligned"),
    cols_of_interest = c("species_scientific_name", "species_common_name",
      "donor_sex", "donor_age", "region_of_interest_name",
      "anatomical_division_label")
  ),
  HMBA_Human = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "HMBA-10xMultiome-BG",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "HMBA-10xMultiome-BG",
      col_classes = c(cell_label = "character")
    ),
    subset = list(
      dir = "HMBA-10xMultiome-BG",
      file = "donor",
      key = "donor_label",
      column = "species_genus",
      value = "Human"
    ),
    taxonomy = "HMBA-BG-taxonomy-CCN20250428",
    gene = list(
      dir = "HMBA-10xMultiome-BG",
      file = "human_gene"
    ),
    cols_of_interest = c(
      "species_genus",
      "species_scientific_name",
      "species_common_name",
      "donor_sex",
      "donor_age",
      "region_of_interest_name",
      "anatomical_division_label"
    )
  ),

  HMBA_Macaque = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "HMBA-10xMultiome-BG",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "HMBA-10xMultiome-BG",
      col_classes = c(cell_label = "character")
    ),
    subset = list(
      dir = "HMBA-10xMultiome-BG",
      file = "donor",
      key = "donor_label",
      column = "species_genus",
      value = "Macaque"
    ),
    taxonomy = "HMBA-BG-taxonomy-CCN20250428",
    gene = list(
      dir = "HMBA-10xMultiome-BG",
      file = "macaque_gene"
    ),
    cols_of_interest = c(
      "species_genus",
      "species_scientific_name",
      "species_common_name",
      "donor_sex",
      "donor_age",
      "region_of_interest_name",
      "anatomical_division_label"
    )
  ),

  HMBA_Marmoset = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "HMBA-10xMultiome-BG",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "HMBA-10xMultiome-BG",
      col_classes = c(cell_label = "character")
    ),
    subset = list(
      dir = "HMBA-10xMultiome-BG",
      file = "donor",
      key = "donor_label",
      column = "species_genus",
      value = "Marmoset"
    ),
    taxonomy = "HMBA-BG-taxonomy-CCN20250428",
    gene = list(
      dir = "HMBA-10xMultiome-BG",
      file = "marmoset_gene"
    ),
    cols_of_interest = c(
      "species_genus",
      "species_scientific_name",
      "species_common_name",
      "donor_sex",
      "donor_age",
      "region_of_interest_name",
      "anatomical_division_label"
    )
  ),

  MERFISH = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "MERFISH-C57BL6J-638850",
    modality = "spatial_transcriptomics",
    spatial = TRUE,
    cell = list(
      dir = "MERFISH-C57BL6J-638850",
      col_classes = c(cell_label = "character")
    ),
    taxonomy = "WMB-taxonomy",
    gene = list(dir = "MERFISH-C57BL6J-638850"),
    cols_of_interest = c(
      "brain_section_label",
      "donor_sex",
      "donor_genotype",
      "x",
      "y",
      "z",
      "class",
      "subclass",
      "supertype",
      "cluster",
      "neurotransmitter"
    )
  ),

  MERFISH_imputed = list(
    default_data_type = "log2",
    available_data_types = "log2",
    expression_directories = "MERFISH-C57BL6J-638850-imputed",
    expression_matrix_map = c(
      "C57BL6J-638850" = "C57BL6J-638850-imputed"
    ),
    modality = "imputed_spatial_transcriptomics",
    spatial = TRUE,
    cell = list(
      dir = "MERFISH-C57BL6J-638850",
      col_classes = c(cell_label = "character")
    ),
    taxonomy = "WMB-taxonomy",
    gene = list(dir = "WMB-10X"),
    cols_of_interest = c(
      "brain_section_label",
      "donor_sex",
      "donor_genotype",
      "x",
      "y",
      "z",
      "class",
      "subclass",
      "supertype",
      "cluster",
      "neurotransmitter"
    )
  ),

  SEAAD_Multiregion = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "SEA-AD-Multiregion-10X",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "SEA-AD-Multiregion-10X",
      col_classes = c(cell_label = "character")
    ),
    auto_joins = TRUE,
    taxonomy = "SEA-AD-Multiregion-taxonomy",
    gene = list(dir = "SEA-AD-Multiregion-10X"),
    cols_of_interest = c(
      "feature_matrix_label",
      "region_of_interest_label",
      "region_of_interest_name",
      "Brain Region",
      "species_scientific_name",
      "species_genus",
      "donor_sex",
      "donor_age",
      "donor_race",
      "Overall AD neuropathological Change",
      "Braak",
      "CERAD score",
      "Cognitive Status",
      "APOE Genotype"
    )
  ),

  SEAAD_CaH = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "SEA-AD-CaH-10X",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(
      dir = "SEA-AD-CaH-10X",
      col_classes = c(cell_label = "character")
    ),
    auto_joins = TRUE,
    taxonomy = "SEA-AD-CaH-taxonomy",
    gene = list(dir = "SEA-AD-CaH-10X"),
    cols_of_interest = c(
      "feature_matrix_label",
      "region_of_interest_label",
      "region_of_interest_name",
      "Brain Region",
      "species_scientific_name",
      "species_genus",
      "donor_sex",
      "donor_age",
      "donor_race",
      "Overall AD neuropathological Change",
      "Braak",
      "CERAD score",
      "Cognitive Status",
      "APOE Genotype"
    )
  ),

  PMDBS = list(
    default_data_type = "raw",
    available_data_types = c("raw", "log2"),
    expression_directories = "ASAP-PMDBS-10X",
    modality = "single_nucleus_transcriptomics",
    spatial = FALSE,
    cell = list(dir = "ASAP-PMDBS-10X", col_classes = c(cell_label = "character")),
    joins = list(
      list(dir = "ASAP-PMDBS-10X", file = "sample", key = "sample_label"),
      list(dir = "ASAP-PMDBS-10X", file = "donor", key = "donor_label")
    ),
    # PMDBS expression matrices use the native 36601-gene human gene universe.
    gene = list(dir = "ASAP-PMDBS-10X"),
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
#' @param release Optional ABC Atlas release version in \code{YYYYMMDD} form.
#'   When supplied, metadata are loaded from exactly that release.
#' @return A list with \code{cell_metadata}, \code{gene_data}, and
#'   \code{unique_values} (the unique entries of each column of interest, used to
#'   build filters for \code{\link{fetch_data}}).
#' @examples
#' \dontrun{
#' WHB <- load_data("WHB")
#' names(WHB$unique_values)
#' }
#' @export
load_data <- function(dataset, download_base = NULL, offline = FALSE,
                      release = .abc_default_release) {
  datasets <- names(.abc_datasets)
  if (missing(dataset)) {
    stop("Specify a dataset, one of: ", paste(datasets, collapse = ", "))
  }
  dataset <- match.arg(dataset, datasets)
  spec <- .abc_datasets[[dataset]]

  download_base <- abc_download_dir(download_base)

  abc_cache <- abc_cache_connect(
  download_base,
  offline = offline,
  release = release
  )
  .build_dataset(spec, abc_cache, dataset)
}

# ---- Generic builder --------------------------------------------------------

# Left many-to-one metadata join that preserves cell count and row order.
# Annotation tables must have a non-missing unique join key; otherwise joining
# them could silently duplicate cells.
.abc_safe_join <- function(left, right, key, suffix = "annotation",
                           column_suffix = paste0("__", suffix)) {
  if (!key %in% names(left) || !key %in% names(right)) {
    stop(suffix, ": join key ", key, " missing.", call. = FALSE)
  }

  left[[key]] <- as.character(left[[key]])
  right[[key]] <- as.character(right[[key]])
  right <- unique(right)

  if (anyNA(right[[key]]) || anyDuplicated(right[[key]])) {
    stop(
      suffix, ": right-hand ", key,
      " is not a non-missing unique key; join would multiply cells.",
      call. = FALSE
    )
  }

  n_before <- nrow(left)
  left$.abc_join_order <- seq_len(n_before)

  overlap <- setdiff(intersect(names(left), names(right)), key)
  if (length(overlap)) {
    replacement <- paste0(overlap, column_suffix)

    if (any(replacement %in% union(names(left), names(right)))) {
      stop(suffix, ": ambiguous annotation column suffix.", call. = FALSE)
    }

    names(right)[match(overlap, names(right))] <- replacement
  }

  out <- dplyr::left_join(left, right, by = key)
  out <- out[order(out$.abc_join_order), , drop = FALSE]
  out$.abc_join_order <- NULL

  if (nrow(out) != n_before) {
    stop("Join changed the number of cells.", call. = FALSE)
  }

  out
}

# Make published gene identifiers unique without dropping or aggregating
# features. Duplicate source IDs retain their original identity and occurrence
# so expression readers can reconcile them one-to-one against the H5AD axis.
.abc_normalize_gene_metadata <- function(gene, context = "gene metadata") {
  if (!"gene_identifier" %in% names(gene)) {
    stop(
      context, ": required column 'gene_identifier' missing.",
      call. = FALSE
    )
  }

  out <- gene
  original <- as.character(out$gene_identifier)

  if (anyNA(original) || any(!nzchar(original))) {
    stop(context, ": missing gene_identifier.", call. = FALSE)
  }

  counts <- table(original, useNA = "no")
  duplicate_groups <- counts[counts > 1L]

  if (!length(duplicate_groups)) {
    out$gene_identifier <- original
    return(out)
  }

  if ("source_gene_identifier" %in% names(out)) {
    stop(
      context,
      ": source_gene_identifier already exists; ",
      "cannot safely disambiguate duplicates.",
      call. = FALSE
    )
  }

  used <- unique(original)
  seen <- integer()
  unique_ids <- character(length(original))
  occurrence <- integer(length(original))

  for (i in seq_along(original)) {
    value <- original[[i]]

    if (counts[[value]] <= 1L) {
      unique_ids[[i]] <- value
      occurrence[[i]] <- 1L
      next
    }

    previous <- if (value %in% names(seen)) seen[[value]] else 0L
    index <- previous + 1L
    candidate <- paste0(value, ".iabdup", index)

    while (candidate %in% used) {
      index <- index + 1L
      candidate <- paste0(value, ".iabdup", index)
    }

    seen[[value]] <- index
    used <- c(used, candidate)
    unique_ids[[i]] <- candidate

    # True ordinal within the original duplicate group. This deliberately
    # differs from the suffix if an existing identifier forced a suffix skip.
    occurrence[[i]] <- sum(original[seq_len(i)] == value)
  }

  group_size <- as.integer(counts[original])
  out$source_gene_identifier <- original
  out$source_gene_occurrence <- occurrence
  out$source_gene_group_size <- group_size
  out$gene_identifier <- unique_ids

  if (anyNA(out$gene_identifier) ||
      any(!nzchar(out$gene_identifier)) ||
      anyDuplicated(out$gene_identifier)) {
    stop(
      context,
      ": duplicate gene identifiers remain after normalization; ",
      "no silent aggregation allowed.",
      call. = FALSE
    )
  }

  out
}

# Turn a dataset spec into the assembled list(cell_metadata, gene_data,
# unique_values). All dataset-specific behaviour comes from `spec`.
.build_dataset <- function(spec, abc_cache, dataset = "") {
  get_df <- function(dir, file, col_classes = NULL) {
    abc_cache$get_metadata_dataframe(directory = dir, file_name = file,
                                     col_classes = col_classes)
  }

  cell <- get_df(spec$cell$dir, spec$cell$file %||% "cell_metadata",
                 spec$cell$col_classes)
  cell_extended <- cell

  # Some logical routes share one upstream cell table but represent a strict
  # subset of it. Resolve that population before expensive taxonomy joins.
  # This is internal route configuration, not a user-facing filter.
  if (!is.null(spec$subset)) {
    subset_spec <- spec$subset

    subset_metadata <- get_df(
      subset_spec$dir,
      subset_spec$file
    )

    required <- c(subset_spec$key, subset_spec$column)
    missing <- setdiff(required, names(subset_metadata))

    if (length(missing)) {
      stop(
        dataset,
        ": subset metadata is missing required column(s): ",
        paste(missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    keep_keys <- unique(
      as.character(
        subset_metadata[
          !is.na(subset_metadata[[subset_spec$column]]) &
            as.character(subset_metadata[[subset_spec$column]]) ==
              subset_spec$value,
          subset_spec$key
        ]
      )
    )

    keep_keys <- keep_keys[!is.na(keep_keys) & nzchar(keep_keys)]

    if (!length(keep_keys)) {
      stop(
        dataset,
        ": no ",
        subset_spec$value,
        " cells found from ",
        subset_spec$file,
        ".",
        call. = FALSE
      )
    }

    cell_keys <- as.character(cell_extended[[subset_spec$key]])
    cell_extended <- cell_extended[
      !is.na(cell_keys) & cell_keys %in% keep_keys,
      ,
      drop = FALSE
    ]

    # Attach the same route-defining metadata after filtering so species and
    # donor variables remain available to users without another manual join.
    cell_extended <- .abc_safe_join(
      cell_extended,
      subset_metadata,
      subset_spec$key,
      subset_spec$file
    )
  }

  # Some datasets publish standard metadata tables separately from the cell
  # table. Join those automatically when requested by the route, matching the
  # reference backend while keeping this complexity out of the user API.
  if (isTRUE(spec$auto_joins)) {
    standard_joins <- list(
      list(file = "library", key = "library_label"),
      list(file = "sample", key = "sample_label"),
      list(file = "donor", key = "donor_label"),
      list(file = "disease", key = "donor_label"),
      list(file = "specimen_metadata", key = "brain_section_label")
    )

    for (j in standard_joins) {
      if (
        j$key %in% names(cell_extended) &&
        .abc_metadata_exists(
          abc_cache,
          spec$cell$dir,
          j$file
        )
      ) {
        aux <- get_df(
          spec$cell$dir,
          j$file
        )

        cell_extended <- .abc_safe_join(
          cell_extended,
          aux,
          j$key,
          j$file
        )
      }
    }
  }

  message(dataset, ": ", nrow(cell_extended), " cells")

  # Taxonomy annotation. Some routes already carry cluster_alias in their
  # cell table; others publish authoritative cell_to_cluster_membership in
  # either the route directory or taxonomy directory.
  if (!is.null(spec$taxonomy)) {
    membership_dir <- NULL

    for (candidate in c(spec$cell$dir, spec$taxonomy)) {
      if (.abc_metadata_exists(
        abc_cache,
        candidate,
        "cell_to_cluster_membership"
      )) {
        membership_dir <- candidate
        break
      }
    }

    if (!is.null(membership_dir)) {
      membership <- get_df(
        membership_dir,
        "cell_to_cluster_membership"
      )

      if ("cluster_alias" %in% names(cell_extended)) {
        names(cell_extended)[
          names(cell_extended) == "cluster_alias"
        ] <- "cluster_alias_original"

        cell_extended$cluster_alias_original <-
          as.character(cell_extended$cluster_alias_original)
      }

      cell_extended <- .abc_safe_join(
        cell_extended,
        membership,
        "cell_label",
        "membership"
      )
    }

    if (!"cluster_alias" %in% names(cell_extended)) {
      stop(
        spec$cell$dir,
        ": no cluster_alias or applicable membership table.",
        call. = FALSE
      )
    }

    cell_extended$cluster_alias <-
      as.character(cell_extended$cluster_alias)

    tax <- .build_taxonomy(abc_cache, spec$taxonomy)

    cell_extended <- .abc_safe_join(
      cell_extended,
      tax$cluster_details,
      "cluster_alias",
      "taxonomy"
    )

    if (!is.null(tax$cluster_colors)) {
      cell_extended <- .abc_safe_join(
        cell_extended,
        tax$cluster_colors,
        "cluster_alias",
        "taxonomy_color",
        column_suffix = "_color"
      )
    }
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

  # Metadata joins are many-to-one onto cells. Keep the existing cell-side
  # columns authoritative and give overlapping annotation columns explicit
  # provenance-bearing suffixes rather than allowing a join to multiply cells.
  for (j in spec$joins %||% list()) {
    aux <- get_df(j$dir, j$file)

    cell_extended <- .abc_safe_join(
      cell_extended,
      aux,
      j$key,
      j$file,
      column_suffix = j$column_suffix %||% paste0("__", j$file)
    )
  }

  gene <- get_df(spec$gene$dir, spec$gene$file %||% "gene")
  gene <- .abc_normalize_gene_metadata(
    gene,
    context = spec$gene$dir
  )
  rownames(gene) <- gene$gene_identifier
  message(dataset, ": ", nrow(gene), " genes")

  cols <- intersect(spec$cols_of_interest, colnames(cell_extended))
  unique_values <- lapply(cell_extended[cols], unique)

  attr(cell_extended, "abc_dataset") <- dataset
  attr(cell_extended, "abc_default_data_type") <- spec$default_data_type
  attr(cell_extended, "abc_available_data_types") <- spec$available_data_types
  attr(cell_extended, "abc_expression_directories") <- spec$expression_directories
  attr(cell_extended, "abc_expression_matrix_map") <- spec$expression_matrix_map
  attr(cell_extended, "abc_modality") <- spec$modality
  attr(cell_extended, "abc_spatial") <- spec$spatial

  list(cell_metadata = cell_extended, gene_data = gene, unique_values = unique_values)
}

# Build the wide cluster_details / cluster_colors lookup tables from a taxonomy
# directory's membership + term-set tables, keyed by cluster_alias. Shared by
# every dataset that declares a `taxonomy` (currently WHB and WMB).
.build_taxonomy <- function(abc_cache, taxonomy_dir) {
  get_df <- function(file) {
    abc_cache$get_metadata_dataframe(
      directory = taxonomy_dir,
      file_name = file
    )
  }

  has <- function(file) {
    .abc_metadata_exists(abc_cache, taxonomy_dir, file)
  }

  validate_alias <- function(x, source) {
    if (!"cluster_alias" %in% names(x)) {
      stop(
        taxonomy_dir, ": ", source, " is missing cluster_alias.",
        call. = FALSE
      )
    }

    x$cluster_alias <- as.character(x$cluster_alias)

    if (anyNA(x$cluster_alias) || anyDuplicated(x$cluster_alias)) {
      stop(
        taxonomy_dir, ": ", source,
        " does not have a non-missing unique cluster_alias.",
        call. = FALSE
      )
    }

    x
  }

  # Load the long membership when available. It is required for taxonomy
  # colours and is the fallback source for annotation names when no upstream
  # pivoted view is published.
  membership <- NULL

  if (has("cluster_to_cluster_annotation_membership")) {
    membership <- get_df("cluster_to_cluster_annotation_membership")

    # Some taxonomy releases store only term labels in the membership table.
    # Recover human-readable names from the accompanying term tables, matching
    # the validated Python backend rather than inventing annotations.
    if (!"cluster_annotation_term_name" %in% names(membership)) {
      if (!has("cluster_annotation_term")) {
        stop(
          taxonomy_dir,
          ": term names absent and no term table available.",
          call. = FALSE
        )
      }

      term <- get_df("cluster_annotation_term")
      names(term)[names(term) == "label"] <-
        "cluster_annotation_term_label"
      names(term)[names(term) == "name"] <-
        "cluster_annotation_term_name"

      membership <- .abc_safe_join(
        membership,
        term[, c(
          "cluster_annotation_term_label",
          "cluster_annotation_term_name"
        ), drop = FALSE],
        "cluster_annotation_term_label",
        "term"
      )
    }

    if (!"cluster_annotation_term_set_name" %in% names(membership)) {
      if (!has("cluster_annotation_term_set")) {
        stop(
          taxonomy_dir,
          ": term-set names absent and no term-set table available.",
          call. = FALSE
        )
      }

      term_sets <- get_df("cluster_annotation_term_set")
      names(term_sets)[names(term_sets) == "label"] <-
        "cluster_annotation_term_set_label"
      names(term_sets)[names(term_sets) == "name"] <-
        "cluster_annotation_term_set_name"

      membership <- .abc_safe_join(
        membership,
        term_sets[, c(
          "cluster_annotation_term_set_label",
          "cluster_annotation_term_set_name"
        ), drop = FALSE],
        "cluster_annotation_term_set_label",
        "termset"
      )
    }

    required <- c(
      "cluster_alias",
      "cluster_annotation_term_set_name",
      "cluster_annotation_term_name"
    )

    missing <- setdiff(required, names(membership))
    if (length(missing)) {
      stop(
        taxonomy_dir, ": missing taxonomy field(s): ",
        paste(missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    membership$cluster_alias <- as.character(membership$cluster_alias)

    if (
      anyNA(membership$cluster_alias) ||
      anyNA(membership$cluster_annotation_term_set_name)
    ) {
      stop(taxonomy_dir, ": missing taxonomy keys.", call. = FALSE)
    }

    # One cluster may have only one annotation within a given term set.
    key <- paste(
      membership$cluster_alias,
      membership$cluster_annotation_term_set_name,
      sep = "\r"
    )

    conflicts <- tapply(
      membership$cluster_annotation_term_name,
      key,
      function(x) length(unique(x)) > 1L
    )

    if (any(conflicts)) {
      stop(
        taxonomy_dir, ": conflicting cluster annotations.",
        call. = FALSE
      )
    }
  }

  # Prefer Allen's already-pivoted taxonomy view when one is published.
  if (has("cluster_to_cluster_annotation_membership_pivoted")) {
    cluster_details <- validate_alias(
      get_df("cluster_to_cluster_annotation_membership_pivoted"),
      "pivoted taxonomy"
    )
  } else {
    if (is.null(membership)) {
      stop(
        taxonomy_dir, ": taxonomy membership file missing.",
        call. = FALSE
      )
    }

    detail <- unique(membership[, c(
      "cluster_alias",
      "cluster_annotation_term_set_name",
      "cluster_annotation_term_name"
    )])

    cluster_details <- stats::reshape(
      detail,
      idvar = "cluster_alias",
      timevar = "cluster_annotation_term_set_name",
      direction = "wide"
    )

    names(cluster_details) <- sub(
      "^cluster_annotation_term_name\\.",
      "",
      names(cluster_details)
    )

    cluster_details <- validate_alias(
      cluster_details,
      "taxonomy annotations"
    )
  }

  # Preserve the established R-only convention that absent hierarchy levels
  # are represented as "Other".
  detail_cols <- setdiff(names(cluster_details), "cluster_alias")
  for (column in detail_cols) {
    cluster_details[[column]][is.na(cluster_details[[column]])] <- "Other"
  }

  cluster_colors <- NULL

  if (!is.null(membership) &&
      "color_hex_triplet" %in% names(membership)) {
    color <- unique(membership[, c(
      "cluster_alias",
      "cluster_annotation_term_set_name",
      "color_hex_triplet"
    )])

    color_key <- paste(
      color$cluster_alias,
      color$cluster_annotation_term_set_name,
      sep = "\r"
    )

    color_conflicts <- tapply(
      color$color_hex_triplet,
      color_key,
      function(x) length(unique(x[!is.na(x)])) > 1L
    )

    if (any(color_conflicts)) {
      stop(
        taxonomy_dir, ": conflicting taxonomy colours.",
        call. = FALSE
      )
    }

    cluster_colors <- stats::reshape(
      color,
      idvar = "cluster_alias",
      timevar = "cluster_annotation_term_set_name",
      direction = "wide"
    )

    names(cluster_colors) <- sub(
      "^color_hex_triplet\\.",
      "",
      names(cluster_colors)
    )

    cluster_colors <- validate_alias(
      cluster_colors,
      "taxonomy colours"
    )
  }

  list(
    cluster_details = cluster_details,
    cluster_colors = cluster_colors
  )
}
