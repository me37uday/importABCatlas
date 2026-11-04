test_that("load_data validates the dataset name", {
  expect_error(load_data(), "Specify a dataset")
  expect_error(load_data("NotADataset"), "should be one of")
})

test_that("every registry entry has the fields the builder needs", {
  for (name in names(.abc_datasets)) {
    spec <- .abc_datasets[[name]]
    expect_true(is.character(spec$cell$dir), info = name)
    expect_true(is.character(spec$gene$dir), info = name)
    expect_true(is.character(spec$cols_of_interest), info = name)
    for (j in spec$joins) expect_true(all(c("dir", "file", "key") %in% names(j)), info = name)
  }
})

test_that(".build_dataset assembles cells, taxonomy, ROI colours and genes", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  msgs <- capture_messages(res <- .build_dataset(fixture_spec, cache, "TEST"))
  expect_match(msgs, "TEST: 6 cells", all = FALSE)
  expect_match(msgs, "TEST: 4 genes", all = FALSE)

  cm <- res$cell_metadata

  expect_identical(attr(cm, "abc_dataset"), "TEST")
  expect_identical(attr(cm, "abc_default_data_type"), "raw")
  expect_identical(
    attr(cm, "abc_available_data_types"),
    c("raw", "log2")
  )
  expect_identical(
    attr(cm, "abc_expression_directories"),
    "TEST-10X"
  )
  expect_identical(
    attr(cm, "abc_modality"),
    "single_cell_transcriptomics"
  )
  expect_identical(attr(cm, "abc_spatial"), FALSE)

  expect_equal(nrow(cm), 6)
  expect_equal(cm$cell_label, fx$cells$cell_label)
  # Taxonomy columns follow the term-set order; missing terms become "Other".
  expect_true(all(c("supercluster", "cluster", "supercluster_color", "cluster_color")
                  %in% names(cm)))
  expect_equal(cm$cluster, c("cl1", "cl2", "cl1", "cl2", "cl2", "cl1"))
  expect_equal(cm$supercluster, c("sup1", "Other", "sup1", "Other", "Other", "sup1"))
  expect_equal(cm$cluster_color[1:2], c("#aa0000", "#cc0000"))
  expect_equal(cm$region_of_interest_color, c("#111111", "#222222", "#111111",
                                              "#222222", "#111111", "#111111"))

  expect_equal(rownames(res$gene_data), fx$genes$gene_identifier)
  expect_equal(res$gene_data$name[1], "alpha ")

  expect_equal(names(res$unique_values), c("donor_sex", "cluster", "supercluster"))
  expect_setequal(res$unique_values$donor_sex, c("M", "F", NA))
})

test_that("load_data runs end to end on a registered dataset (offline)", {
  fx <- make_fixture_cache()
  local_mocked_bindings(.abc_datasets = list(TEST = fixture_spec))
  res <- suppressMessages(load_data(
   "TEST",
   download_base = fx$dir,
   offline = TRUE,
   release = fx$version
  ))
  expect_named(res, c("cell_metadata", "gene_data", "unique_values"))
  expect_equal(nrow(res$cell_metadata), 6)
})

test_that("registered datasets declare expression routing metadata", {
  allowed_modalities <- c(
    "single_cell_transcriptomics",
    "single_nucleus_transcriptomics",
    "spatial_transcriptomics",
    "imputed_spatial_transcriptomics",
    "patch_seq_transcriptomics"
  )

  for (dataset in names(.abc_datasets)) {
    spec <- .abc_datasets[[dataset]]

    expect_true(
      length(spec$expression_directories) >= 1L,
      info = dataset
    )
    expect_true(
      all(nzchar(spec$expression_directories)),
      info = dataset
    )
    expect_true(
      spec$modality %in% allowed_modalities,
      info = dataset
    )
    expect_true(
      is.logical(spec$spatial) && length(spec$spatial) == 1L,
      info = dataset
    )
    expect_true(
      length(spec$available_data_types) >= 1L,
      info = dataset
    )
    expect_true(
      all(spec$available_data_types %in% c("raw", "log2")),
      info = dataset
    )
    expect_true(
      spec$default_data_type %in% spec$available_data_types,
      info = dataset
    )
  }
})

test_that("taxonomy membership supplies the authoritative cluster alias", {
  fx <- make_fixture_cache()

  # Add cell_to_cluster_membership to the taxonomy directory. Deliberately
  # reverse the fixture's original aliases so the test proves that membership,
  # rather than cell_metadata$cluster_alias, drives taxonomy annotation.
  membership <- data.frame(
    cell_label = fx$cells$cell_label,
    cluster_alias = ifelse(fx$cells$cluster_alias == 1L, "2", "1"),
    stringsAsFactors = FALSE
  )

  rel <- file.path(
    "metadata",
    "TEST-taxonomy",
    fx$version,
    "cell_to_cluster_membership.csv"
  )
  path <- file.path(fx$dir, rel)

  utils::write.csv(
    membership,
    path,
    row.names = FALSE,
    na = ""
  )

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  entry <- list(
    version = fx$version,
    relative_path = rel,
    url = paste0("https://example.invalid/", rel),
    size = file.size(path),
    file_hash = unname(tools::md5sum(path))
  )

  manifest$file_listing[["TEST-taxonomy"]][["metadata"]][["cell_to_cluster_membership"]] <- list(
    files = stats::setNames(list(entry), tools::file_ext(rel))
  )

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  res <- suppressMessages(
    .build_dataset(
      fixture_spec,
      cache,
      "TEST"
    )
  )

  cm <- res$cell_metadata

  expect_equal(nrow(cm), nrow(fx$cells))
  expect_equal(cm$cell_label, fx$cells$cell_label)

  # Existing aliases are retained for provenance.
  expect_equal(
    cm$cluster_alias_original,
    as.character(fx$cells$cluster_alias)
  )

  # Membership is authoritative.
  expect_equal(
    cm$cluster_alias,
    membership$cluster_alias
  )

  # Taxonomy annotations follow the authoritative membership alias.
  expect_equal(
    cm$cluster,
    c("cl2", "cl1", "cl2", "cl1", "cl1", "cl2")
  )
  expect_equal(
    cm$supercluster,
    c("Other", "sup1", "Other", "sup1", "sup1", "Other")
  )
})

test_that("safe metadata joins reject duplicate right-hand keys", {
  left <- data.frame(
    cell_label = c("c1", "c2"),
    value = c(1, 2)
  )

  right <- data.frame(
    cell_label = c("c1", "c1"),
    annotation = c("A", "B")
  )

  expect_error(
    .abc_safe_join(
      left,
      right,
      "cell_label",
      "membership"
    ),
    "right-hand cell_label is not a non-missing unique key"
  )
})

test_that("taxonomy prefers an upstream pivoted annotation view", {
  fx <- make_fixture_cache()

  pivoted <- data.frame(
    cluster_alias = c("1", "2"),
    supercluster = c("pivot_sup1", "pivot_sup2"),
    cluster = c("pivot_cl1", "pivot_cl2"),
    stringsAsFactors = FALSE
  )

  rel <- file.path(
    "metadata",
    "TEST-taxonomy",
    fx$version,
    "cluster_to_cluster_annotation_membership_pivoted.csv"
  )
  path <- file.path(fx$dir, rel)

  utils::write.csv(
    pivoted,
    path,
    row.names = FALSE,
    na = ""
  )

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  entry <- list(
    version = fx$version,
    relative_path = rel,
    url = paste0("https://example.invalid/", rel),
    size = file.size(path),
    file_hash = unname(tools::md5sum(path))
  )

  manifest$file_listing[["TEST-taxonomy"]][["metadata"]][["cluster_to_cluster_annotation_membership_pivoted"]] <- list(
    files = stats::setNames(list(entry), tools::file_ext(rel))
  )

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  tax <- .build_taxonomy(cache, "TEST-taxonomy")

  expect_equal(
    tax$cluster_details$cluster_alias,
    c("1", "2")
  )
  expect_equal(
    tax$cluster_details$supercluster,
    c("pivot_sup1", "pivot_sup2")
  )
  expect_equal(
    tax$cluster_details$cluster,
    c("pivot_cl1", "pivot_cl2")
  )

  # Colours still come from the long membership table.
  expect_equal(
    tax$cluster_colors$cluster,
    c("#aa0000", "#cc0000")
  )
})

test_that("taxonomy reconstructs names from term tables when membership has labels only", {
  fx <- make_fixture_cache()

  membership_rel <- file.path(
    "metadata",
    "TEST-taxonomy",
    fx$version,
    "cluster_to_cluster_annotation_membership.csv"
  )
  membership_path <- file.path(fx$dir, membership_rel)

  membership <- utils::read.csv(
    membership_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Replace embedded names with stable labels, matching taxonomy releases
  # where human-readable names must be recovered from the term tables.
  membership$cluster_annotation_term_label <- c("T1", "T2", "T3", "T4")
  membership$cluster_annotation_term_set_label <- c("L2", "L1", "L2", "L1")

  membership$cluster_annotation_term_name <- NULL
  membership$cluster_annotation_term_set_name <- NULL

  utils::write.csv(
    membership,
    membership_path,
    row.names = FALSE,
    na = ""
  )

  term <- data.frame(
    label = c("T1", "T2", "T3", "T4"),
    name = c("cl1", "sup1", "cl2", NA),
    stringsAsFactors = FALSE
  )

  term_rel <- file.path(
    "metadata",
    "TEST-taxonomy",
    fx$version,
    "cluster_annotation_term.csv"
  )
  term_path <- file.path(fx$dir, term_rel)

  utils::write.csv(
    term,
    term_path,
    row.names = FALSE,
    na = ""
  )

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  # The existing membership file changed, so refresh its manifest metadata.
  membership_entry <-
    manifest$file_listing[["TEST-taxonomy"]][["metadata"]][["cluster_to_cluster_annotation_membership"]][["files"]][["csv"]]

  membership_entry$size <- file.size(membership_path)
  membership_entry$file_hash <-
    unname(tools::md5sum(membership_path))

  manifest$file_listing[["TEST-taxonomy"]][["metadata"]][["cluster_to_cluster_annotation_membership"]][["files"]][["csv"]] <- membership_entry

  term_entry <- list(
    version = fx$version,
    relative_path = term_rel,
    url = paste0("https://example.invalid/", term_rel),
    size = file.size(term_path),
    file_hash = unname(tools::md5sum(term_path))
  )

  manifest$file_listing[["TEST-taxonomy"]][["metadata"]][["cluster_annotation_term"]] <- list(
    files = stats::setNames(list(term_entry), tools::file_ext(term_rel))
  )

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  tax <- .build_taxonomy(cache, "TEST-taxonomy")

  expect_equal(
    tax$cluster_details$cluster,
    c("cl1", "cl2")
  )
  expect_equal(
    tax$cluster_details$supercluster,
    c("sup1", "Other")
  )
  expect_equal(
    tax$cluster_colors$cluster,
    c("#aa0000", "#cc0000")
  )
})

test_that("taxonomy rejects conflicting annotations for one hierarchy level", {
  fx <- make_fixture_cache()

  rel <- file.path(
    "metadata",
    "TEST-taxonomy",
    fx$version,
    "cluster_to_cluster_annotation_membership.csv"
  )
  path <- file.path(fx$dir, rel)

  membership <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  conflict <- membership[1, , drop = FALSE]
  conflict$cluster_annotation_term_name <- "different-cluster"
  membership <- rbind(membership, conflict)

  utils::write.csv(
    membership,
    path,
    row.names = FALSE,
    na = ""
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  expect_error(
    .build_taxonomy(cache, "TEST-taxonomy"),
    "conflicting cluster annotations"
  )
})

test_that("sequential annotation joins preserve authoritative columns and provenance", {
  cells <- data.frame(
    cell_label = c("c1", "c2", "c3"),
    cluster_alias = c("10", "20", "30"),
    wmb_cluster_alias = c(100L, 200L, NA_integer_),
    stringsAsFactors = FALSE
  )

  native <- data.frame(
    cell_label = c("c1", "c2", "c3"),
    cluster_alias = c("10", "20", "30"),
    cluster_name = c("native-1", "native-2", "native-3"),
    cluster_color = c("#111111", "#222222", "#333333"),
    stringsAsFactors = FALSE
  )

  cross <- data.frame(
    cell_label = c("c1", "c2"),
    cluster_alias = c("10", "20"),
    cluster_name = c("native-1", "native-2"),
    cluster_color = c("#111111", "#222222"),
    wmb_cluster_alias = c(100L, 200L),
    wmb_cluster_name = c("wmb-1", "wmb-2"),
    stringsAsFactors = FALSE
  )

  out <- .abc_safe_join(
    cells,
    native,
    "cell_label",
    "cell_cluster_annotations"
  )

  out <- .abc_safe_join(
    out,
    cross,
    "cell_label",
    "cell_cross_mapping_annotations"
  )

  expect_equal(out$cluster_alias, c("10", "20", "30"))

  expect_equal(
    out$cluster_alias__cell_cluster_annotations,
    c("10", "20", "30")
  )

  expect_equal(
    out$cluster_name,
    c("native-1", "native-2", "native-3")
  )

  expect_equal(
    out$cluster_name__cell_cross_mapping_annotations,
    c("native-1", "native-2", NA)
  )

  expect_equal(
    out$wmb_cluster_alias,
    c(100L, 200L, NA_integer_)
  )

  expect_equal(
    out$wmb_cluster_alias__cell_cross_mapping_annotations,
    c(100L, 200L, NA_integer_)
  )

  expect_equal(
    out$wmb_cluster_name,
    c("wmb-1", "wmb-2", NA)
  )

  expect_equal(out$cell_label, c("c1", "c2", "c3"))
})

test_that("PMDBS uses its native gene universe", {
  expect_identical(
    .abc_datasets$PMDBS$gene$dir,
    "ASAP-PMDBS-10X"
  )

  expect_identical(
    .abc_datasets$PMDBS$gene$file %||% "gene",
    "gene"
  )
})

test_that("Batch 1 transcriptomic routes declare reference-compatible sources", {
  aibs <- .abc_datasets$ConsensusMouseAIBS

  expect_identical(
    aibs$cell$dir,
    "Consensus-WMB-AIBS-10X"
  )
  expect_identical(
    aibs$gene$dir,
    "WMB-10X"
  )
  expect_identical(
    aibs$taxonomy,
    "Consensus-WMB-integrated-taxonomy"
  )
  expect_identical(
    aibs$expression_directories,
    c("WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti")
  )
  expect_identical(
    aibs$modality,
    "single_cell_transcriptomics"
  )

  macosko <- .abc_datasets$ConsensusMouseMacosko

  expect_identical(
    macosko$cell$dir,
    "Consensus-WMB-Macosko-10X"
  )
  expect_identical(
    macosko$gene$dir,
    "Consensus-WMB-Macosko-10X"
  )
  expect_identical(
    macosko$taxonomy,
    "Consensus-WMB-integrated-taxonomy"
  )
  expect_identical(
    macosko$expression_directories,
    "Consensus-WMB-Macosko-10X"
  )
  expect_identical(
    macosko$modality,
    "single_nucleus_transcriptomics"
  )

  developing <- .abc_datasets$DevelopingMouse

  expect_identical(
    developing$cell$dir,
    "Developing-Mouse-Vis-Cortex-10X"
  )
  expect_identical(
    developing$gene$dir,
    "Developing-Mouse-Vis-Cortex-10X"
  )
  expect_identical(
    developing$taxonomy,
    "Developing-Mouse-Vis-Cortex-taxonomy"
  )
  expect_identical(
    developing$expression_directories,
    "Developing-Mouse-Vis-Cortex-10X"
  )
  expect_identical(
    developing$modality,
    "single_cell_transcriptomics"
  )
})

test_that("gene normalization leaves unique identifiers unchanged", {
  gene <- data.frame(
    gene_identifier = c("G1", "G2", "G3"),
    gene_symbol = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  out <- .abc_normalize_gene_metadata(gene)

  expect_s3_class(out, "data.frame")
  expect_identical(out$gene_identifier, c("G1", "G2", "G3"))
  expect_false("source_gene_identifier" %in% names(out))
})

test_that("gene normalization preserves duplicate source occurrence", {
  gene <- data.frame(
    gene_identifier = c("G1", "G2", "G1", "G1"),
    gene_symbol = c("A1", "B", "A2", "A3"),
    stringsAsFactors = FALSE
  )

  out <- .abc_normalize_gene_metadata(gene)

  expect_s3_class(out, "data.frame")
  expect_identical(
    out$gene_identifier,
    c("G1.iabdup1", "G2", "G1.iabdup2", "G1.iabdup3")
  )
  expect_identical(
    out$source_gene_identifier,
    c("G1", "G2", "G1", "G1")
  )
  expect_identical(
    out$source_gene_occurrence,
    c(1L, 1L, 2L, 3L)
  )
  expect_identical(
    out$source_gene_group_size,
    c(3L, 1L, 3L, 3L)
  )
  expect_identical(anyDuplicated(out$gene_identifier), 0L)
})

test_that("gene normalization avoids collisions with published identifiers", {
  gene <- data.frame(
    gene_identifier = c("G1", "G1.iabdup1", "G1"),
    gene_symbol = c("A1", "existing", "A2"),
    stringsAsFactors = FALSE
  )

  out <- .abc_normalize_gene_metadata(gene)

  expect_s3_class(out, "data.frame")
  expect_identical(
    out$gene_identifier,
    c("G1.iabdup2", "G1.iabdup1", "G1.iabdup3")
  )
  expect_identical(
    out$source_gene_occurrence,
    c(1L, 1L, 2L)
  )
  expect_identical(
    out$source_gene_group_size,
    c(2L, 1L, 2L)
  )
  expect_identical(anyDuplicated(out$gene_identifier), 0L)
})

test_that("HMBA species routes declare automatic species-specific sources", {
  expected <- list(
    HMBA_Human = list(
      species = "Human",
      gene = "human_gene"
    ),
    HMBA_Macaque = list(
      species = "Macaque",
      gene = "macaque_gene"
    ),
    HMBA_Marmoset = list(
      species = "Marmoset",
      gene = "marmoset_gene"
    )
  )

  for (route in names(expected)) {
    spec <- .abc_datasets[[route]]
    want <- expected[[route]]

    expect_identical(spec$cell$dir, "HMBA-10xMultiome-BG")
    expect_identical(
      spec$expression_directories,
      "HMBA-10xMultiome-BG"
    )
    expect_identical(
      spec$taxonomy,
      "HMBA-BG-taxonomy-CCN20250428"
    )
    expect_identical(spec$gene$dir, "HMBA-10xMultiome-BG")
    expect_identical(spec$gene$file, want$gene)

    expect_identical(spec$subset$dir, "HMBA-10xMultiome-BG")
    expect_identical(spec$subset$file, "donor")
    expect_identical(spec$subset$key, "donor_label")
    expect_identical(spec$subset$column, "species_genus")
    expect_identical(spec$subset$value, want$species)

    expect_identical(
      spec$modality,
      "single_nucleus_transcriptomics"
    )
    expect_false(spec$spatial)
  }
})

test_that("route subsets are resolved internally before annotation", {
  fx <- make_fixture_cache()

  # Give fixture cells donor identities without changing the shared helper.
  cell_path <- file.path(
    fx$dir,
    "metadata",
    "TEST-10X",
    fx$version,
    "cell_metadata.csv"
  )

  cells <- utils::read.csv(
    cell_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cells$donor_label <- c("d1", "d2", "d1", "d2", "d3", "d3")

  utils::write.csv(
    cells,
    cell_path,
    row.names = FALSE,
    na = ""
  )

  donor <- data.frame(
    donor_label = c("d1", "d2", "d3"),
    species_genus = c("A", "B", "A"),
    donor_age = c("young", "old", "middle"),
    stringsAsFactors = FALSE
  )

  donor_rel <- file.path(
    "metadata",
    "TEST-10X",
    fx$version,
    "donor.csv"
  )

  donor_path <- file.path(fx$dir, donor_rel)

  utils::write.csv(
    donor,
    donor_path,
    row.names = FALSE,
    na = ""
  )

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  # Refresh cell_metadata because the fixture file was modified.
  cell_entry <-
    manifest$file_listing[["TEST-10X"]][["metadata"]][["cell_metadata"]][["files"]][["csv"]]

  cell_entry$size <- file.size(cell_path)
  cell_entry$file_hash <- unname(tools::md5sum(cell_path))

  manifest$file_listing[["TEST-10X"]][["metadata"]][["cell_metadata"]][["files"]][["csv"]] <-
    cell_entry

  donor_entry <- list(
    version = fx$version,
    relative_path = donor_rel,
    url = paste0("https://example.invalid/", donor_rel),
    size = file.size(donor_path),
    file_hash = unname(tools::md5sum(donor_path))
  )

  manifest$file_listing[["TEST-10X"]][["metadata"]][["donor"]] <- list(
    files = stats::setNames(
      list(donor_entry),
      tools::file_ext(donor_rel)
    )
  )

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  spec <- fixture_spec
  spec$subset <- list(
    dir = "TEST-10X",
    file = "donor",
    key = "donor_label",
    column = "species_genus",
    value = "A"
  )

  msgs <- capture_messages(
    out <- .build_dataset(
      spec,
      cache,
      "TEST_A"
    )
  )

  cm <- out$cell_metadata

  expect_match(
    msgs,
    "TEST_A: 4 cells",
    all = FALSE
  )

  expect_identical(
    cm$cell_label,
    c("c1", "c3", "c5", "c6")
  )

  expect_true(all(cm$species_genus == "A"))

  expect_identical(
    cm$donor_label,
    c("d1", "d1", "d3", "d3")
  )

  expect_identical(
    cm$donor_age,
    c("young", "young", "middle", "middle")
  )

  expect_equal(nrow(cm), 4L)
})

test_that("SEA-AD routes declare reference-compatible sources", {
  expected <- list(
    SEAAD_Multiregion = list(
      directory = "SEA-AD-Multiregion-10X",
      taxonomy = "SEA-AD-Multiregion-taxonomy"
    ),
    SEAAD_CaH = list(
      directory = "SEA-AD-CaH-10X",
      taxonomy = "SEA-AD-CaH-taxonomy"
    )
  )

  for (route in names(expected)) {
    spec <- .abc_datasets[[route]]
    want <- expected[[route]]

    expect_identical(spec$cell$dir, want$directory)
    expect_identical(spec$gene$dir, want$directory)
    expect_identical(
      spec$expression_directories,
      want$directory
    )
    expect_identical(spec$taxonomy, want$taxonomy)
    expect_true(spec$auto_joins)
    expect_identical(
      spec$default_data_type,
      "raw"
    )
    expect_identical(
      spec$available_data_types,
      c("raw", "log2")
    )
    expect_identical(
      spec$modality,
      "single_nucleus_transcriptomics"
    )
    expect_false(spec$spatial)
  }
})

test_that("standard metadata auto-joins attach available route annotations", {
  fx <- make_fixture_cache()

  cell_path <- file.path(
    fx$dir,
    "metadata",
    "TEST-10X",
    fx$version,
    "cell_metadata.csv"
  )

  cells <- utils::read.csv(
    cell_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cells$library_label <- c("l1", "l1", "l2", "l2", "l3", "l3")
  cells$donor_label <- c("d1", "d1", "d2", "d2", "d3", "d3")

  utils::write.csv(
    cells,
    cell_path,
    row.names = FALSE,
    na = ""
  )

  library_table <- data.frame(
    library_label = c("l1", "l2", "l3"),
    region = c("R1", "R2", "R3"),
    stringsAsFactors = FALSE
  )

  donor <- data.frame(
    donor_label = c("d1", "d2", "d3"),
    donor_age = c("young", "middle", "old"),
    diagnosis = c("control", "case", "case"),
    stringsAsFactors = FALSE
  )

  disease <- data.frame(
    donor_label = c("d1", "d2", "d3"),
    diagnosis = c("none", "AD", "AD"),
    pathology = c("low", "medium", "high"),
    stringsAsFactors = FALSE
  )

  tables <- list(
    library = library_table,
    donor = donor,
    disease = disease
  )

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  # Refresh the modified cell table entry.
  cell_entry <-
    manifest$file_listing[["TEST-10X"]][["metadata"]][["cell_metadata"]][["files"]][["csv"]]

  cell_entry$size <- file.size(cell_path)
  cell_entry$file_hash <- unname(tools::md5sum(cell_path))

  manifest$file_listing[["TEST-10X"]][["metadata"]][["cell_metadata"]][["files"]][["csv"]] <-
    cell_entry

  for (name in names(tables)) {
    rel <- file.path(
      "metadata",
      "TEST-10X",
      fx$version,
      paste0(name, ".csv")
    )

    path <- file.path(fx$dir, rel)

    utils::write.csv(
      tables[[name]],
      path,
      row.names = FALSE,
      na = ""
    )

    entry <- list(
      version = fx$version,
      relative_path = rel,
      url = paste0("https://example.invalid/", rel),
      size = file.size(path),
      file_hash = unname(tools::md5sum(path))
    )

    manifest$file_listing[["TEST-10X"]][["metadata"]][[name]] <- list(
      files = stats::setNames(
        list(entry),
        tools::file_ext(rel)
      )
    )
  }

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  spec <- fixture_spec
  spec$auto_joins <- TRUE

  out <- suppressMessages(
    .build_dataset(
      spec,
      cache,
      "TEST_AUTO"
    )
  )

  cm <- out$cell_metadata

  expect_identical(
    cm$region,
    c("R1", "R1", "R2", "R2", "R3", "R3")
  )

  expect_identical(
    cm$donor_age,
    c("young", "young", "middle", "middle", "old", "old")
  )

  # Donor metadata arrives first and remains authoritative.
  expect_identical(
    cm$diagnosis,
    c("control", "control", "case", "case", "case", "case")
  )

  # The overlapping disease value is preserved with provenance.
  expect_identical(
    cm$diagnosis__disease,
    c("none", "none", "AD", "AD", "AD", "AD")
  )

  expect_identical(
    cm$pathology,
    c("low", "low", "medium", "medium", "high", "high")
  )

  expect_equal(nrow(cm), 6L)
})

test_that("MERFISH routes declare measured and imputed expression correctly", {
  measured <- .abc_datasets$MERFISH

  expect_identical(
    measured$cell$dir,
    "MERFISH-C57BL6J-638850"
  )
  expect_identical(
    measured$gene$dir,
    "MERFISH-C57BL6J-638850"
  )
  expect_identical(
    measured$expression_directories,
    "MERFISH-C57BL6J-638850"
  )
  expect_identical(measured$taxonomy, "WMB-taxonomy")
  expect_identical(measured$default_data_type, "raw")
  expect_identical(
    measured$available_data_types,
    c("raw", "log2")
  )
  expect_identical(
    measured$modality,
    "spatial_transcriptomics"
  )
  expect_true(measured$spatial)

  imputed <- .abc_datasets$MERFISH_imputed

  expect_identical(
    imputed$cell$dir,
    "MERFISH-C57BL6J-638850"
  )
  expect_identical(imputed$gene$dir, "WMB-10X")
  expect_identical(
    imputed$expression_directories,
    "MERFISH-C57BL6J-638850-imputed"
  )
  expect_identical(
    imputed$expression_matrix_map,
    c("C57BL6J-638850" = "C57BL6J-638850-imputed")
  )
  expect_identical(imputed$taxonomy, "WMB-taxonomy")
  expect_identical(imputed$default_data_type, "log2")
  expect_identical(imputed$available_data_types, "log2")
  expect_identical(
    imputed$modality,
    "imputed_spatial_transcriptomics"
  )
  expect_true(imputed$spatial)
})

test_that("Zhuang spatial routes declare native sources", {
  expected <- c(
    Zhuang1 = "Zhuang-ABCA-1",
    Zhuang2 = "Zhuang-ABCA-2",
    Zhuang3 = "Zhuang-ABCA-3",
    Zhuang4 = "Zhuang-ABCA-4"
  )

  for (route in names(expected)) {
    directory <- unname(expected[[route]])
    spec <- .abc_datasets[[route]]

    expect_identical(spec$cell$dir, directory)
    expect_identical(spec$gene$dir, directory)
    expect_identical(
      spec$expression_directories,
      directory
    )
    expect_identical(spec$taxonomy, "WMB-taxonomy")

    expect_identical(
      spec$default_data_type,
      "raw"
    )
    expect_identical(
      spec$available_data_types,
      c("raw", "log2")
    )
    expect_identical(
      spec$modality,
      "spatial_transcriptomics"
    )
    expect_true(spec$spatial)

    # These routes already publish x/y/z in cell_metadata;
    # no separate spatial metadata join is required.
    expect_null(spec$joins)
  }
})

test_that("HMBA spatial routes declare native spatial sources", {
  expected <- c(
    HMBA_spatial_Human = "HMBA-MERSCOPE-H22.30.001-BG",
    HMBA_spatial_Macaque = "HMBA-MERSCOPE-QM23.50.001-BG",
    HMBA_spatial_Marmoset = "HMBA-Xenium-CJ23.56.004-BG"
  )

  for (route in names(expected)) {
    directory <- unname(expected[[route]])
    spec <- .abc_datasets[[route]]

    expect_identical(spec$cell$dir, directory)
    expect_identical(spec$gene$dir, directory)
    expect_identical(spec$expression_directories, directory)

    expect_identical(
      spec$taxonomy,
      "HMBA-BG-taxonomy-CCN20250428"
    )

    expect_identical(spec$default_data_type, "raw")
    expect_identical(
      spec$available_data_types,
      c("raw", "log2")
    )
    expect_identical(
      spec$modality,
      "spatial_transcriptomics"
    )
    expect_true(spec$spatial)
  }
})

test_that("spatial routes automatically attach available spatial metadata", {
  fx <- make_fixture_cache()

  manifest_path <- file.path(
    fx$dir,
    "releases",
    fx$version,
    "manifest.json"
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  tables <- list(
    slab_plane_coordinates = data.frame(
      cell_label = fx$cells$cell_label,
      slab_plane_label = rep("plane1", nrow(fx$cells)),
      x_slab_mm = seq_len(nrow(fx$cells)),
      y_slab_mm = seq_len(nrow(fx$cells)) + 10,
      stringsAsFactors = FALSE
    ),
    cell_anatomical_annotations = data.frame(
      cell_label = fx$cells$cell_label,
      parcellation_term_symbol = rep("STR", nrow(fx$cells)),
      stringsAsFactors = FALSE
    )
  )

  for (name in names(tables)) {
    rel <- file.path(
      "metadata",
      "TEST-10X",
      fx$version,
      paste0(name, ".csv")
    )

    path <- file.path(fx$dir, rel)

    utils::write.csv(
      tables[[name]],
      path,
      row.names = FALSE,
      na = ""
    )

    entry <- list(
      version = fx$version,
      relative_path = rel,
      url = paste0("https://example.invalid/", rel),
      size = file.size(path),
      file_hash = unname(tools::md5sum(path))
    )

    manifest$file_listing[["TEST-10X"]][["metadata"]][[name]] <- list(
      files = stats::setNames(
        list(entry),
        tools::file_ext(rel)
      )
    )
  }

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE
  )

  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  spec <- fixture_spec
  spec$spatial <- TRUE

  out <- suppressMessages(
    .build_dataset(
      spec,
      cache,
      "TEST_SPATIAL"
    )
  )

  cm <- out$cell_metadata

  expect_identical(
    cm$slab_plane_label,
    rep("plane1", nrow(fx$cells))
  )
  expect_identical(
    cm$x_slab_mm,
    as.integer(seq_len(nrow(fx$cells)))
  )
  expect_identical(
    cm$parcellation_term_symbol,
    rep("STR", nrow(fx$cells))
  )
  expect_equal(nrow(cm), nrow(fx$cells))
})

test_that("HMBA PatchSeq declares native log2-only sources", {
  spec <- .abc_datasets$HMBA_PatchSeq

  expect_identical(
    spec$cell$dir,
    "HMBA-Macaque-PatchSeq"
  )
  expect_identical(
    spec$gene$dir,
    "HMBA-Macaque-PatchSeq"
  )
  expect_identical(
    spec$expression_directories,
    "HMBA-Macaque-PatchSeq"
  )
  expect_identical(
    spec$taxonomy,
    "HMBA-BG-taxonomy-CCN20250428"
  )

  expect_true(spec$auto_joins)
  expect_identical(spec$default_data_type, "log2")
  expect_identical(spec$available_data_types, "log2")
  expect_identical(
    spec$modality,
    "patch_seq_transcriptomics"
  )
  expect_false(spec$spatial)
})
