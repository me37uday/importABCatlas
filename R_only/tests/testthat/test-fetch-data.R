dense <- function(m) as.matrix(m)

test_that("expression file planning matches cell routing without downloading", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  files <- .abc_expression_files(
    cache,
    fx$cells,
    data_type = "raw"
  )

  expect_equal(
    files$feature_matrix_label,
    c("TEST-A", "TEST-B")
  )

  expect_equal(files$n_cells, c(3, 3))
  expect_true(all(files$data_type == "raw"))
  expect_true(all(files$size > 0))
  expect_true(all(files$cached))
  expect_equal(files$partial_bytes, c(0, 0))
})

test_that("expression file planning reports resumable partial bytes", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  attrs <- .abc_file_attributes(
    cache,
    "TEST-10X",
    "TEST-A/raw"
  )

  part <- paste0(attrs$local_path, ".part")
  writeBin(as.raw(1:10), part)

  files <- .abc_expression_files(
    cache,
    fx$cells,
    data_type = "raw"
  )

  expect_equal(
    files$partial_bytes[
      files$feature_matrix_label == "TEST-A"
    ],
    10
  )
})

test_that("expression file planning validates routing columns", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  bad <- fx$cells
  bad$feature_matrix_label <- NULL

  expect_error(
    .abc_expression_files(cache, bad),
    "feature_matrix_label"
  )
})

test_that("expression matrix matches the truth for all cells and genes", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(.abc_expression_matrix(cache, cells, fx$genes, fx$genes$gene_symbol))
  expect_s4_class(m, "dgCMatrix")
  expect_equal(rownames(m), c("Aaa", "Bbb", "Aaa.1", "Ddd"))
  expect_equal(colnames(m), cells$cell_label)
  expect_equal(unname(dense(m)), unname(fx$truth))
})

test_that("gene and cell selection, order and block size are respected", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells[c(5, 1, 6, 3), ]          # mixed files, not on-disk order
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(.abc_expression_matrix(cache, cells, fx$genes, c("Aaa", "Ddd"),
                                               block_size = 1L))
  expect_equal(rownames(m), c("Aaa", "Aaa.1", "Ddd"))
  expect_equal(colnames(m), c("c5", "c1", "c6", "c3"))
  expect_equal(unname(dense(m)), unname(fx$truth[c("G1", "G3", "G4"), c("c5", "c1", "c6", "c3")]))
})

test_that("cells absent from their expression file give empty columns", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- rbind(fx$cells[1:2, ], transform(fx$cells[1, ], cell_label = "ghost"))
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(.abc_expression_matrix(cache, cells, fx$genes, fx$genes$gene_symbol))
  expect_equal(colnames(m), c("c1", "c2", "ghost"))
  expect_equal(sum(m[, "ghost"]), 0)
  expect_equal(unname(dense(m[, 1:2])), unname(fx$truth[, 1:2]))
})

test_that("a gene table from another dataset is reported", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells
  rownames(cells) <- cells$cell_label
  other <- transform(fx$genes, gene_identifier = paste0("X", gene_identifier))
  expect_error(suppressMessages(.abc_expression_matrix(cache, cells, other, "Aaa")),
               "missing from")
})

test_that("malformed identifiers in the file fall back to positional alignment", {
  # As in Zeng-Aging-Mouse-10Xv3: a printed pandas Series instead of the ID.
  garbled <- "gene_symbol\nAaa    G1\nAaa    G3\nName: gene_identifier, dtype: object"
  fx <- make_fixture_cache(var_ids = c("G1", "G2", garbled, "G4"))
  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  warns <- capture_warnings(
    m <- suppressMessages(.abc_expression_matrix(cache, cells, fx$genes, "Aaa")))
  expect_length(warns, 2)   # one per expression file
  expect_match(warns, "1 gene identifiers in TEST-[AB]-raw.h5ad are malformed")
  expect_equal(unname(dense(m)), unname(fx$truth[c("G1", "G3"), ]))
  # Genes whose identifiers are intact need no fallback.
  expect_no_warning(suppressMessages(.abc_expression_matrix(cache, cells, fx$genes, "Ddd")))
})

test_that("fetch_data filters cells and returns a Seurat object", {
  fx <- make_fixture_cache()
  obj <- suppressMessages(suppressWarnings(fetch_data(
    download_base = fx$dir, metadata = fx$cells, filters = list(donor_sex = "M"),
    gene_data = fx$genes, genes = c("Aaa", "Ddd"), offline = TRUE,
    release = fx$version)))

  expect_s4_class(obj, "Seurat")
  expect_equal(colnames(obj), c("c1", "c3", "c6"))
  expect_equal(obj$donor_sex, c(c1 = "M", c3 = "M", c6 = "M"))
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  expect_equal(unname(dense(counts)),
               unname(fx$truth[c("G1", "G3", "G4"), c("c1", "c3", "c6")]))
})

test_that("fetch_data validates filters", {
  fx <- make_fixture_cache()
  expect_error(fetch_data(download_base = fx$dir, metadata = fx$cells,
                          filters = list(nope = 1), gene_data = fx$genes,
                          offline = TRUE, release = fx$version),
               "not found in metadata")
  expect_error(fetch_data(download_base = fx$dir, metadata = fx$cells,
                          filters = list(donor_sex = "X"), gene_data = fx$genes,
                          offline = TRUE, release = fx$version),
               "No cells left")
})

test_that("plan_fetch reports the request without downloading expression files", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    filters = list(donor_sex = "M"),
    gene_data = fx$genes,
    genes = c("Aaa", "Ddd"),
    offline = TRUE,
    release = fx$version
  )

  expect_s3_class(plan, "abc_fetch_plan")
  expect_equal(plan$data_type, "raw")
  expect_equal(plan$n_cells, 3)
  expect_equal(plan$n_genes, 3)
  expect_equal(plan$release, fx$version)
  expect_equal(plan$files$feature_matrix_label, c("TEST-A", "TEST-B"))
  expect_equal(plan$files$n_cells, c(2, 1))
  expect_equal(plan$total_download, 0)
})

test_that("plan_fetch reports remaining bytes without downloading", {
  fx <- make_fixture_cache()

  cache <- abc_cache_connect(fx$dir, offline = TRUE, release = fx$version)
  attrs <- .abc_file_attributes(cache, "TEST-10X", "TEST-A/raw")

  unlink(attrs$local_path)
  part <- paste0(attrs$local_path, ".part")
  writeBin(as.raw(1:10), part)

  size_before <- file.info(part)$size

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "Aaa",
    offline = TRUE,
    release = fx$version
  )

  row <- plan$files[plan$files$feature_matrix_label == "TEST-A", ]

  expect_false(row$cached)
  expect_equal(row$partial_bytes, 10)
  expect_equal(row$additional_download, row$size - 10)
  expect_equal(file.info(part)$size, size_before)
  expect_false(file.exists(attrs$local_path))
})

test_that("abc_fetch_plan prints a readable preflight summary", {
  fx <- make_fixture_cache()

  plan <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    filters = list(donor_sex = "M"),
    gene_data = fx$genes,
    genes = c("Aaa", "Ddd"),
    offline = TRUE,
    release = fx$version
  )

  expect_output(print(plan), "<abc_fetch_plan>")
  expect_output(print(plan), "3 genes x 3 cells")
  expect_output(print(plan), "TEST-10X/TEST-A/raw")
  expect_output(print(plan), "additional download")
  expect_output(print(plan), "No expression files were downloaded")
})

test_that("data type resolution uses the dataset default", {
  fx <- make_fixture_cache()

  expect_equal(
    .abc_resolve_data_type(NULL, fx$cells),
    "raw"
  )
})

test_that("data type resolution accepts published representations", {
  fx <- make_fixture_cache()

  expect_equal(.abc_resolve_data_type("raw", fx$cells), "raw")
  expect_equal(.abc_resolve_data_type("log2", fx$cells), "log2")
})

test_that("data type resolution rejects invalid representations", {
  fx <- make_fixture_cache()

  expect_error(
    .abc_resolve_data_type("banana", fx$cells),
    "data_type must be one of"
  )
})

test_that("data type resolution rejects unavailable representations", {
  fx <- make_fixture_cache()
  attr(fx$cells, "abc_default_data_type") <- "log2"
  attr(fx$cells, "abc_available_data_types") <- "log2"

  expect_error(
    .abc_resolve_data_type("raw", fx$cells),
    "does not publish"
  )
})

test_that("plan_fetch routes raw and log2 to the requested representation", {
  fx <- make_fixture_cache()

  raw <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "Aaa",
    data_type = "raw",
    offline = TRUE,
    release = fx$version
  )

  log2 <- plan_fetch(
    download_base = fx$dir,
    metadata = fx$cells,
    gene_data = fx$genes,
    genes = "Aaa",
    data_type = "log2",
    offline = TRUE,
    release = fx$version
  )

  expect_equal(raw$data_type, "raw")
  expect_true(all(raw$files$data_type == "raw"))
  expect_match(raw$files$relative_path, "-raw\\.h5ad$")

  expect_equal(log2$data_type, "log2")
  expect_true(all(log2$files$data_type == "log2"))
  expect_match(log2$files$relative_path, "-log2\\.h5ad$")
})

test_that("expression matrix reads log2 representation without transforming it", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(
    .abc_expression_matrix(
      abc_cache = cache,
      cells = cells,
      gene_data = fx$genes,
      selected_genes = fx$genes$gene_symbol,
      data_type = "log2"
    )
  )

  expect_s4_class(m, "dgCMatrix")
  expect_equal(
    unname(dense(m)),
    unname(fx$truth_log2),
    tolerance = 1e-12
  )
})

test_that("fetch_data stores raw expression in the counts layer", {
  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      gene_data = fx$genes,
      genes = c("Aaa", "Ddd"),
      data_type = "raw",
      offline = TRUE,
      release = fx$version
    )
  )

  counts <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "counts"
  )

  expect_equal(
    unname(dense(counts)),
    unname(fx$truth[c("G1", "G3", "G4"), ])
  )
})

test_that("fetch_data stores log2 expression in the data layer", {
  fx <- make_fixture_cache()

  warns <- capture_warnings(
    obj <- suppressMessages(
      fetch_data(
        download_base = fx$dir,
        metadata = fx$cells,
        gene_data = fx$genes,
        genes = c("Aaa", "Ddd"),
        data_type = "log2",
        offline = TRUE,
        release = fx$version
      )
    )
  )

  expect_length(warns, 1)
  expect_match(
    warns,
    "upstream log2-scale data",
    fixed = TRUE
  )

  data <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "data"
  )

  expect_equal(
    unname(dense(data)),
    unname(fx$truth_log2[c("G1", "G3", "G4"), ]),
    tolerance = 1e-12
  )

  expect_false("counts" %in% SeuratObject::Layers(obj[["RNA"]]))
})

test_that("fetch_data uses the route default representation", {
  fx <- make_fixture_cache()

  obj <- suppressMessages(
    fetch_data(
      download_base = fx$dir,
      metadata = fx$cells,
      gene_data = fx$genes,
      genes = "Ddd",
      offline = TRUE,
      release = fx$version
    )
  )

  expect_true("counts" %in% SeuratObject::Layers(obj[["RNA"]]))
  expect_false("data" %in% SeuratObject::Layers(obj[["RNA"]]))
})

test_that("expression inventory reports published representations", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  inventory <- .abc_expression_inventory(cache, "TEST-10X")

  expect_equal(nrow(inventory), 4L)
  expect_setequal(
    paste(
      inventory$feature_matrix_label,
      inventory$data_type,
      sep = "/"
    ),
    c(
      "TEST-A/raw",
      "TEST-A/log2",
      "TEST-B/raw",
      "TEST-B/log2"
    )
  )
})


test_that("route provenance overrides legacy dataset_label routing", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  cells$dataset_label <- "NOT-AN-EXPRESSION-DIRECTORY"

  files <- .abc_expression_files(
    cache,
    cells,
    data_type = "raw"
  )

  expect_equal(
    files$dataset_label,
    c("TEST-10X", "TEST-10X")
  )
  expect_equal(
    files$feature_matrix_label,
    c("TEST-A", "TEST-B")
  )
  expect_equal(files$n_cells, c(3L, 3L))
})


test_that("route-aware expression routing rejects unknown matrix labels", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  cells$feature_matrix_label[1] <- "OUTSIDE-ROUTE"

  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "Cannot route feature_matrix_label='OUTSIDE-ROUTE'"
  )
})


test_that("route-aware expression routing rejects missing matrix labels", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  cells$feature_matrix_label[1] <- NA_character_

  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "missing feature_matrix_label"
  )
})


test_that("legacy metadata still routes through dataset_label", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(
    fx$dir,
    offline = TRUE,
    release = fx$version
  )

  cells <- fx$cells
  attr(cells, "abc_expression_directories") <- NULL

  files <- .abc_expression_files(
    cache,
    cells,
    data_type = "raw"
  )

  expect_equal(
    files$dataset_label,
    c("TEST-10X", "TEST-10X")
  )
  expect_equal(
    files$feature_matrix_label,
    c("TEST-A", "TEST-B")
  )
})

test_that("route-aware expression routing rejects ambiguous matrices", {
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

  manifest$file_listing[["TEST-ALT"]] <-
    manifest$file_listing[["TEST-10X"]]

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

  cells <- fx$cells
  attr(cells, "abc_expression_directories") <-
    c("TEST-10X", "TEST-ALT")

  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "Ambiguous expression routing"
  )
})

test_that("expression files outside declared route are not accessible", {
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

  outside <- manifest$file_listing[["TEST-10X"]]
  names(outside$expression_matrices) <-
    paste0(
      names(outside$expression_matrices),
      "-OUTSIDE"
    )

  manifest$file_listing[["TEST-OUTSIDE"]] <- outside

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

  cells <- fx$cells[1, , drop = FALSE]
  cells$feature_matrix_label <- "TEST-A-OUTSIDE"

  # TEST-OUTSIDE really contains this matrix, but this route only permits
  # TEST-10X. Expression discovery must not escape the declared route.
  expect_error(
    .abc_expression_files(
      cache,
      cells,
      data_type = "raw"
    ),
    "Cannot route feature_matrix_label='TEST-A-OUTSIDE'"
  )
})
