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
