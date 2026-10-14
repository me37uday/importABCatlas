dense <- function(m) as.matrix(m)

test_that("expression matrix matches the truth for all cells and genes", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
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
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
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
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
  cells <- rbind(fx$cells[1:2, ], transform(fx$cells[1, ], cell_label = "ghost"))
  rownames(cells) <- cells$cell_label

  m <- suppressMessages(.abc_expression_matrix(cache, cells, fx$genes, fx$genes$gene_symbol))
  expect_equal(colnames(m), c("c1", "c2", "ghost"))
  expect_equal(sum(m[, "ghost"]), 0)
  expect_equal(unname(dense(m[, 1:2])), unname(fx$truth[, 1:2]))
})

test_that("a gene table from another dataset is reported", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
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
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
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
    gene_data = fx$genes, genes = c("Aaa", "Ddd"), offline = TRUE)))

  expect_s4_class(obj, "Seurat")
  expect_equal(colnames(obj), c("c1", "c3", "c6"))
  expect_equal(obj$donor_sex, c(c1 = "M", c3 = "M", c6 = "M"))
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  expect_equal(unname(dense(counts)), unname(fx$truth[c("G1", "G3", "G4"), c("c1", "c3", "c6")]))
})

test_that("fetch_data validates filters", {
  fx <- make_fixture_cache()
  expect_error(fetch_data(download_base = fx$dir, metadata = fx$cells,
                          filters = list(nope = 1), gene_data = fx$genes, offline = TRUE),
               "not found in metadata")
  expect_error(fetch_data(download_base = fx$dir, metadata = fx$cells,
                          filters = list(donor_sex = "X"), gene_data = fx$genes, offline = TRUE),
               "No cells left")
})
