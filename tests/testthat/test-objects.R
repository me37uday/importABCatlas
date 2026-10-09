# These tests are supplied but were not run in the delivery environment (no R).
# Python numerical dependencies must be installed before running them.
test_that("Python bridge preserves a synthetic sparse matrix and enormous IDs", {
  x <- example_data("matrix")
  expect_s4_class(x$expression, "sparseMatrix")
  expect_equal(dim(x$expression), c(4L, 3L))
  expect_identical(colnames(x$expression)[1], "182941331246012878296807398333956011710")
  expect_equal(unname(as.matrix(x$expression)), t(matrix(c(1,0,3,0, 0,2,0,4, 5,6,0,0), nrow=3, byrow=TRUE)))
  expect_identical(rownames(x$gene_data), rownames(x$expression))
  expect_true(x$provenance$synthetic_fixture)
})

test_that("Seurat preserves exact axes and feature metadata", {
  skip_if_not_installed("SeuratObject")
  x <- example_data("seurat")
  expect_s4_class(x, "Seurat")
  expect_equal(dim(x), c(4L, 3L))
  expect_true(x@misc$importABCatlas$synthetic_fixture)
  expect_identical(unname(x$cell_label), unname(colnames(x)))
})

test_that("SCE and SpatialExperiment retain matrix and anatomical coordinates", {
  skip_if_not_installed("SingleCellExperiment")
  x <- example_data("sce")
  expect_s4_class(x, "SingleCellExperiment")
  expect_equal(dim(x), c(4L, 3L))
  expect_identical(SummarizedExperiment::assayNames(x), "counts")
  skip_if_not_installed("SpatialExperiment")
  y <- example_data("spatial")
  expect_equal(unname(SpatialExperiment::spatialCoords(y)), cbind(c(1,2,3),c(2,3,4)))
})
