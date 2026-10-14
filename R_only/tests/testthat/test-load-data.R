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
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
  msgs <- capture_messages(res <- .build_dataset(fixture_spec, cache, "TEST"))
  expect_match(msgs, "TEST: 6 cells", all = FALSE)
  expect_match(msgs, "TEST: 4 genes", all = FALSE)

  cm <- res$cell_metadata
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
  res <- suppressMessages(load_data("TEST", download_base = fx$dir, offline = TRUE))
  expect_named(res, c("cell_metadata", "gene_data", "unique_values"))
  expect_equal(nrow(res$cell_metadata), 6)
})
