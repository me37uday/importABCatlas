test_that("list_datasets exposes every native route exactly once", {
  out <- list_datasets()

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 23L)

  expect_identical(
    out$dataset,
    names(.abc_datasets)
  )

  expect_false(anyDuplicated(out$dataset) > 0L)

  expect_setequal(
    out$dataset,
    names(.abc_datasets)
  )
})

test_that("list_datasets exposes route declarations faithfully", {
  out <- list_datasets()

  expect_identical(
    names(out),
    c(
      "dataset",
      "modality",
      "spatial",
      "imputed",
      "default_data_type",
      "available_data_types",
      "cell_directory",
      "gene_directory",
      "taxonomy",
      "expression_directories"
    )
  )

  for (i in seq_len(nrow(out))) {
    dataset <- out$dataset[[i]]
    spec <- .abc_datasets[[dataset]]

    expect_identical(
      out$modality[[i]],
      as.character(spec$modality)
    )

    expect_identical(
      out$spatial[[i]],
      isTRUE(spec$spatial)
    )

    expect_identical(
      out$imputed[[i]],
      identical(
        spec$modality,
        "imputed_spatial_transcriptomics"
      )
    )

    expect_identical(
      out$default_data_type[[i]],
      as.character(spec$default_data_type[[1L]])
    )

    expect_identical(
      out$available_data_types[[i]],
      paste(
        as.character(spec$available_data_types),
        collapse = ","
      )
    )

    expect_identical(
      out$cell_directory[[i]],
      as.character(spec$cell$dir)
    )

    expect_identical(
      out$gene_directory[[i]],
      as.character(spec$gene$dir)
    )

    expect_identical(
      out$expression_directories[[i]],
      paste(
        as.character(spec$expression_directories),
        collapse = ","
      )
    )

    expected_taxonomy <- if (
      is.null(spec$taxonomy)
    ) {
      NA_character_
    } else {
      as.character(spec$taxonomy[[1L]])
    }

    expect_identical(
      out$taxonomy[[i]],
      expected_taxonomy
    )
  }
})

test_that("list_datasets identifies spatial and imputed routes", {
  out <- list_datasets()

  expect_setequal(
    out$dataset[out$spatial],
    c(
      "HMBA_spatial_Human",
      "HMBA_spatial_Macaque",
      "HMBA_spatial_Marmoset",
      "Zhuang1",
      "Zhuang2",
      "Zhuang3",
      "Zhuang4",
      "MERFISH",
      "MERFISH_imputed"
    )
  )

  expect_identical(
    out$dataset[out$imputed],
    "MERFISH_imputed"
  )

  imputed <- out[
    out$dataset == "MERFISH_imputed",
    ,
    drop = FALSE
  ]

  expect_identical(
    imputed$default_data_type,
    "log2"
  )

  expect_identical(
    imputed$available_data_types,
    "log2"
  )
})

test_that("list_datasets is dependency-light and deterministic", {
  first <- list_datasets()
  second <- list_datasets()

  expect_identical(first, second)
  expect_equal(nrow(first), 23L)
})

test_that("list_datasets returns a printable catalog object", {
  out <- list_datasets()

  expect_s3_class(
    out,
    "abc_dataset_catalog"
  )
  expect_s3_class(
    out,
    "data.frame"
  )

  expect_output(
    print(out),
    "ABC Atlas datasets: 23 routes \\(9 spatial, 1 imputed\\)"
  )

  expect_output(
    print(out),
    "MERFISH_imputed"
  )
})

test_that("catalog remains an ordinary subsettable data frame", {
  out <- list_datasets()

  spatial <- out[out$spatial, , drop = FALSE]

  expect_equal(nrow(spatial), 9L)
  expect_true(all(spatial$spatial))

  expect_setequal(
    spatial$dataset,
    out$dataset[out$spatial]
  )
})
