# Network tests against the public ABC Atlas bucket. Skipped on CRAN and when
# offline; they only touch the release listing and a tiny metadata file.

test_that("releases can be listed from S3", {
  skip_on_cran()
  skip_if_offline("allen-brain-cell-atlas.s3.us-west-2.amazonaws.com")
  m <- .abc_remote_manifests()
  expect_true(length(m) > 0)
  expect_true(all(grepl("^releases/.+\\.json$", m)))
  expect_false(is.unsorted(m))
})

test_that("a small file is downloaded, verified, recorded and reused", {
  skip_on_cran()
  skip_if_offline("allen-brain-cell-atlas.s3.us-west-2.amazonaws.com")
  dir <- withr::local_tempdir()
  cache <- abc_cache_connect(dir)
  expect_false(cache$offline)

  expect_message(path <- cache$get_file_path("ASAP-PMDBS-10X", "value_sets"), "Downloading")
  attrs <- .abc_file_attributes(cache, "ASAP-PMDBS-10X", "value_sets")
  expect_equal(unname(tools::md5sum(path)), attrs$file_hash)
  expect_true(.abc_is_recorded(dir, path))
  expect_false(file.exists(paste0(path, ".part")))

  # Second call reuses the file; offline mode finds it too.
  expect_no_message(cache$get_file_path("ASAP-PMDBS-10X", "value_sets"))
  offline <- abc_cache_connect(dir, offline = TRUE)
  df <- offline$get_metadata_dataframe("ASAP-PMDBS-10X", "value_sets")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})
