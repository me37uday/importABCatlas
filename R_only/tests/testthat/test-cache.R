test_that("offline connect fails fast without a local manifest", {
  expect_error(abc_cache_connect(withr::local_tempdir(), offline = TRUE),
               "no local manifest")
})

test_that("offline connect loads the local manifest and pins it", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
  expect_s3_class(cache, "abc_cache")
  expect_equal(cache$manifest_name, "releases/20250101/manifest.json")
  expect_equal(cache$version, "20250101")
  # Written without a newline, as abc_atlas_access reads it verbatim.
  expect_identical(readChar(file.path(fx$dir, "_manifest_last_used.txt"), 1000),
                   "releases/20250101/manifest.json")
  expect_output(print(cache), "abc_cache")
})

test_that("file attributes resolve plain names and name/kind", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE)

  a <- .abc_file_attributes(cache, "TEST-10X", "gene")
  expect_equal(a$relative_path, "metadata/TEST-10X/20250101/gene.csv")
  expect_equal(a$file_type, "csv")
  expect_equal(a$local_path, file.path(fx$dir, a$relative_path))

  b <- .abc_file_attributes(cache, "TEST-10X", "TEST-A/raw")
  expect_equal(b$file_type, "h5ad")

  expect_error(.abc_file_attributes(cache, "TEST-10X", "TEST-A"), "TEST-A/raw")
  expect_error(.abc_file_attributes(cache, "TEST-10X", "nope"), "not found")
  expect_error(.abc_file_attributes(cache, "NOPE", "gene"), "not found in manifest")
})

test_that("offline get_file_path errors for files that were never downloaded", {
  fx <- make_fixture_cache()
  cache <- abc_cache_connect(fx$dir, offline = TRUE)
  unlink(file.path(fx$dir, "metadata/TEST-10X/20250101/gene.csv"))
  expect_error(cache$get_file_path("TEST-10X", "gene"), "not available locally")
})

test_that("manifest choice follows the last-used / latest rules", {
  dir <- withr::local_tempdir()
  avail <- c("releases/20240101/manifest.json", "releases/20250101/manifest.json")

  expect_equal(.abc_choose_manifest(dir, avail), avail[2])

  cat(avail[1], file = file.path(dir, "_manifest_last_used.txt"))
  expect_warning(m <- .abc_choose_manifest(dir, avail), "more up to date")
  expect_equal(m, avail[1])

  cat("releases/19990101/manifest.json", file = file.path(dir, "_manifest_last_used.txt"))
  expect_warning(m <- .abc_choose_manifest(dir, avail), "not a valid manifest")
  expect_equal(m, avail[2])

  unlink(file.path(dir, "_manifest_last_used.txt"))
  dir.create(file.path(dir, "releases", "20240101"), recursive = TRUE)
  file.create(file.path(dir, avail[1]))
  expect_warning(m <- .abc_choose_manifest(dir, avail), "newer than")
  expect_equal(m, avail[2])
})

test_that("manifest versions outside the supported window are rejected", {
  expect_silent(.abc_check_manifest_version("20250101"))
  expect_error(.abc_check_manifest_version("20220101"), "supports")
  expect_error(.abc_check_manifest_version("20310101"), "supports")
})

test_that("download record is a JSON array of absolute paths", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "a.csv")
  writeLines("x", f)
  expect_false(.abc_is_recorded(dir, f))
  .abc_record(dir, f)
  .abc_record(dir, f)
  expect_true(.abc_is_recorded(dir, f))
  expect_equal(jsonlite::read_json(file.path(dir, "_downloaded_data.json"),
                                   simplifyVector = TRUE), normalizePath(f))
  .abc_unrecord(dir, f)
  expect_false(.abc_is_recorded(dir, f))
})

test_that("metadata CSVs are read like pandas.read_csv", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("id,txt,num", "007,a ,1", "008,NaN,", "009,None,3"), path)
  df <- .abc_read_csv(path, col_classes = c(id = "character"))
  expect_s3_class(df, "data.frame")
  expect_equal(df$id, c("007", "008", "009"))   # col_classes respected
  expect_equal(df$txt, c("a ", NA, NA))         # whitespace kept, pandas NA strings
  expect_equal(df$num, c(1L, NA, 3L))
})
