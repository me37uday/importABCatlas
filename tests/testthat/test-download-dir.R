test_that("abc_download_dir resolves argument > option > env var > user cache", {
  arg <- withr::local_tempdir()
  opt <- withr::local_tempdir()
  env <- withr::local_tempdir()

  withr::local_options(importABCatlas.download_base = opt)
  withr::local_envvar(ABC_ATLAS_CACHE = env)
  expect_equal(abc_download_dir(arg), normalizePath(arg))
  expect_equal(abc_download_dir(), normalizePath(opt))

  withr::local_options(importABCatlas.download_base = NULL)
  expect_equal(abc_download_dir(), normalizePath(env))

  withr::local_envvar(ABC_ATLAS_CACHE = "")
  expect_equal(abc_download_dir(),
               normalizePath(tools::R_user_dir("importABCatlas", "cache"), mustWork = FALSE))
})

test_that("abc_download_dir creates missing directories", {
  target <- file.path(withr::local_tempdir(), "a", "b")
  expect_false(dir.exists(target))
  abc_download_dir(target)
  expect_true(dir.exists(target))
})
