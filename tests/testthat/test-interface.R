test_that("the registry has the intended routes without needing Python", {
  x <- list_datasets()
  expect_equal(nrow(x), 23L)
  expect_true(all(c("MERFISH", "MERFISH_imputed", "SEAAD_Multiregion", "SEAAD_CaH") %in% x$dataset))
  expect_false(anyDuplicated(x$dataset) > 0L)
  expect_true(all(c("modality", "available_data_types", "default_data_type",
                    "validation_release", "validation") %in% names(x)))
  expect_identical(x$available_data_types[x$dataset == "MERFISH_imputed"], "log2")
  expect_identical(x$available_data_types[x$dataset == "HMBA_PatchSeq"], "log2")
  expect_identical(x$default_data_type[x$dataset == "HMBA_PatchSeq"], "log2")
  expect_true(all(x$validation_release == "20260711"))
  expect_true(all(list_datasets(release = "latest")$validation ==
                    "DYNAMIC_RELEASE_NOT_VALIDATED"))
  expect_true(all(grepl("^NOT_VALIDATED_FOR_RELEASE_20261031$",
                        list_datasets(release = "20261031")$validation)))
})

test_that("legacy and new wrappers dispatch to the correct dataset", {
  local_mocked_bindings(load_data = function(dataset, ...) dataset)
  expect_identical(load_data_WMB(), "WMB")
  expect_identical(load_data_WHB(), "WHB")
  expect_identical(load_data_AgingMouse(), "AgingMouse")
  expect_identical(load_data_PMDBS(), "PMDBS")
  expect_identical(load_data_HMBA(), "HMBA")
  expect_identical(load_data_MERFISH(), "MERFISH")
  expect_identical(load_data_MERFISH(imputed = TRUE), "MERFISH_imputed")
  for (i in 1:4) expect_identical(load_data_MERFISH(paste0("Zhuang", i)), paste0("Zhuang", i))
  expect_error(load_data_MERFISH("Zhuang1", imputed = TRUE))
  expect_identical(load_data_ConsensusMouse("AIBS"), "ConsensusMouseAIBS")
  expect_identical(load_data_ConsensusMouse("Macosko"), "ConsensusMouseMacosko")
  expect_identical(load_data_DevelopingMouse(), "DevelopingMouse")
  for (s in c("Human", "Macaque", "Marmoset")) {
    expect_identical(load_data_HMBA_species(s), paste0("HMBA_", s))
    expect_identical(load_data_HMBA_spatial(s), paste0("HMBA_spatial_", s))
  }
  expect_identical(load_data_HMBA_PatchSeq(), "HMBA_PatchSeq")
  expect_identical(load_data_SEAAD_Multiregion(), "SEAAD_Multiregion")
  expect_identical(load_data_SEAAD_CaH(), "SEAAD_CaH")
})

test_that("supplemental annotations do not multiply cells", {
  d <- structure(list(cell_metadata = data.frame(cell_label = c("b", "a", "c")), provenance = list()), class = "abc_data")
  a <- data.frame(cell_label = c("a", "b"), region = c("A", "B"))
  x <- join_metadata(d, a)
  expect_identical(x$cell_metadata$annotation_region, c("B", "A", NA_character_))
  expect_identical(x$cell_metadata$cell_label, c("b", "a", "c"))
  expect_error(join_metadata(d, rbind(a, a[1, ])), "duplicated")
})

test_that("no-download and numeric-ID safety checks occur before the backend", {
  expect_error(download_atlas_file("d", "f"), "allow_downloads")
  expect_error(load_data("WMB", cell_ids = 9007199254740993), "character")
})


test_that("supplemental annotations cannot silently cross releases", {
  x <- structure(list(cell_metadata = data.frame(cell_label = "a"),
    release = "releases/20260711/manifest.json", provenance = list()), class = "abc_data")
  y <- data.frame(cell_label = "a", location = "test")
  attr(y, "manifest") <- "releases/20260415/manifest.json"
  expect_error(join_metadata(x, y), "different manifests")
  attr(y, "manifest") <- "20260711"
  expect_equal(join_metadata(x, y)$cell_metadata$annotation_location, "test")
})

test_that("explicit Python argument always wins", {
  explicit_python <- Sys.which("R")
  expect_true(nzchar(explicit_python))
  expect_identical(
    .resolve_abc_python(python = explicit_python),
    path.expand(explicit_python)
  )
})


test_that("managed Python is the default and ignores ambient paths", {
  managed_python <- Sys.which("R")
  ambient_python <- Sys.which("Rscript")
  expect_true(nzchar(managed_python))
  expect_true(nzchar(ambient_python))

  old_mode <- getOption("importABCatlas.python_mode")
  old_env <- Sys.getenv(
    c("IMPORTABCATLAS_PYTHON_MODE", "IMPORTABCATLAS_PYTHON", "RETICULATE_PYTHON"),
    unset = NA_character_
  )
  on.exit({
    options(importABCatlas.python_mode = old_mode)
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
    }
  }, add = TRUE)

  options(importABCatlas.python_mode = "managed")
  Sys.unsetenv("IMPORTABCATLAS_PYTHON_MODE")
  Sys.setenv(IMPORTABCATLAS_PYTHON = ambient_python, RETICULATE_PYTHON = ambient_python)

  local_mocked_bindings(
    .iab_managed_python_status = function(probe = TRUE) {
      list(healthy = TRUE, python = managed_python)
    }
  )

  expect_identical(.resolve_abc_python(), path.expand(managed_python))
})


test_that("advanced external mode uses the deliberate external path", {
  external_python <- Sys.which("Rscript")
  expect_true(nzchar(external_python))

  old_mode <- Sys.getenv("IMPORTABCATLAS_PYTHON_MODE", unset = NA_character_)
  old_python <- Sys.getenv("IMPORTABCATLAS_PYTHON", unset = NA_character_)
  on.exit({
    if (is.na(old_mode)) Sys.unsetenv("IMPORTABCATLAS_PYTHON_MODE") else Sys.setenv(IMPORTABCATLAS_PYTHON_MODE = old_mode)
    if (is.na(old_python)) Sys.unsetenv("IMPORTABCATLAS_PYTHON") else Sys.setenv(IMPORTABCATLAS_PYTHON = old_python)
  }, add = TRUE)

  Sys.setenv(
    IMPORTABCATLAS_PYTHON_MODE = "external",
    IMPORTABCATLAS_PYTHON = external_python
  )
  expect_identical(.resolve_abc_python(), path.expand(external_python))
})


test_that("unsupported expression representations fail before the backend", {
  p <- list(dataset = "HMBA_PatchSeq", default_data_type = "log2",
            available_data_types = "log2")
  expect_error(.resolve_data_type("raw", p), "does not publish")
  expect_error(.resolve_data_type("invalid", p), "data_type must be one of")
  expect_identical(.resolve_data_type(NULL, p), "log2")
})


test_that("managed Python is prepared automatically when missing", {
  managed_python <- Sys.which("R")
  expect_true(nzchar(managed_python))
  calls <- 0L
  old_mode <- getOption("importABCatlas.python_mode")
  old_auto <- getOption("importABCatlas.auto_setup")
  old_env_mode <- Sys.getenv("IMPORTABCATLAS_PYTHON_MODE", unset = NA_character_)
  on.exit({
    options(importABCatlas.python_mode = old_mode,
            importABCatlas.auto_setup = old_auto)
    if (is.na(old_env_mode)) Sys.unsetenv("IMPORTABCATLAS_PYTHON_MODE") else Sys.setenv(IMPORTABCATLAS_PYTHON_MODE = old_env_mode)
  }, add = TRUE)

  options(importABCatlas.python_mode = "managed",
          importABCatlas.auto_setup = TRUE)
  Sys.unsetenv("IMPORTABCATLAS_PYTHON_MODE")

  local_mocked_bindings(
    .iab_managed_python_status = function(probe = TRUE) {
      calls <<- calls + 1L
      if (calls == 1L) list(healthy = FALSE, python = NA_character_) else list(healthy = TRUE, python = managed_python)
    },
    setup_importABCatlas = function(repair = FALSE, quiet = FALSE) invisible(managed_python)
  )

  expect_identical(.resolve_abc_python(), path.expand(managed_python))
  expect_gte(calls, 2L)
})


test_that("automatic setup resets stale external mode", {
  managed_python <- Sys.which("R")
  expect_true(nzchar(managed_python))

  old_mode <- Sys.getenv("IMPORTABCATLAS_PYTHON_MODE", unset = NA_character_)
  old_python <- Sys.getenv("IMPORTABCATLAS_PYTHON", unset = NA_character_)
  old_option_mode <- getOption("importABCatlas.python_mode")
  old_option_python <- getOption("importABCatlas.python")
  on.exit({
    options(importABCatlas.python_mode = old_option_mode,
            importABCatlas.python = old_option_python)
    if (is.na(old_mode)) Sys.unsetenv("IMPORTABCATLAS_PYTHON_MODE") else Sys.setenv(IMPORTABCATLAS_PYTHON_MODE = old_mode)
    if (is.na(old_python)) Sys.unsetenv("IMPORTABCATLAS_PYTHON") else Sys.setenv(IMPORTABCATLAS_PYTHON = old_python)
  }, add = TRUE)

  Sys.setenv(
    IMPORTABCATLAS_PYTHON_MODE = "external",
    IMPORTABCATLAS_PYTHON = "/stale/missing/python"
  )
  options(importABCatlas.python_mode = "external",
          importABCatlas.python = "/stale/missing/python")

  local_mocked_bindings(
    .iab_managed_python_status = function(probe = TRUE) {
      list(healthy = TRUE, python = managed_python)
    }
  )

  expect_silent(setup_importABCatlas(quiet = TRUE))
  expect_identical(Sys.getenv("IMPORTABCATLAS_PYTHON_MODE"), "managed")
  expect_identical(getOption("importABCatlas.python_mode"), "managed")
  expect_identical(getOption("importABCatlas.python"), path.expand(managed_python))
})
