test_that("the registry has the intended routes without needing Python", {
  x <- list_datasets()
  expect_equal(nrow(x), 23L)
  expect_true(all(c("MERFISH", "MERFISH_imputed", "SEAAD_Multiregion", "SEAAD_CaH") %in% x$dataset))
  expect_false(anyDuplicated(x$dataset) > 0L)
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

test_that("Python discovery honors documented precedence", {
  argument_python <- Sys.which("R")
  env_python <- Sys.which("Rscript")
  reticulate_python <- Sys.which("sh")

  expect_true(nzchar(argument_python))
  expect_true(nzchar(env_python))
  expect_true(nzchar(reticulate_python))

  old_option <- getOption("importABCatlas.python")
  old_env <- Sys.getenv(
    c("IMPORTABCATLAS_PYTHON", "RETICULATE_PYTHON", "ABCATLAS_PYTHON"),
    unset = NA_character_
  )

  on.exit({
    options(importABCatlas.python = old_option)
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
      }
    }
  }, add = TRUE)

  options(importABCatlas.python = reticulate_python)
  Sys.setenv(
    IMPORTABCATLAS_PYTHON = env_python,
    RETICULATE_PYTHON = reticulate_python
  )

  expect_identical(.abc_python(argument_python), path.expand(argument_python))
  expect_identical(.abc_python(), path.expand(env_python))
})


test_that("explicit Python argument cannot be overridden", {
  explicit_python <- Sys.which("R")
  env_python <- Sys.which("Rscript")
  reticulate_python <- Sys.which("sh")

  expect_true(nzchar(explicit_python))
  expect_true(nzchar(env_python))
  expect_true(nzchar(reticulate_python))

  old_option <- getOption("importABCatlas.python")
  old_env <- Sys.getenv(
    c("IMPORTABCATLAS_PYTHON", "RETICULATE_PYTHON", "ABCATLAS_PYTHON"),
    unset = NA_character_
  )

  on.exit({
    options(importABCatlas.python = old_option)
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
      }
    }
  }, add = TRUE)

  options(importABCatlas.python = reticulate_python)
  Sys.setenv(
    IMPORTABCATLAS_PYTHON = env_python,
    RETICULATE_PYTHON = reticulate_python,
    ABCATLAS_PYTHON = reticulate_python
  )

  expect_identical(
    .resolve_abc_python(python = explicit_python),
    path.expand(explicit_python)
  )
  expect_identical(
    .abc_python(python = explicit_python),
    path.expand(explicit_python)
  )
})


test_that("documented environment variable wins over R option", {
  env_python <- Sys.which("R")
  option_python <- Sys.which("Rscript")

  expect_true(nzchar(env_python))
  expect_true(nzchar(option_python))

  old_option <- getOption("importABCatlas.python")
  old_env <- Sys.getenv("IMPORTABCATLAS_PYTHON", unset = NA_character_)

  on.exit({
    options(importABCatlas.python = old_option)
    if (is.na(old_env)) {
      Sys.unsetenv("IMPORTABCATLAS_PYTHON")
    } else {
      Sys.setenv(IMPORTABCATLAS_PYTHON = old_env)
    }
  }, add = TRUE)

  options(importABCatlas.python = option_python)
  Sys.setenv(IMPORTABCATLAS_PYTHON = env_python)

  expect_identical(.resolve_abc_python(), path.expand(env_python))
})
