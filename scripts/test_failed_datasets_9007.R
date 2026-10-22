# Focused live regression for the two routes that exposed feature-axis
# reconciliation defects during the 0.2.0.9004 large-file validation.
#
# Routes:
#   * HMBA_Macaque       -- native macaque gene metadata vs H5AD feature axis
#   * MERFISH_imputed    -- WMB metadata superset vs 8,460 imputed H5AD genes
#
# The script reuses already downloaded H5AD files when available. By default it
# will not start another multi-gigabyte download. Set download_missing=TRUE only
# after reviewing the exact plan and free-space report.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

.iab9007_try <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = force(expr), error = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL,
                             error = conditionMessage(e))
  )
}

.iab9007_elapsed <- function(start) {
  unname(proc.time()[["elapsed"]] - start)
}

.iab9007_human_bytes <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("unknown")
  units <- c("B", "KiB", "MiB", "GiB", "TiB")
  value <- as.numeric(x)
  index <- 1L
  while (value >= 1024 && index < length(units)) {
    value <- value / 1024
    index <- index + 1L
  }
  sprintf("%.2f %s", value, units[[index]])
}

.iab9007_plan_rows <- function(plan, dataset, query) {
  files <- plan$files %||% list()
  if (!length(files)) return(data.frame())
  do.call(rbind, lapply(files, function(file) {
    data.frame(
      dataset = dataset,
      query = query,
      data_type = as.character(plan$data_type %||% NA_character_),
      directory = as.character(file$directory %||% NA_character_),
      file_name = as.character(file$file_name %||% NA_character_),
      feature_matrix_label = as.character(
        file$feature_matrix_label %||% NA_character_
      ),
      n_cells = as.integer(file$n_cells %||% NA_integer_),
      file_bytes = suppressWarnings(as.numeric(file$file_bytes %||% NA_real_)),
      size_source = as.character(file$size_source %||% "unknown"),
      cached = if (is.null(file$cached)) NA else isTRUE(file$cached),
      partial_bytes = suppressWarnings(as.numeric(
        file$partial_bytes %||% NA_real_
      )),
      download_required_bytes = suppressWarnings(as.numeric(
        file$download_required_bytes %||% NA_real_
      )),
      relative_path = as.character(file$relative_path %||% NA_character_),
      size_error = as.character(file$size_error %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }))
}

.iab9007_validate_plan <- function(plan) {
  files <- plan$files %||% list()
  if (!length(files)) return("plan contains no expression files")
  sizes <- vapply(files, function(x) {
    suppressWarnings(as.numeric(x$file_bytes %||% NA_real_))
  }, numeric(1))
  sources <- vapply(files, function(x) {
    as.character(x$size_source %||% "unknown")
  }, character(1))
  paths <- vapply(files, function(x) {
    as.character(x$relative_path %||% "")
  }, character(1))
  if (anyNA(sizes) || any(!is.finite(sizes)) || any(sizes <= 0)) {
    return("one or more selected expression files lack an exact positive size")
  }
  if (any(grepl("directory|folder", sources, ignore.case = TRUE))) {
    return("a directory/folder total was reported as an expression-file size")
  }
  if (any(!nzchar(paths)) || any(!grepl("^expression_matrices/", paths))) {
    return("one or more selected files lack a safe expression_matrices path")
  }
  reported <- suppressWarnings(as.numeric(plan$total_file_bytes %||% NA_real_))
  if (is.na(reported) || !isTRUE(all.equal(reported, sum(sizes), tolerance = 0))) {
    return("total_file_bytes does not equal the exact selected-file sizes")
  }
  NA_character_
}

.iab9007_python_has_anndata <- function(python) {
  if (is.null(python) || !length(python) || !nzchar(python)) return(FALSE)
  status <- suppressWarnings(system2(
    python,
    c("-c", shQuote("import anndata")),
    stdout = FALSE,
    stderr = FALSE
  ))
  identical(as.integer(status), 0L)
}

.iab9007_validate_matrix <- function(result, cells, genes, data_type) {
  if (!is.list(result) || is.null(result$expression)) {
    return("fetch_data did not return the documented matrix-list structure")
  }
  if (!inherits(result$expression, "Matrix")) {
    return("expression is not a sparse Matrix object")
  }
  if (!identical(colnames(result$expression), as.character(cells))) {
    return("returned cell identifiers/order differ from the request")
  }
  if (!identical(rownames(result$expression), as.character(genes))) {
    return("returned gene identifiers/order differ from the request")
  }
  if (!identical(dim(result$expression), c(length(genes), length(cells)))) {
    return("returned dimensions do not match requested genes x cells")
  }
  if (length(result$expression@x) &&
      any(!is.finite(result$expression@x))) {
    return("returned expression contains non-finite values")
  }
  provenance_type <- as.character(result$provenance$data_type %||% "")
  if (!identical(provenance_type, data_type)) {
    return("provenance data_type differs from the requested representation")
  }
  NA_character_
}

.iab9007_overlap_status <- function(simple, complex) {
  common_genes <- intersect(rownames(simple$expression),
                            rownames(complex$expression))
  common_cells <- intersect(colnames(simple$expression),
                            colnames(complex$expression))
  if (!length(common_genes) || !length(common_cells)) {
    return("FAIL_NO_OVERLAP")
  }
  left <- as.matrix(simple$expression[common_genes, common_cells, drop = FALSE])
  right <- as.matrix(complex$expression[common_genes, common_cells, drop = FALSE])
  if (isTRUE(all.equal(left, right, tolerance = 0, check.attributes = TRUE))) {
    "PASS_IDENTICAL"
  } else {
    "FAIL_VALUES_DIFFER"
  }
}

.iab9007_safe_expression_path <- function(cache_root, relative_path) {
  relative_path <- as.character(relative_path[[1L]])
  if (!nzchar(relative_path) ||
      !grepl("^expression_matrices/", relative_path)) {
    stop("Refusing cleanup for an unsafe expression path: ", relative_path,
         call. = FALSE)
  }
  root <- normalizePath(path.expand(cache_root), mustWork = FALSE,
                        winslash = "/")
  candidate <- normalizePath(file.path(root, relative_path), mustWork = FALSE,
                             winslash = "/")
  prefix <- paste0(sub("/+$", "", root), "/")
  if (!startsWith(candidate, prefix)) {
    stop("Refusing cleanup outside cache_root: ", candidate, call. = FALSE)
  }
  candidate
}

.iab9007_checkpoint <- function(summary, files, output_prefix) {
  dir.create(dirname(path.expand(output_prefix)), recursive = TRUE,
             showWarnings = FALSE)
  utils::write.csv(summary, paste0(output_prefix, "_summary.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(files, paste0(output_prefix, "_files.csv"),
                   row.names = FALSE, na = "")
  saveRDS(list(summary = summary, files = files,
               updated_at = as.character(Sys.time())),
          paste0(output_prefix, ".rds"))
}

#' Re-test the two routes that failed 0.2.0.9004 feature-axis validation.
#'
#' @param package_lib Optional R library containing importABCatlas 0.2.0.9007 or later.
#' @param cache_root Cache containing the large H5AD files. The existing 9004
#'   large-fetch cache may be reused. Only exact expression files are eligible
#'   for optional cleanup; metadata are left intact.
#' @param release Pinned Allen release.
#' @param preview_n Number of metadata rows retained by load_data().
#' @param simple_cells,simple_genes Tiny explicit-query dimensions.
#' @param complex_cells,complex_genes Larger same-source query dimensions.
#' @param download_missing Permit an absent H5AD to be downloaded. FALSE by
#'   default so this regression test cannot unexpectedly transfer 7-50 GiB.
#' @param confirm_download Ask for YES before an absent H5AD is downloaded.
#' @param cleanup_expression_after_success Delete the exact H5AD and .part after
#'   its route passes, freeing the large-file space. Failed files are preserved.
#' @param test_object_types Test SCE and Seurat when their packages are present.
#' @param verify_reference Compare the first measured subset to AnnData when the
#'   Python environment contains anndata.
#' @param output_prefix Checkpoint/output prefix; keep it outside cache_root.
#' @return A list with summary, files, and all_selected_routes_pass.
test_failed_datasets_9007 <- function(
    package_lib = NULL,
    cache_root = "~/importABCatlas_9004_large_fetch_cache",
    release = "20260711",
    preview_n = 1000L,
    simple_cells = 3L,
    simple_genes = 2L,
    complex_cells = 8L,
    complex_genes = 8L,
    download_missing = FALSE,
    confirm_download = interactive(),
    cleanup_expression_after_success = TRUE,
    test_object_types = TRUE,
    verify_reference = TRUE,
    output_prefix = "~/Desktop/importABCatlas_9007_failed_routes") {

  if (!is.null(package_lib) && length(package_lib) && nzchar(package_lib[[1L]])) {
    package_lib <- normalizePath(path.expand(package_lib[[1L]]),
                                 mustWork = TRUE, winslash = "/")
    .libPaths(unique(c(package_lib, .libPaths())))
  }
  if ("importABCatlas" %in% loadedNamespaces()) {
    loaded_version <- utils::packageVersion("importABCatlas")
    loaded_path <- normalizePath(
      getNamespaceInfo(asNamespace("importABCatlas"), "path"),
      mustWork = TRUE, winslash = "/"
    )
    if (loaded_version < numeric_version("0.2.0.9007")) {
      stop(
        "An older importABCatlas namespace is already loaded (", loaded_version,
        ") from ", loaded_path,
        ". Restart R before running the 9007 regression.", call. = FALSE
      )
    }
    if (!is.null(package_lib) && length(package_lib) && nzchar(package_lib[[1L]]) &&
        !startsWith(loaded_path, paste0(package_lib, "/"))) {
      stop(
        "importABCatlas is already loaded from a different library: ",
        loaded_path, ". Restart R and source this script before library(importABCatlas).",
        call. = FALSE
      )
    }
  } else {
    loaded <- tryCatch(
      loadNamespace("importABCatlas", lib.loc = package_lib %||% NULL),
      error = function(e) NULL
    )
    if (is.null(loaded)) {
      stop("importABCatlas is not installed in the active R libraries.",
           call. = FALSE)
    }
  }
  installed <- utils::packageVersion("importABCatlas",
                                     lib.loc = package_lib %||% NULL)
  if (installed < numeric_version("0.2.0.9007")) {
    stop("This script requires importABCatlas >= 0.2.0.9007; installed: ",
         installed, ". Restart R after the clean installation.", call. = FALSE)
  }
  importABCatlas::setup_importABCatlas()
  python <- importABCatlas:::.abc_python()
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required.", call. = FALSE)
  }

  cache_root <- normalizePath(path.expand(cache_root), mustWork = FALSE,
                              winslash = "/")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  output_prefix <- path.expand(output_prefix)
  output_dir <- normalizePath(dirname(output_prefix), mustWork = FALSE,
                              winslash = "/")
  cache_prefix <- paste0(sub("/+$", "", cache_root), "/")
  if (identical(output_dir, cache_root) || startsWith(output_dir, cache_prefix)) {
    stop("output_prefix must be outside cache_root.", call. = FALSE)
  }

  old_options <- options(
    importABCatlas.check_disk_space = TRUE,
    importABCatlas.disk_reserve_gb = 1,
    importABCatlas.disk_safety_factor = 1.05,
    importABCatlas.allow_unknown_size = FALSE,
    importABCatlas.network_retries = 8L,
    importABCatlas.network_read_timeout_seconds = 900L
  )
  on.exit(options(old_options), add = TRUE)

  smoke <- importABCatlas::example_data("matrix")
  if (!identical(dim(smoke$expression), c(4L, 3L))) {
    stop("Synthetic R/Python smoke test failed.", call. = FALSE)
  }

  route_specs <- list(
    HMBA_Macaque = list(
      data_type = "raw",
      expected_genes = 35219L,
      previous_problem_genes = c("NCBIGene:114678393", "NCBIGene:114674980")
    ),
    MERFISH_imputed = list(
      data_type = "log2",
      expected_genes = 8460L,
      previous_problem_genes = c("ENSMUSG00000025902", "ENSMUSG00000102343")
    )
  )

  summary <- data.frame()
  files <- data.frame()
  has_anndata <- isTRUE(verify_reference) &&
    .iab9007_python_has_anndata(python)

  for (dataset in names(route_specs)) {
    spec <- route_specs[[dataset]]
    cat("\n============================================================\n")
    cat("9007 focused feature-axis regression: ", dataset, "\n", sep = "")
    cat("============================================================\n")
    row <- data.frame(
      dataset = dataset,
      data_type = spec$data_type,
      route_status = "NOT_RUN",
      load_status = "NOT_RUN",
      load_seconds = NA_real_,
      preview_cells = NA_integer_,
      metadata_gene_rows = NA_integer_,
      plan_status = "NOT_RUN",
      exact_file_bytes = NA_real_,
      cached_at_plan = NA,
      download_required_bytes = NA_real_,
      measured_axis_status = "NOT_RUN",
      measured_axis_seconds = NA_real_,
      measured_genes = NA_integer_,
      expected_measured_genes = spec$expected_genes,
      simple_fetch_status = "NOT_RUN",
      simple_warm_status = "NOT_RUN",
      complex_fetch_status = "NOT_RUN",
      overlap_status = "NOT_RUN",
      previous_problem_status = "NOT_RUN",
      duplicate_identifier_status = if (identical(dataset, "HMBA_Macaque")) "NOT_RUN" else "SKIP_NOT_APPLICABLE",
      independent_reference_status = if (has_anndata) "NOT_RUN" else "SKIP_MISSING_ANNDATA",
      sce_status = "NOT_RUN",
      seurat_status = "NOT_RUN",
      cleanup_status = "NOT_RUN",
      error = NA_character_,
      updated_at = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )

    start <- proc.time()[["elapsed"]]
    loaded <- .iab9007_try(importABCatlas::load_data(
      dataset = dataset,
      download_base = cache_root,
      release = release,
      preview_n = as.integer(preview_n)
    ))
    row$load_seconds <- .iab9007_elapsed(start)
    if (!loaded$ok) {
      row$load_status <- "FAIL"
      row$route_status <- "FAIL_LOAD"
      row$error <- loaded$error
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    data <- loaded$value
    row$load_status <- "PASS"
    row$preview_cells <- nrow(data$cell_metadata)
    row$metadata_gene_rows <- nrow(data$gene_data)

    cell_ids <- head(as.character(data$cell_metadata$cell_label),
                     as.integer(simple_cells))
    plan_result <- .iab9007_try(importABCatlas::plan_fetch(
      data = data,
      cell_ids = cell_ids,
      genes = NULL,
      data_type = spec$data_type
    ))
    if (!plan_result$ok) {
      row$plan_status <- "FAIL"
      row$route_status <- "FAIL_PLAN"
      row$error <- plan_result$error
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    plan <- plan_result$value
    plan_error <- .iab9007_validate_plan(plan)
    if (!is.na(plan_error)) {
      row$plan_status <- "FAIL"
      row$route_status <- "FAIL_PLAN_VALIDATION"
      row$error <- plan_error
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    row$plan_status <- "PASS_EXACT_FILE_SIZES"
    plan_rows <- .iab9007_plan_rows(plan, dataset, "all_measured")
    files <- rbind(files, plan_rows)
    row$exact_file_bytes <- sum(plan_rows$file_bytes)
    row$cached_at_plan <- all(plan_rows$cached %in% TRUE)
    row$download_required_bytes <- sum(plan_rows$download_required_bytes,
                                       na.rm = TRUE)

    if (!isTRUE(row$cached_at_plan) && !isTRUE(download_missing)) {
      row$route_status <- "SKIP_SOURCE_NOT_CACHED"
      row$error <- paste0(
        "Exact source is not fully cached. Required transfer: ",
        .iab9007_human_bytes(row$download_required_bytes),
        ". Re-run with download_missing=TRUE after reviewing plan_fetch()."
      )
      row$cleanup_status <- "PRESERVED"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    if (!isTRUE(row$cached_at_plan) && isTRUE(download_missing) &&
        isTRUE(confirm_download)) {
      answer <- readline(paste0(
        "Download/resume ", .iab9007_human_bytes(row$download_required_bytes),
        " for ", dataset, "? Type YES to continue: "
      ))
      if (!identical(answer, "YES")) {
        row$route_status <- "SKIP_USER_DECLINED_DOWNLOAD"
        row$cleanup_status <- "PRESERVED"
        summary <- rbind(summary, row)
        .iab9007_checkpoint(summary, files, output_prefix)
        next
      }
    }

    start <- proc.time()[["elapsed"]]
    measured_result <- .iab9007_try(importABCatlas::fetch_data(
      data = data,
      cell_ids = cell_ids,
      genes = NULL,
      data_type = spec$data_type,
      return_type = "matrix",
      allow_downloads = isTRUE(download_missing),
      offline = !isTRUE(download_missing),
      verify_reference = has_anndata,
      max_output_mb = 1024
    ))
    row$measured_axis_seconds <- .iab9007_elapsed(start)
    if (!measured_result$ok) {
      row$measured_axis_status <- "FAIL"
      row$route_status <- "FAIL_MEASURED_AXIS"
      row$error <- measured_result$error
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    measured <- measured_result$value
    measured_ids <- as.character(measured$gene_data$gene_identifier)
    row$measured_genes <- length(measured_ids)
    axis_error <- .iab9007_validate_matrix(
      measured, cell_ids, measured_ids, spec$data_type
    )
    if (!is.na(axis_error)) {
      row$measured_axis_status <- "FAIL"
      row$route_status <- "FAIL_MEASURED_AXIS_VALIDATION"
      row$error <- axis_error
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    if (!identical(as.integer(length(measured_ids)),
                   as.integer(spec$expected_genes))) {
      row$measured_axis_status <- "FAIL_COUNT"
      row$route_status <- "FAIL_MEASURED_GENE_COUNT"
      row$error <- paste0(
        "Expected ", spec$expected_genes,
        " measured genes but reconciled ", length(measured_ids), "."
      )
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    row$measured_axis_status <- "PASS"
    if (has_anndata) row$independent_reference_status <- "PASS"

    simple_gene_ids <- head(measured_ids, as.integer(simple_genes))
    complex_gene_ids <- rev(head(measured_ids, as.integer(complex_genes)))
    complex_cells_ids <- rev(head(as.character(data$cell_metadata$cell_label),
                                  as.integer(complex_cells)))

    simple_result <- .iab9007_try(importABCatlas::fetch_data(
      data = data,
      cell_ids = cell_ids,
      genes = simple_gene_ids,
      data_type = spec$data_type,
      return_type = "matrix",
      allow_downloads = FALSE,
      offline = TRUE
    ))
    if (!simple_result$ok) {
      row$simple_fetch_status <- "FAIL"
      row$route_status <- "FAIL_SIMPLE_FETCH"
      row$error <- simple_result$error
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    simple <- simple_result$value
    simple_error <- .iab9007_validate_matrix(
      simple, cell_ids, simple_gene_ids, spec$data_type
    )
    if (!is.na(simple_error)) {
      row$simple_fetch_status <- "FAIL"
      row$route_status <- "FAIL_SIMPLE_AXES"
      row$error <- simple_error
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    row$simple_fetch_status <- "PASS"

    warm_result <- .iab9007_try(importABCatlas::fetch_data(
      data = data,
      cell_ids = cell_ids,
      genes = simple_gene_ids,
      data_type = spec$data_type,
      return_type = "matrix",
      allow_downloads = FALSE,
      offline = TRUE
    ))
    if (!warm_result$ok ||
        !identical(simple$expression, warm_result$value$expression)) {
      row$simple_warm_status <- "FAIL"
      row$route_status <- "FAIL_WARM_REPEAT"
      row$error <- if (!warm_result$ok) warm_result$error else
        "cold/warm expression differs"
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    row$simple_warm_status <- "PASS_IDENTICAL"

    complex_result <- .iab9007_try(importABCatlas::fetch_data(
      data = data,
      cell_ids = complex_cells_ids,
      genes = complex_gene_ids,
      data_type = spec$data_type,
      return_type = "matrix",
      allow_downloads = FALSE,
      offline = TRUE
    ))
    if (!complex_result$ok) {
      row$complex_fetch_status <- "FAIL"
      row$route_status <- "FAIL_COMPLEX_FETCH"
      row$error <- complex_result$error
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    complex <- complex_result$value
    complex_error <- .iab9007_validate_matrix(
      complex, complex_cells_ids, complex_gene_ids, spec$data_type
    )
    if (!is.na(complex_error)) {
      row$complex_fetch_status <- "FAIL"
      row$route_status <- "FAIL_COMPLEX_AXES"
      row$error <- complex_error
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }
    row$complex_fetch_status <- "PASS"
    row$overlap_status <- .iab9007_overlap_status(simple, complex)
    if (!identical(row$overlap_status, "PASS_IDENTICAL")) {
      row$route_status <- "FAIL_OVERLAP"
      row$error <- row$overlap_status
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }

    problem <- intersect(spec$previous_problem_genes,
                         as.character(data$gene_data$gene_identifier))
    if (!length(problem)) {
      row$previous_problem_status <- "PASS_NOT_IN_METADATA"
    } else {
      problem_fetch <- .iab9007_try(importABCatlas::fetch_data(
        data = data,
        cell_ids = head(cell_ids, 1L),
        genes = problem,
        data_type = spec$data_type,
        return_type = "matrix",
        allow_downloads = FALSE,
        offline = TRUE
      ))
      if (all(problem %in% measured_ids)) {
        row$previous_problem_status <- if (
          problem_fetch$ok &&
          identical(rownames(problem_fetch$value$expression), problem)
        ) "PASS_RESOLVED_BY_ALIAS_RECONCILIATION" else "FAIL_EXPECTED_MEASURED"
      } else {
        row$previous_problem_status <- if (
          !problem_fetch$ok &&
          grepl("unmeasured/absent", problem_fetch$error, ignore.case = TRUE)
        ) "PASS_STRICT_UNMEASURED_ERROR" else "FAIL_EXPECTED_STRICT_ERROR"
      }
    }
    if (startsWith(row$previous_problem_status, "FAIL")) {
      row$route_status <- "FAIL_PREVIOUS_PROBLEM_CHECK"
      row$error <- row$previous_problem_status
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }

    # Native HMBA Macaque contains duplicated source NCBI identifiers.  The
    # package must expose deterministic unique identifiers and be able to fetch
    # them explicitly without collapsing rows or falling back to an ambiguous
    # gene symbol.
    if (identical(dataset, "HMBA_Macaque")) {
      duplicate_ids <- measured_ids[
        grepl("\\.iabdup[1-9][0-9]*$", measured_ids)
      ]
      if (!length(duplicate_ids)) {
        row$duplicate_identifier_status <- "FAIL_NO_DISAMBIGUATED_IDS"
      } else {
        duplicate_ids <- head(duplicate_ids, min(4L, length(duplicate_ids)))
        duplicate_fetch <- .iab9007_try(importABCatlas::fetch_data(
          data = data,
          cell_ids = head(cell_ids, 2L),
          genes = duplicate_ids,
          data_type = spec$data_type,
          return_type = "matrix",
          allow_downloads = FALSE,
          offline = TRUE
        ))
        duplicate_error <- if (duplicate_fetch$ok) {
          .iab9007_validate_matrix(
            duplicate_fetch$value,
            head(cell_ids, 2L),
            duplicate_ids,
            spec$data_type
          )
        } else {
          duplicate_fetch$error
        }
        row$duplicate_identifier_status <- if (is.na(duplicate_error)) {
          "PASS_EXPLICIT_DISAMBIGUATED_FETCH"
        } else {
          paste0("FAIL: ", duplicate_error)
        }
      }
      if (startsWith(row$duplicate_identifier_status, "FAIL")) {
        row$route_status <- "FAIL_DUPLICATE_IDENTIFIER_FETCH"
        row$error <- row$duplicate_identifier_status
        row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
        summary <- rbind(summary, row)
        .iab9007_checkpoint(summary, files, output_prefix)
        next
      }
    }

    if (isTRUE(test_object_types) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      sce <- .iab9007_try(importABCatlas::fetch_data(
        data = data, cell_ids = cell_ids, genes = simple_gene_ids,
        data_type = spec$data_type, return_type = "sce",
        allow_downloads = FALSE, offline = TRUE
      ))
      row$sce_status <- if (sce$ok) "PASS" else paste0("FAIL: ", sce$error)
    } else {
      row$sce_status <- "SKIP_DEPENDENCY_OR_DISABLED"
    }
    if (isTRUE(test_object_types) &&
        requireNamespace("SeuratObject", quietly = TRUE) &&
        !any(grepl("[_|]", simple_gene_ids))) {
      seurat <- .iab9007_try(importABCatlas::fetch_data(
        data = data, cell_ids = cell_ids, genes = simple_gene_ids,
        data_type = spec$data_type, return_type = "seurat",
        allow_downloads = FALSE, offline = TRUE
      ))
      row$seurat_status <- if (seurat$ok) "PASS" else paste0("FAIL: ", seurat$error)
    } else {
      row$seurat_status <- "SKIP_DEPENDENCY_ID_OR_DISABLED"
    }
    object_failures <- c(row$sce_status, row$seurat_status)
    object_failures <- object_failures[startsWith(object_failures, "FAIL")]
    if (length(object_failures)) {
      row$route_status <- "FAIL_OBJECT_CONSTRUCTION"
      row$error <- paste(object_failures, collapse = "; ")
      row$cleanup_status <- "PRESERVED_AFTER_FAILURE"
      summary <- rbind(summary, row)
      .iab9007_checkpoint(summary, files, output_prefix)
      next
    }

    row$route_status <- "PASS"
    row$cleanup_status <- "PRESERVED_BY_REQUEST"
    if (isTRUE(cleanup_expression_after_success)) {
      unique_paths <- unique(plan_rows$relative_path)
      cleanup_errors <- character()
      for (relative in unique_paths) {
        local <- .iab9007_safe_expression_path(cache_root, relative)
        for (candidate in c(local, paste0(local, ".part"))) {
          if (file.exists(candidate) && !isTRUE(unlink(candidate, force = TRUE) == 0L)) {
            cleanup_errors <- c(cleanup_errors, candidate)
          }
        }
      }
      row$cleanup_status <- if (length(cleanup_errors)) {
        paste0("FAIL: ", paste(cleanup_errors, collapse = "; "))
      } else {
        "PASS_DELETED_EXPRESSION_FILES"
      }
    }

    row$updated_at <- as.character(Sys.time())
    summary <- rbind(summary, row)
    .iab9007_checkpoint(summary, files, output_prefix)
    rm(data, measured, simple, complex)
    invisible(gc())
  }

  all_pass <- nrow(summary) == length(route_specs) &&
    all(summary$route_status == "PASS")
  result <- list(
    summary = summary,
    files = files,
    package_version = as.character(installed),
    cache_root = cache_root,
    release = release,
    all_selected_routes_pass = all_pass
  )
  saveRDS(result, paste0(output_prefix, ".rds"))
  cat("\nFocused 9007 result: ", all_pass, "\n", sep = "")
  print(summary)
  invisible(result)
}
