# Comprehensive live validation for importABCatlas.
#
# This script exercises the public catalogue, diagnostics, synthetic example,
# load_data(), get_metadata(), download_atlas_file(), join_metadata(),
# list_files(), plan_fetch(), fetch_data(), cache reuse, object construction,
# and expected error paths. Each dataset uses its own disposable cache. The
# complete per-dataset cache (metadata, manifest, and expression files) is
# deleted after that dataset, even when a test fails.
#
# IMPORTANT:
# * Allen expression assets are complete H5AD files. A tiny returned matrix can
#   still require a large source-file download.
# * The safe default fetches only when the exact additional download is <= 2 GiB.
# * Use max_fetch_gb = Inf only when you deliberately want to download every
#   planned source file one dataset at a time.
# * No pre-existing user cache is touched. A dedicated test root with a marker
#   file is mandatory.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.ia_sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

.ia_elapsed <- function(start) unname(proc.time()[["elapsed"]] - start)

.ia_error <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = force(expr), error = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
}

.ia_bytes <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(NA_real_)
  as.numeric(x[[1L]])
}

.ia_human_bytes <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NA_character_)
  units <- c("B", "KiB", "MiB", "GiB", "TiB")
  value <- as.numeric(x); i <- 1L
  while (is.finite(value) && value >= 1024 && i < length(units)) {
    value <- value / 1024; i <- i + 1L
  }
  sprintf("%.2f %s", value, units[[i]])
}

.ia_dir_info <- function(path) {
  if (!dir.exists(path)) return(list(files = 0L, bytes = 0))
  fs <- list.files(path, recursive = TRUE, full.names = TRUE,
                   all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  info <- if (length(fs)) file.info(fs) else NULL
  list(files = length(fs), bytes = if (is.null(info)) 0 else sum(info$size, na.rm = TRUE))
}

.ia_expression_cache_files <- function(path) {
  if (!dir.exists(path)) return(character())
  fs <- list.files(path, recursive = TRUE, full.names = TRUE,
                   all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  normalized <- gsub("\\\\", "/", fs)
  sort(fs[grepl("/expression_matrices/|\\.h5ad($|\\.)", normalized, ignore.case = TRUE)])
}

.ia_plan_rows <- function(plan, dataset, query_name) {
  fs <- plan$files %||% list()
  if (!length(fs)) return(data.frame())
  do.call(rbind, lapply(fs, function(f) data.frame(
    dataset = dataset,
    query = query_name,
    data_type = as.character(plan$data_type %||% NA_character_),
    directory = as.character(f$directory %||% NA_character_),
    file_name = as.character(f$file_name %||% NA_character_),
    feature_matrix_label = as.character(f$feature_matrix_label %||% NA_character_),
    n_cells = as.integer(f$n_cells %||% NA_integer_),
    file_bytes = .ia_bytes(f$file_bytes),
    size_source = as.character(f$size_source %||% "unknown"),
    cached = if (is.null(f$cached)) NA else isTRUE(f$cached),
    partial_bytes = .ia_bytes(f$partial_bytes),
    download_required_bytes = .ia_bytes(f$download_required_bytes),
    relative_path = as.character(f$relative_path %||% NA_character_),
    size_error = as.character(f$size_error %||% NA_character_),
    stringsAsFactors = FALSE
  )))
}

.ia_validate_plan <- function(plan) {
  fs <- plan$files %||% list()
  if (!length(fs)) return(list(ok = FALSE, error = "plan contains no expression files"))
  keys <- vapply(fs, function(f) paste(f$directory, f$file_name, sep = "/"), character(1))
  unique_fs <- fs[!duplicated(keys)]
  sizes <- vapply(unique_fs, function(f) .ia_bytes(f$file_bytes), numeric(1))
  sources <- vapply(unique_fs, function(f) as.character(f$size_source %||% "unknown"), character(1))
  if (anyNA(sizes) || any(!is.finite(sizes)) || any(sizes < 0)) {
    return(list(ok = FALSE, error = "one or more selected expression files lack an exact non-negative file size"))
  }
  if (any(grepl("directory|folder", sources, ignore.case = TRUE))) {
    return(list(ok = FALSE, error = "plan reported a directory/folder size instead of an exact selected-file size"))
  }
  total <- sum(sizes)
  reported <- .ia_bytes(plan$total_file_bytes)
  if (is.na(reported) || !isTRUE(all.equal(total, reported, tolerance = 0))) {
    return(list(ok = FALSE, error = sprintf("total_file_bytes mismatch: calculated %.0f, reported %s",
                                           total, as.character(reported))))
  }
  list(ok = TRUE, error = NA_character_)
}

.ia_simple_cells <- function(cm, n) {
  n <- min(as.integer(n), nrow(cm))
  if (n < 1L) return(character())
  if ("feature_matrix_label" %in% names(cm)) {
    labels <- as.character(cm$feature_matrix_label)
    good <- !is.na(labels) & nzchar(labels)
    if (any(good)) {
      counts <- sort(table(labels[good]), decreasing = TRUE)
      label <- names(counts)[[1L]]
      ids <- cm$cell_label[good & labels == label]
      return(head(as.character(ids), n))
    }
  }
  head(as.character(cm$cell_label), n)
}

.ia_complex_cells <- function(cm, n) {
  n <- min(as.integer(n), nrow(cm))
  if (n < 1L) return(list(ids = character(), mode = "none"))
  if ("feature_matrix_label" %in% names(cm)) {
    labels <- as.character(cm$feature_matrix_label)
    groups <- unique(labels[!is.na(labels) & nzchar(labels)])
    if (length(groups) >= 2L) {
      groups <- groups[seq_len(min(2L, length(groups)))]
      per_group <- max(1L, ceiling(n / length(groups)))
      selections <- lapply(groups, function(g)
        head(as.character(cm$cell_label[labels == g]), per_group))
      # Interleave groups so the requested order is deliberately different from
      # the file-reading order.
      ids <- as.vector(t(do.call(cbind, lapply(selections, function(x) {
        length(x) <- max(lengths(selections)); x
      }))))
      ids <- unique(ids[!is.na(ids) & nzchar(ids)])
      return(list(ids = head(ids, n), mode = "multi_source_if_manifest_routes_labels"))
    }
  }
  ids <- head(as.character(cm$cell_label), n)
  list(ids = rev(ids), mode = "single_source_larger_mixed_query")
}

.ia_gene_query <- function(gd, n, mix_symbols = FALSE) {
  n <- min(as.integer(n), nrow(gd))
  idx <- seq_len(n)
  expected <- as.character(gd$gene_identifier[idx])
  query <- expected
  if (isTRUE(mix_symbols) && "gene_symbol" %in% names(gd)) {
    symbols <- as.character(gd$gene_symbol)
    counts <- table(symbols[!is.na(symbols) & nzchar(symbols)])
    for (j in idx[seq_along(idx) %% 2L == 0L]) {
      symbol <- symbols[[j]]
      count <- if (!is.na(symbol) && nzchar(symbol)) unname(counts[symbol]) else NA_integer_
      if (!is.na(symbol) && nzchar(symbol) && length(count) == 1L && !is.na(count) && count == 1L)
        query[[j]] <- symbol
    }
  }
  list(query = query, expected = expected)
}

.ia_first_metadata_file <- function(data) {
  sources <- data$provenance$metadata_files %||% list()
  if (!length(sources)) return(NULL)
  item <- sources[[1L]]
  if (is.list(item) && length(item) >= 2L)
    return(list(directory = as.character(item[[1L]]), file_name = as.character(item[[2L]])))
  NULL
}

.ia_validate_matrix_result <- function(x, expected_cells, expected_genes = NULL) {
  if (!is.list(x) || is.null(x$expression)) return("fetch_data did not return a matrix list")
  if (!inherits(x$expression, "Matrix")) return("expression is not a Matrix sparse object")
  if (!identical(colnames(x$expression), as.character(expected_cells)))
    return("returned cell axis/order differs from the request")
  observed_genes <- rownames(x$expression)
  if (is.null(observed_genes) || anyNA(observed_genes) || any(!nzchar(observed_genes)) ||
      anyDuplicated(observed_genes))
    return("returned measured-gene axis is missing, blank, or duplicated")
  if (!is.null(expected_genes) &&
      !identical(observed_genes, as.character(expected_genes)))
    return("returned gene axis/order differs from resolved gene identifiers")
  expected_n_genes <- if (is.null(expected_genes)) length(observed_genes) else length(expected_genes)
  if (nrow(x$expression) != expected_n_genes || ncol(x$expression) != length(expected_cells))
    return("returned matrix dimensions are incorrect")
  if (any(!is.finite(x$expression@x))) return("returned sparse values contain non-finite values")
  NA_character_
}

.ia_checkpoint <- function(validation, output_prefix) {
  if (!is.null(validation$summary) && nrow(validation$summary))
    utils::write.csv(validation$summary, paste0(output_prefix, "_summary.csv"), row.names = FALSE, na = "")
  if (!is.null(validation$files) && nrow(validation$files))
    utils::write.csv(validation$files, paste0(output_prefix, "_files.csv"), row.names = FALSE, na = "")
  saveRDS(validation, paste0(output_prefix, ".rds"))
}

#' Run comprehensive live validation across every bundled dataset route.
#'
#' @param test_root Dedicated disposable cache root. It must be empty or already
#'   contain this script's marker file. Existing user caches are never used.
#' @param release Pinned Allen release.
#' @param datasets NULL for every bundled route, or a character vector of route IDs to test.
#' @param preview_n Metadata rows retained in each abc_data object. The upstream
#'   metadata file may still be downloaded in full.
#' @param max_fetch_gb Maximum exact additional expression download permitted per
#'   simple or complex fetch. Use Inf for exhaustive downloads.
#' @param run_fetch Whether to execute expression fetches after planning.
#' @param test_alternate_plan Try the alternate raw/log2 plan without downloading.
#' @param fetch_both_data_types Also fetch the alternate representation when it
#'   exists and is below max_fetch_gb. Disabled by default to conserve bandwidth.
#' @param test_object_types Test cached SCE/Seurat/spatial materialization where
#'   dependencies and modality permit.
#' @param cleanup Delete the complete per-dataset cache after each iteration.
#' @param network_retries Number of package-level retries for transient Allen/S3 failures.
#' @param network_read_timeout_seconds Per-request read timeout used by the resumable downloader.
#' @return A list with summary, exact file plans, global checks, and a named
#'   metadata_runtime vector.
test_all_functions <- function(
    test_root = "~/importABCatlas_full_validation_cache",
    release = "20260711",
    datasets = NULL,
    preview_n = 1000L,
    simple_cells = 3L,
    simple_genes = 2L,
    complex_cells = 8L,
    complex_genes = 6L,
    run_fetch = TRUE,
    max_fetch_gb = 2,
    test_alternate_plan = TRUE,
    fetch_both_data_types = FALSE,
    test_object_types = TRUE,
    cleanup = TRUE,
    network_retries = 5L,
    network_read_timeout_seconds = 600L,
    output_prefix = "importABCatlas_full_validation") {

  if (!requireNamespace("importABCatlas", quietly = TRUE))
    stop("Install importABCatlas before running this script.", call. = FALSE)
  required_version <- numeric_version("0.2.0.9007")
  installed_version <- utils::packageVersion("importABCatlas")
  if (installed_version < required_version)
    stop(sprintf("This validation script requires importABCatlas >= %s; installed: %s. Restart R after installing the new candidate.",
                 required_version, installed_version), call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Matrix is required for validation.", call. = FALSE)
  for (value in list(preview_n, simple_cells, simple_genes, complex_cells, complex_genes))
    if (length(value) != 1L || is.na(value) || value < 1L)
      stop("preview/cell/gene counts must be positive scalars.", call. = FALSE)
  if (length(max_fetch_gb) != 1L || is.na(max_fetch_gb) || max_fetch_gb < 0)
    stop("max_fetch_gb must be non-negative or Inf.", call. = FALSE)

  root <- normalizePath(path.expand(test_root), mustWork = FALSE, winslash = "/")
  marker <- file.path(root, ".importABCatlas_disposable_validation_root")
  if (dir.exists(root)) {
    contents <- list.files(root, all.files = TRUE, no.. = TRUE)
    if (length(contents) && !file.exists(marker))
      stop("test_root is not empty and lacks the validation marker. Choose a new disposable directory; no files were deleted.", call. = FALSE)
  }
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(marker)) writeLines("Created by scripts/test_all_functions.R; contents are disposable.", marker)

  max_fetch_bytes <- if (is.infinite(max_fetch_gb)) Inf else as.numeric(max_fetch_gb) * 1024^3
  old_options <- options(
    importABCatlas.check_disk_space = TRUE,
    importABCatlas.disk_reserve_gb = 1,
    importABCatlas.disk_safety_factor = 1.05,
    importABCatlas.allow_unknown_size = FALSE,
    importABCatlas.network_retries = as.integer(network_retries),
    importABCatlas.network_backoff_seconds = 2,
    importABCatlas.network_connect_timeout_seconds = 30L,
    importABCatlas.network_read_timeout_seconds = as.integer(network_read_timeout_seconds),
    importABCatlas.download_chunk_mb = 8
  )
  on.exit(options(old_options), add = TRUE)

  full_catalogue <- importABCatlas::list_datasets()
  catalogue <- full_catalogue
  if (!is.null(datasets)) {
    datasets <- as.character(datasets)
    if (!length(datasets) || anyNA(datasets) || any(!nzchar(datasets)) || anyDuplicated(datasets))
      stop("datasets must be NULL or a nonempty character vector of unique route IDs.", call. = FALSE)
    missing_routes <- setdiff(datasets, full_catalogue$dataset)
    if (length(missing_routes))
      stop("Unknown dataset route(s): ", paste(missing_routes, collapse = ", "), call. = FALSE)
    catalogue <- full_catalogue[match(datasets, full_catalogue$dataset), , drop = FALSE]
    rownames(catalogue) <- NULL
  }
  global_rows <- list()
  add_global <- function(name, result, detail = NA_character_) {
    global_rows[[length(global_rows) + 1L]] <<- data.frame(
      check = name, status = if (result$ok) "PASS" else "FAIL",
      detail = if (result$ok) detail else result$error, stringsAsFactors = FALSE)
  }
  setup_result <- .ia_error(importABCatlas::setup_importABCatlas())
  add_global("setup_importABCatlas", setup_result,
             if (setup_result$ok) "managed Python ready" else NA_character_)
  if (!setup_result$ok) {
    stop("Automatic Python setup failed: ", setup_result$error, call. = FALSE)
  }
  diag_result <- .ia_error(importABCatlas::diagnose())
  add_global("diagnose", diag_result,
             if (diag_result$ok) as.character(diag_result$value$python$executable %||% NA_character_) else NA_character_)
  add_global("list_datasets", list(ok = nrow(full_catalogue) > 0L, error = "empty catalogue"),
             sprintf("%d bundled routes; %d selected for this run", nrow(full_catalogue), nrow(catalogue)))
  demo <- .ia_error(importABCatlas::example_data("matrix"))
  add_global("example_data_matrix", demo,
             if (demo$ok) paste(dim(demo$value$expression), collapse = " x ") else NA_character_)
  if (isTRUE(test_object_types)) {
    for (type in c("sce", "seurat", "spatial")) {
      needed <- switch(type, sce = "SingleCellExperiment", seurat = "SeuratObject", spatial = "SpatialExperiment")
      if (!requireNamespace(needed, quietly = TRUE)) {
        global_rows[[length(global_rows) + 1L]] <- data.frame(
          check = paste0("example_data_", type), status = "SKIP_DEPENDENCY",
          detail = paste("missing", needed), stringsAsFactors = FALSE)
      } else {
        result <- .ia_error(importABCatlas::example_data(type))
        add_global(paste0("example_data_", type), result,
                   if (result$ok) paste(class(result$value), collapse = "/") else NA_character_)
      }
    }
  }
  global <- do.call(rbind, global_rows)

  summary_rows <- list(); file_rows <- list()
  checkpoint <- function() {
    validation <- list(
      summary = if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame(),
      files = if (length(file_rows)) do.call(rbind, file_rows) else data.frame(),
      global = global,
      metadata_runtime = if (length(summary_rows)) {
        x <- do.call(rbind, summary_rows)
        stats::setNames(x$load_seconds, x$dataset)
      } else numeric()
    )
    .ia_checkpoint(validation, output_prefix)
    validation
  }

  for (i in seq_len(nrow(catalogue))) {
    dataset <- as.character(catalogue$dataset[[i]])
    is_spatial <- isTRUE(catalogue$spatial[[i]])
    primary_type <- if ("default_data_type" %in% names(catalogue)) {
      as.character(catalogue$default_data_type[[i]])
    } else if (isTRUE(catalogue$imputed[[i]]) || identical(dataset, "MERFISH_imputed")) {
      "log2"
    } else {
      "raw"
    }
    supported_types <- if ("available_data_types" %in% names(catalogue)) {
      strsplit(as.character(catalogue$available_data_types[[i]]), ",", fixed = TRUE)[[1L]]
    } else c("raw", "log2")
    supported_types <- supported_types[nzchar(supported_types)]
    alternate_type <- if (identical(primary_type, "raw")) "log2" else "raw"
    alternate_is_published <- alternate_type %in% supported_types
    cache_dir <- file.path(root, sprintf("%02d_%s", i, .ia_sanitize(dataset)))
    if (dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE, force = TRUE)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    cat(sprintf("\n[%d/%d] %s | primary=%s | isolated cache=%s\n",
                i, nrow(catalogue), dataset, primary_type, cache_dir))

    row <- as.list(stats::setNames(rep(NA, 55), c(
      "dataset", "primary_data_type", "load_status", "load_seconds", "preview_cells", "candidate_genes",
      "list_files_status", "get_metadata_status", "download_file_status", "join_metadata_status",
      "simple_plan_status", "simple_plan_seconds", "simple_plan_size_status", "simple_files",
      "simple_total_file_bytes", "simple_download_bytes", "simple_disk_enough",
      "alternate_plan_status", "alternate_data_type", "alternate_fetch_status", "alternate_fetch_seconds", "alternate_fetch_error",
      "error_no_consent", "error_insufficient_disk", "error_tiny_output", "error_missing_gene",
      "error_duplicate_cells", "error_invalid_data_type", "no_expression_download_before_fetch",
      "simple_fetch_status", "simple_fetch_seconds", "simple_warm_status", "simple_warm_seconds",
      "simple_axes_status", "complex_mode", "complex_plan_status", "complex_plan_seconds",
      "complex_fetch_status", "complex_fetch_seconds", "complex_axes_status", "complex_sources",
      "sce_status", "seurat_status", "spatial_status", "fetch_skip_reason", "cache_files_deleted",
      "cache_bytes_deleted", "cleanup_status", "load_error", "simple_plan_error", "alternate_plan_error",
      "simple_fetch_error", "complex_plan_error", "complex_fetch_error", "notes"
    )))
    row$dataset <- dataset; row$primary_data_type <- primary_type; row$alternate_data_type <- alternate_type
    row$load_status <- "NOT_RUN"; row$simple_plan_status <- "NOT_RUN"; row$alternate_plan_status <- "NOT_RUN"
    row$simple_fetch_status <- "NOT_RUN"; row$complex_fetch_status <- "NOT_RUN"
    row$cleanup_status <- "PENDING"; row$notes <- ""

    dat <- simple_result <- warm_result <- complex_result <- alternate_result <- NULL
    simple_plan <- complex_plan <- alternate_plan <- NULL

    run <- tryCatch({
      start <- proc.time()[["elapsed"]]
      loaded <- .ia_error(importABCatlas::load_data(
        dataset = dataset, download_base = cache_dir, release = release,
        preview_n = as.integer(preview_n)))
      row$load_seconds <- .ia_elapsed(start)
      if (!loaded$ok) {
        row$load_status <- "FAIL"; row$load_error <- loaded$error
        stop(structure(list(message = loaded$error), class = c("ia_load_stop", "error", "condition")))
      }
      dat <- loaded$value
      row$load_status <- "PASS"; row$preview_cells <- nrow(dat$cell_metadata); row$candidate_genes <- nrow(dat$gene_data)
      if (!all(c("cell_label") %in% names(dat$cell_metadata)) || !all(c("gene_identifier") %in% names(dat$gene_data)))
        stop("load_data returned no exact cell_label/gene_identifier axes")

      # Exercise catalogue access for the selected directory only.
      directory <- as.character(dat$provenance$spec$directory)
      lf <- .ia_error(importABCatlas::list_files(directory = directory, download_base = cache_dir,
                                                 release = release, offline = TRUE))
      row$list_files_status <- if (lf$ok) "PASS" else paste0("FAIL: ", lf$error)

      # Exercise direct metadata retrieval and cached asset-path access.
      source <- .ia_first_metadata_file(dat)
      if (is.null(source)) {
        row$get_metadata_status <- row$download_file_status <- "SKIP_NO_SOURCE"
      } else {
        gm <- .ia_error(importABCatlas::get_metadata(
          source$directory, source$file_name, download_base = cache_dir,
          release = release, offline = TRUE, preview_n = 2L))
        row$get_metadata_status <- if (gm$ok && is.data.frame(gm$value)) "PASS" else paste0("FAIL: ", gm$error)
        df <- .ia_error(importABCatlas::download_atlas_file(
          source$directory, source$file_name, download_base = cache_dir,
          release = release, offline = TRUE))
        row$download_file_status <- if (df$ok && file.exists(as.character(df$value))) "PASS" else paste0("FAIL: ", df$error)
      }

      # Exercise validated supplemental cell join without touching source data.
      annotation <- data.frame(
        cell_label = head(as.character(dat$cell_metadata$cell_label), min(2L, nrow(dat$cell_metadata))),
        validation_flag = "joined", stringsAsFactors = FALSE)
      attr(annotation, "manifest") <- dat$release
      joined <- .ia_error(importABCatlas::join_metadata(dat, annotation, prefix = "test_"))
      row$join_metadata_status <- if (joined$ok && nrow(joined$value$cell_metadata) == nrow(dat$cell_metadata)) "PASS" else paste0("FAIL: ", joined$error)

      simple_ids <- .ia_simple_cells(dat$cell_metadata, simple_cells)
      authoritative_axis_route <- dataset %in% c("HMBA_Macaque", "MERFISH_imputed")
      simple_gene <- if (authoritative_axis_route) {
        list(query = NULL, expected = NULL)
      } else {
        .ia_gene_query(dat$gene_data, simple_genes, mix_symbols = FALSE)
      }
      start <- proc.time()[["elapsed"]]
      planned <- .ia_error(importABCatlas::plan_fetch(
        data = dat, cell_ids = simple_ids, genes = simple_gene$query, data_type = primary_type))
      row$simple_plan_seconds <- .ia_elapsed(start)
      if (!planned$ok) {
        row$simple_plan_status <- "FAIL"; row$simple_plan_error <- planned$error
        stop(structure(list(message = planned$error), class = c("ia_plan_stop", "error", "condition")))
      }
      simple_plan <- planned$value; row$simple_plan_status <- "PASS"
      check <- .ia_validate_plan(simple_plan)
      row$simple_plan_size_status <- if (check$ok) "PASS_EXACT_FILE_SIZES" else paste0("FAIL: ", check$error)
      row$simple_files <- paste(unique(vapply(simple_plan$files, function(f)
        paste(f$directory, f$file_name, sep = "/"), character(1))), collapse = " | ")
      row$simple_total_file_bytes <- .ia_bytes(simple_plan$total_file_bytes)
      row$simple_download_bytes <- .ia_bytes(simple_plan$total_download_bytes)
      row$simple_disk_enough <- if (is.null(simple_plan$disk$enough_space)) NA else isTRUE(simple_plan$disk$enough_space)
      file_rows[[length(file_rows) + 1L]] <- .ia_plan_rows(simple_plan, dataset, "simple_primary")

      if (isTRUE(test_alternate_plan)) {
        alt <- .ia_error(importABCatlas::plan_fetch(
          data = dat, cell_ids = simple_ids, genes = simple_gene$query, data_type = alternate_type))
        if (isTRUE(alternate_is_published)) {
          row$alternate_plan_status <- if (alt$ok) "PASS" else "FAIL"
          row$alternate_plan_error <- alt$error
          if (alt$ok) {
            alternate_plan <- alt$value
            file_rows[[length(file_rows) + 1L]] <- .ia_plan_rows(alternate_plan, dataset, "simple_alternate")
          }
        } else {
          row$alternate_plan_status <- if (!alt$ok && grepl(
            "does not publish|Available data_type values", alt$error, ignore.case = TRUE
          )) "PASS_UNAVAILABLE_DECLARED" else paste0("FAIL_UNPUBLISHED_ACCEPTED: ", alt$error)
          row$alternate_plan_error <- alt$error
        }
      } else row$alternate_plan_status <- "SKIP_DISABLED"

      # Expected error handling. These calls must not download expression.
      before_expression <- .ia_expression_cache_files(cache_dir)
      consent <- .ia_error(importABCatlas::fetch_data(
        data = dat, cell_ids = simple_ids, genes = simple_gene$query,
        data_type = primary_type, return_type = "matrix", allow_downloads = FALSE))
      row$error_no_consent <- if (!consent$ok && grepl("consent|allow_downloads", consent$error, ignore.case = TRUE)) "PASS" else paste0("FAIL: ", consent$error)

      tiny <- .ia_error(importABCatlas::fetch_data(
        data = dat, cell_ids = simple_ids, genes = simple_gene$query,
        data_type = primary_type, return_type = "matrix", allow_downloads = TRUE,
        max_output_mb = .Machine$double.eps))
      row$error_tiny_output <- if (!tiny$ok && grepl("max_output_mb|exceeds", tiny$error, ignore.case = TRUE)) "PASS" else paste0("FAIL: ", tiny$error)

      missing <- .ia_error(importABCatlas::fetch_data(
        data = dat, cell_ids = simple_ids, genes = "__IMPORTABCATLAS_MISSING_GENE__",
        data_type = primary_type, return_type = "matrix", allow_downloads = TRUE))
      row$error_missing_gene <- if (!missing$ok && grepl("gene|matches", missing$error, ignore.case = TRUE)) "PASS" else paste0("FAIL: ", missing$error)

      duplicate <- .ia_error(importABCatlas::fetch_data(
        data = dat, cell_ids = rep(simple_ids[[1L]], 2L), genes = simple_gene$query,
        data_type = primary_type, return_type = "matrix", allow_downloads = TRUE))
      row$error_duplicate_cells <- if (!duplicate$ok && grepl("duplicate", duplicate$error, ignore.case = TRUE)) "PASS" else paste0("FAIL: ", duplicate$error)

      invalid <- .ia_error(importABCatlas::fetch_data(
        data = dat, cell_ids = simple_ids, genes = simple_gene$query,
        data_type = "invalid", return_type = "matrix", allow_downloads = TRUE))
      row$error_invalid_data_type <- if (!invalid$ok && grepl("must be one of|raw.*log2|data_type", invalid$error, ignore.case = TRUE)) "PASS" else paste0("FAIL: ", invalid$error)

      if (isTRUE(check$ok) && !is.na(row$simple_download_bytes) && row$simple_download_bytes > 0) {
        old_reserve <- getOption("importABCatlas.disk_reserve_gb")
        options(importABCatlas.disk_reserve_gb = 1e12)
        disk <- .ia_error(importABCatlas::fetch_data(
          data = dat, cell_ids = simple_ids, genes = simple_gene$query,
          data_type = primary_type, return_type = "matrix", allow_downloads = TRUE))
        options(importABCatlas.disk_reserve_gb = old_reserve)
        row$error_insufficient_disk <- if (!disk$ok && grepl("Insufficient disk space", disk$error, ignore.case = TRUE)) "PASS" else paste0("FAIL: ", disk$error)
      } else row$error_insufficient_disk <- "SKIP_NO_EXACT_COLD_DOWNLOAD"

      after_expression <- .ia_expression_cache_files(cache_dir)
      row$no_expression_download_before_fetch <- if (identical(before_expression, after_expression)) "PASS" else "FAIL_EXPRESSION_APPEARED_DURING_PLAN_OR_ERROR_TEST"

      # Complex plan is always tested, even when the source file is above the
      # configured download cap. It mixes identifiers/symbols and uses two
      # feature-matrix groups when those are visible in the metadata preview.
      complex_cell <- .ia_complex_cells(dat$cell_metadata, complex_cells)
      complex_gene <- if (authoritative_axis_route) {
        list(query = NULL, expected = NULL)
      } else {
        .ia_gene_query(dat$gene_data, complex_genes, mix_symbols = TRUE)
      }
      row$complex_mode <- if (authoritative_axis_route) {
        paste0(complex_cell$mode, "_authoritative_measured_axis")
      } else {
        complex_cell$mode
      }
      start <- proc.time()[["elapsed"]]
      cp <- .ia_error(importABCatlas::plan_fetch(
        data = dat, cell_ids = complex_cell$ids, genes = complex_gene$query,
        data_type = primary_type))
      row$complex_plan_seconds <- .ia_elapsed(start)
      complex_can_fetch <- FALSE; cp_download <- NA_real_; cp_check <- list(ok = FALSE, error = "plan failed")
      if (!cp$ok) {
        row$complex_plan_status <- "FAIL"; row$complex_plan_error <- cp$error
      } else {
        complex_plan <- cp$value; row$complex_plan_status <- "PASS"
        file_rows[[length(file_rows) + 1L]] <- .ia_plan_rows(complex_plan, dataset, "complex_primary")
        cp_check <- .ia_validate_plan(complex_plan)
        cp_download <- .ia_bytes(complex_plan$total_download_bytes)
        row$complex_sources <- length(unique(vapply(complex_plan$files, function(f)
          paste(f$directory, f$file_name, sep = "/"), character(1))))
        complex_can_fetch <- isTRUE(run_fetch) && isTRUE(cp_check$ok) &&
          !is.na(cp_download) && cp_download <= max_fetch_bytes
        if (!cp_check$ok) row$complex_fetch_status <- paste0("SKIP_PLAN_SIZE: ", cp_check$error)
        else if (!isTRUE(run_fetch)) row$complex_fetch_status <- "SKIP_DISABLED"
        else if (is.na(cp_download) || cp_download > max_fetch_bytes) row$complex_fetch_status <- "SKIP_SIZE_CAP"
      }

      # Safe fetch decision based only on exact additional selected-file bytes.
      can_fetch <- isTRUE(run_fetch) && isTRUE(check$ok) &&
        !is.na(row$simple_download_bytes) && row$simple_download_bytes <= max_fetch_bytes
      if (!isTRUE(run_fetch)) row$fetch_skip_reason <- "run_fetch=FALSE"
      else if (!isTRUE(check$ok)) row$fetch_skip_reason <- "exact selected-file size unavailable/invalid"
      else if (is.na(row$simple_download_bytes)) row$fetch_skip_reason <- "additional download size unknown"
      else if (row$simple_download_bytes > max_fetch_bytes)
        row$fetch_skip_reason <- sprintf("%.2f GiB exceeds max_fetch_gb=%.2f",
                                         row$simple_download_bytes / 1024^3, max_fetch_gb)

      if (can_fetch) {
        start <- proc.time()[["elapsed"]]
        fetched <- .ia_error(importABCatlas::fetch_data(
          data = dat, cell_ids = simple_ids, genes = simple_gene$query,
          data_type = primary_type, return_type = "matrix", allow_downloads = TRUE))
        row$simple_fetch_seconds <- .ia_elapsed(start)
        if (!fetched$ok) {
          row$simple_fetch_status <- "FAIL"; row$simple_fetch_error <- fetched$error
        } else {
          simple_result <- fetched$value
          axis_error <- .ia_validate_matrix_result(simple_result, simple_ids, simple_gene$expected)
          row$simple_axes_status <- if (is.na(axis_error)) "PASS" else paste0("FAIL: ", axis_error)
          row$simple_fetch_status <- if (is.na(axis_error)) "PASS" else "FAIL_AXES"

          start <- proc.time()[["elapsed"]]
          warm <- .ia_error(importABCatlas::fetch_data(
            data = dat, cell_ids = simple_ids, genes = simple_gene$query,
            data_type = primary_type, return_type = "matrix", allow_downloads = FALSE))
          row$simple_warm_seconds <- .ia_elapsed(start)
          if (!warm$ok) row$simple_warm_status <- paste0("FAIL: ", warm$error)
          else {
            warm_result <- warm$value
            row$simple_warm_status <- if (isTRUE(all.equal(simple_result$expression, warm_result$expression,
                                                            check.attributes = TRUE))) "PASS_IDENTICAL" else "FAIL_NOT_IDENTICAL"
          }

          # Optional R object construction uses the already cached expression file.
          if (isTRUE(test_object_types)) {
            if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
              sce <- .ia_error(importABCatlas::fetch_data(
                data = dat, cell_ids = simple_ids, genes = simple_gene$query,
                data_type = primary_type, return_type = "sce", allow_downloads = FALSE))
              row$sce_status <- if (sce$ok && inherits(sce$value, "SingleCellExperiment")) "PASS" else paste0("FAIL: ", sce$error)
            } else row$sce_status <- "SKIP_DEPENDENCY"

            if (requireNamespace("SeuratObject", quietly = TRUE) &&
                (is.null(simple_gene$expected) ||
                 !any(grepl("[_|]", simple_gene$expected)))) {
              seu <- .ia_error(importABCatlas::fetch_data(
                data = dat, cell_ids = simple_ids, genes = simple_gene$query,
                data_type = primary_type, return_type = "seurat", allow_downloads = FALSE))
              row$seurat_status <- if (seu$ok && inherits(seu$value, "Seurat")) "PASS" else paste0("FAIL: ", seu$error)
            } else row$seurat_status <- "SKIP_DEPENDENCY_OR_UNSAFE_FEATURE_IDS"

            if (is_spatial && requireNamespace("SpatialExperiment", quietly = TRUE)) {
              sp <- .ia_error(importABCatlas::fetch_data(
                data = dat, cell_ids = simple_ids, genes = simple_gene$query,
                data_type = primary_type, return_type = "spatial", allow_downloads = FALSE))
              row$spatial_status <- if (sp$ok && inherits(sp$value, "SpatialExperiment")) "PASS" else paste0("FAIL: ", sp$error)
            } else row$spatial_status <- if (is_spatial) "SKIP_DEPENDENCY" else "SKIP_NONSPATIAL"
          }

        }
      } else row$simple_fetch_status <- "SKIP_SIZE_OR_CONFIGURATION"

      if (isTRUE(complex_can_fetch)) {
        start <- proc.time()[["elapsed"]]
        cf <- .ia_error(importABCatlas::fetch_data(
          data = dat, cell_ids = complex_cell$ids, genes = complex_gene$query,
          data_type = primary_type, return_type = "matrix", allow_downloads = TRUE))
        row$complex_fetch_seconds <- .ia_elapsed(start)
        if (!cf$ok) {
          row$complex_fetch_status <- "FAIL"; row$complex_fetch_error <- cf$error
        } else {
          complex_result <- cf$value
          axis_error <- .ia_validate_matrix_result(complex_result, complex_cell$ids, complex_gene$expected)
          row$complex_axes_status <- if (is.na(axis_error)) "PASS" else paste0("FAIL: ", axis_error)
          row$complex_fetch_status <- if (is.na(axis_error)) "PASS" else "FAIL_AXES"
        }
      }

      # Alternate fetch is optional because it may require a second large file.
      if (isTRUE(fetch_both_data_types) && identical(row$alternate_plan_status, "PASS")) {
        alt_check <- .ia_validate_plan(alternate_plan)
        alt_download <- .ia_bytes(alternate_plan$total_download_bytes)
        if (!alt_check$ok) {
          row$alternate_fetch_status <- paste0("SKIP_PLAN_SIZE: ", alt_check$error)
        } else if (is.na(alt_download) || alt_download > max_fetch_bytes) {
          row$alternate_fetch_status <- "SKIP_SIZE_CAP"
        } else {
          start <- proc.time()[["elapsed"]]
          af <- .ia_error(importABCatlas::fetch_data(
            data = dat, cell_ids = simple_ids, genes = simple_gene$query,
            data_type = alternate_type, return_type = "matrix", allow_downloads = TRUE))
          row$alternate_fetch_seconds <- .ia_elapsed(start)
          if (!af$ok) {
            row$alternate_fetch_status <- "FAIL"; row$alternate_fetch_error <- af$error
          } else {
            alternate_result <- af$value
            axis_error <- .ia_validate_matrix_result(alternate_result, simple_ids, simple_gene$expected)
            row$alternate_fetch_status <- if (is.na(axis_error)) "PASS" else paste0("FAIL: ", axis_error)
          }
        }
      } else if (isTRUE(fetch_both_data_types)) {
        row$alternate_fetch_status <- "SKIP_ALTERNATE_UNAVAILABLE"
      } else row$alternate_fetch_status <- "SKIP_DISABLED"
      TRUE
    }, ia_load_stop = function(e) FALSE, ia_plan_stop = function(e) FALSE,
       error = function(e) {
         if (is.na(row$load_error) || !nzchar(row$load_error)) row$load_error <<- conditionMessage(e)
         FALSE
       })

    # Release R references before deleting the complete isolated cache.
    dat <- simple_result <- warm_result <- complex_result <- alternate_result <- simple_plan <- complex_plan <- alternate_plan <- NULL
    invisible(gc())
    info <- .ia_dir_info(cache_dir)
    row$cache_files_deleted <- if (isTRUE(cleanup)) info$files else 0L
    row$cache_bytes_deleted <- if (isTRUE(cleanup)) info$bytes else 0
    if (isTRUE(cleanup)) {
      unlink(cache_dir, recursive = TRUE, force = TRUE)
      row$cleanup_status <- if (!dir.exists(cache_dir)) "PASS_DELETED_DATASET_CACHE" else "FAIL_DIRECTORY_REMAINS"
    } else row$cleanup_status <- "SKIP_RETAINED"

    # Normalize list fields to scalar data-frame values.
    scalar <- function(x) {
      if (is.null(x) || !length(x)) return(NA)
      if (length(x) > 1L) return(paste(x, collapse = " | "))
      x
    }
    row <- lapply(row, scalar)
    summary_rows[[length(summary_rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
    validation <- checkpoint()
    cat(sprintf("  load=%s %.2fs | plan=%s | fetch=%s | cleanup=%s (%s)\n",
                row$load_status, as.numeric(row$load_seconds), row$simple_plan_status,
                row$simple_fetch_status, row$cleanup_status, .ia_human_bytes(as.numeric(row$cache_bytes_deleted))))
  }

  validation <- checkpoint()
  validation$metadata_runtime <- stats::setNames(validation$summary$load_seconds,
                                                  validation$summary$dataset)
  saveRDS(validation, paste0(output_prefix, ".rds"))
  cat("\nValidation complete. Results:\n",
      paste0(output_prefix, "_summary.csv\n"),
      paste0(output_prefix, "_files.csv\n"),
      paste0(output_prefix, ".rds\n"), sep = "")
  validation
}
