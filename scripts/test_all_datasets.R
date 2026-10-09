# importABCatlas catalogue smoke test: load_data + plan_fetch ONLY.
# Never calls fetch_data. Uses an isolated cache and deletes only NEW metadata
# files created by each dataset. Existing files and expression files are preserved.
# Run source('scripts/test_all_datasets.R') after library(importABCatlas).

.test_cache_files <- function(root) {
  if (!dir.exists(root)) return(character())
  list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
             include.dirs = FALSE, no.. = TRUE)
}

.test_is_metadata <- function(path, root) {
  # Restrict cleanup to the metadata tree; no expression, manifests, or other files.
  root <- normalizePath(root, mustWork = TRUE, winslash = "/")
  candidate <- normalizePath(path, mustWork = TRUE, winslash = "/")
  startsWith(candidate, paste0(root, "/metadata/"))
}

.test_plan_files <- function(plan) {
  if (is.null(plan$files) || !length(plan$files)) return(NA_character_)
  paste(vapply(plan$files, function(f) paste0(f$directory, "/", f$file_name), ""), collapse = " | ")
}

.test_plan_bytes <- function(plan) {
  if (is.null(plan$files) || !length(plan$files)) return(NA_real_)
  fs <- plan$files
  keys <- vapply(fs, function(f) paste0(f$directory, "/", f$file_name), "")
  fs <- fs[!duplicated(keys)]
  vals <- vapply(fs, function(f) if (is.null(f$file_bytes)) NA_real_ else as.numeric(f$file_bytes), 0.0)
  if (anyNA(vals)) NA_real_ else sum(vals)
}

test_all_datasets <- function(
    download_base = "~/abc_atlas_validation_isolated",
    release = "20260711", preview_n = 20L,
    plan_cells = 5L, plan_genes = 3L, data_type = "raw",
    cleanup = TRUE, output_prefix = "importABCatlas_dataset_validation") {
  stopifnot(requireNamespace("importABCatlas", quietly = TRUE))
  data_type <- match.arg(data_type, c("raw", "log2"))
  stopifnot(preview_n > 0L, plan_cells > 0L, plan_genes > 0L)
  download_base <- path.expand(download_base)
  dir.create(download_base, recursive = TRUE, showWarnings = FALSE)
  catalogue <- importABCatlas::list_datasets()  # deliberately NOT live=TRUE
  if (!"dataset" %in% names(catalogue)) stop("Catalogue missing dataset column")
  rows <- vector("list", nrow(catalogue))
  for (i in seq_len(nrow(catalogue))) {
    dataset <- as.character(catalogue$dataset[[i]])
    representation <- if (identical(dataset, "MERFISH_imputed")) "log2" else data_type
    cat(sprintf("\n[%d/%d] %s (%s)\n", i, nrow(catalogue), dataset, representation))
    before <- .test_cache_files(download_base)
    load_status <- "FAIL"; plan_status <- "NOT_RUN"
    load_error <- NA_character_; plan_error <- NA_character_
    elapsed <- NA_real_; plan_elapsed <- NA_real_
    n_cells <- NA_integer_; n_genes <- NA_integer_
    files <- NA_character_; file_bytes <- NA_real_
    download_bytes <- NA_real_; cleaned_files <- 0L; cleaned_bytes <- 0
    dat <- NULL; plan <- NULL
    tryCatch({
      start <- proc.time()[["elapsed"]]
      # Generic dispatch covers Consensus, Zhuang and HMBA species/spatial routes.
      dat <- tryCatch(importABCatlas::load_data(
        dataset, download_base = download_base, release = release,
        preview_n = preview_n), error = function(e) {
          elapsed <<- proc.time()[["elapsed"]] - start
          stop(e)
        })
      elapsed <- proc.time()[["elapsed"]] - start
      n_cells <- nrow(dat$cell_metadata); n_genes <- nrow(dat$gene_data)
      load_status <- "PASS"
      if (!all(c("cell_label") %in% names(dat$cell_metadata)) ||
          !all(c("gene_identifier") %in% names(dat$gene_data)) ||
          n_cells < 1L || n_genes < 1L) {
        plan_status <- "SKIP"; plan_error <- "No valid cell_label/gene_identifier selection"
      } else {
        start_plan <- proc.time()[["elapsed"]]
        plan <- tryCatch(importABCatlas::plan_fetch(
          data = dat,
          cell_ids = head(dat$cell_metadata$cell_label, plan_cells),
          genes = head(dat$gene_data$gene_identifier, plan_genes),
          data_type = representation), error = function(e) {
            plan_elapsed <<- proc.time()[["elapsed"]] - start_plan
            plan_error <<- conditionMessage(e)
            NULL
          })
        if (!is.null(plan)) {
          plan_elapsed <- proc.time()[["elapsed"]] - start_plan
          plan_status <- "PASS"; files <- .test_plan_files(plan)
          file_bytes <- .test_plan_bytes(plan)
          if (!is.null(plan$total_download_bytes))
            download_bytes <- as.numeric(plan$total_download_bytes)
        } else plan_status <- "FAIL"
      }
    }, error = function(e) { load_error <<- conditionMessage(e) }, finally = {
      # Cleanup happens even on failure; only files created by this iteration,
      # only under the isolated metadata tree. No pre-existing cache is removed.
      after <- .test_cache_files(download_base)
      new <- setdiff(after, before)
      if (isTRUE(cleanup) && length(new)) {
        for (f in new) {
          if (!file.exists(f) || !.test_is_metadata(f, download_base)) next
          bytes <- file.info(f)$size
          if (file.remove(f)) {
            cleaned_files <- cleaned_files + 1L
            cleaned_bytes <- cleaned_bytes + if (is.na(bytes)) 0 else bytes
          }
        }
      }
    })
    rows[[i]] <- data.frame(
      dataset = dataset, data_type = representation,
      load_status = load_status, metadata_elapsed_seconds = elapsed,
      preview_cells = n_cells, genes = n_genes,
      plan_status = plan_status, plan_elapsed_seconds = plan_elapsed,
      expression_files = files, selected_file_bytes = file_bytes,
      additional_download_bytes = download_bytes,
      metadata_files_deleted = cleaned_files, metadata_bytes_deleted = cleaned_bytes,
      load_error = load_error, plan_error = plan_error,
      stringsAsFactors = FALSE)
    results <- do.call(rbind, rows[seq_len(i)])
    write.csv(results, paste0(output_prefix, ".csv"), row.names = FALSE, na = "")
    saveRDS(results, paste0(output_prefix, ".rds"))
    cat(sprintf("  load=%s %.2fs | plan=%s | deleted %d new metadata files\n",
                load_status, elapsed, plan_status, cleaned_files))
    rm(dat, plan); invisible(gc())
  }
  results <- do.call(rbind, rows)
  attr(results, "metadata_runtime") <- stats::setNames(results$metadata_elapsed_seconds, results$dataset)
  results
}
