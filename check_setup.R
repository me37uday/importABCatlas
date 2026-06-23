# =============================================================================
# check_setup.R  --  Smoke test for importABCatlas after the recent changes.
#
# Run it stage by stage (source the file, or step through in RStudio). The
# stages get progressively heavier: 0-1 are instant and offline, 2-3 hit S3 for
# a few small metadata files, and 4 (opt-in) downloads a large expression
# matrix. Flip the flags in the CONFIG block to control how far it goes.
# =============================================================================

## ---- CONFIG -----------------------------------------------------------------

# Path to the cloned repo (so devtools::load_all() can find the package).
REPO <- "/u/f/fcaretti/importABCatlas"

# Where atlas data lives.
DOWNLOAD_BASE <- "/data/abc_atlas_2"

# Stage 4 downloads a full expression matrix (can be several GB). Off by default.
RUN_FETCH <- FALSE

## ---- Stage 0: load the package ---------------------------------------------

library(devtools)
load_all(REPO)

stopifnot(
  is.function(load_data),
  is.function(fetch_data),
  is.function(abc_download_dir),
  is.function(abc_cache_connect)
)
cat("Stage 0 OK: package loaded, all functions exported.\n\n")

## ---- Stage 1: download-folder resolver (instant, no network) ---------------

# Precedence: explicit arg > option > ABC_ATLAS_CACHE env var > per-user cache.
cat("Default cache dir:        ", abc_download_dir(), "\n")
cat("Explicit arg resolves to: ", abc_download_dir(DOWNLOAD_BASE), "\n")

options(importABCatlas.download_base = DOWNLOAD_BASE)
stopifnot(abc_download_dir() == normalizePath(DOWNLOAD_BASE))
cat("Stage 1 OK: option override works; resolved path is absolute.\n\n")
# (Leaving the option set means the rest of the script uses DOWNLOAD_BASE
#  even without passing download_base explicitly.)

## ---- Stage 2: load metadata (small downloads on first run) -----------------

# First call: provisions Python (reticulate/uv) if needed, downloads any
# metadata CSVs not already present (abc_atlas_access caches them), assembles it.
t <- system.time(WHB <- load_data("WHB"))
cat(sprintf("Stage 2: load_data(\"WHB\") took %.1fs\n", t["elapsed"]))

stopifnot(
  is.list(WHB),
  all(c("cell_metadata", "gene_data", "unique_values") %in% names(WHB)),
  nrow(WHB$cell_metadata) > 0,
  "gene_symbol" %in% colnames(WHB$gene_data)
)
cat(sprintf("Stage 2 OK: %d cells, %d genes.\n\n",
            nrow(WHB$cell_metadata), nrow(WHB$gene_data)))

## ---- Stage 3: offline mode reuses already-downloaded data ------------------

# Now that stage 2 has cached the raw metadata files, offline = TRUE must rebuild
# the same table without contacting S3. (Each call re-runs the joins; there is no
# separate processed cache -- we rely on abc_atlas_access's raw-file cache.)
WHB_off <- load_data("WHB", offline = TRUE)

# Compare data content, not attributes (each build carries its own reticulate
# pandas.index reference, so identical() on the raw frames would differ).
stopifnot(
  nrow(WHB_off$cell_metadata) == nrow(WHB$cell_metadata),
  identical(colnames(WHB_off$cell_metadata), colnames(WHB$cell_metadata)),
  identical(WHB_off$cell_metadata$cell_label, WHB$cell_metadata$cell_label)
)
cat("Stage 3 OK: offline = TRUE rebuilt the same table with no S3 access.\n\n")

## ---- Stage 4 (opt-in): fetch a Seurat object & check sparsity --------------

if (isTRUE(RUN_FETCH)) {
  cat("Stage 4: fetching expression (this downloads a large matrix)...\n")
  filters <- list(
    feature_matrix_label = "WHB-10Xv3-Nonneurons",
    cluster              = "VendV_17"
  )
  obj <- fetch_data(
    metadata  = WHB$cell_metadata,
    filters   = filters,
    gene_data = WHB$gene_data
  )

  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  cat("  Seurat object:", ncol(counts), "cells x", nrow(counts), "genes\n")
  cat("  counts class:", class(counts), "\n")

  # The whole point of the fix: counts must be a SPARSE matrix, never dense.
  stopifnot(inherits(counts, "dgCMatrix"))
  cat("Stage 4 OK: counts is a sparse dgCMatrix (no densification).\n")
} else {
  cat("Stage 4 SKIPPED (set RUN_FETCH <- TRUE to test the full fetch).\n")
}

cat("\nAll requested stages passed.\n")
