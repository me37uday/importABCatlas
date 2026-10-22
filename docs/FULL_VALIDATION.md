# Comprehensive live validation

`TEST_ALL_FUNCTIONS.R` and `scripts/test_all_functions.R` contain the same release-scoped validation harness. The harness uses the package's static registry and a pinned Allen manifest; it does not depend on `live = TRUE` discovery.

## Evidence carried into 0.2.0.9007

For Allen release `20260711`:

- all 23 routes completed `load_data()`;
- all 23 routes completed `plan_fetch()` with exact selected S3 objects and exact file sizes;
- 21 routes completed real simple and complex expression retrieval, exact axis/order checks, warm-cache equality, and available object-construction checks;
- HMBA Macaque and imputed MERFISH downloaded and opened their real source H5ADs but exposed one-to-one metadata/H5AD feature-axis defects;
- no complete dataset route is absent upstream or removed from the registry;
- raw imputed MERFISH and raw HMBA Patch-seq are not published representations and remain excluded from advertised choices.

Candidate 0.2.0.9007 retains the 0.2.0.9006 feature-axis fixes and adds automatic managed Python. The two route statuses remain pending until the focused live script succeeds.

## One-to-one measured-axis behavior

When `genes = NULL`, metadata rows are matched against the authoritative feature axis of every selected H5AD source. The package selects at most one metadata row per measured feature and then takes the conservative intersection across sources. Stable Ensembl/NCBI/transcript identifiers outrank gene symbols. Provenance records unmatched rows, alias rows collapsed to an already represented feature, selected rows, feature counts, and match methods.

Explicit gene requests remain strict. If two requested identifiers resolve to the same measured feature, or if a requested identifier is absent, the operation stops rather than duplicating or zero-filling expression.

## Safe all-route run

```r
source("TEST_ALL_FUNCTIONS.R")

validation <- test_all_functions(
  test_root = "~/importABCatlas_9007_safe_validation_cache",
  release = "20260711",
  preview_n = 1000L,
  run_fetch = TRUE,
  max_fetch_gb = 2,
  cleanup = TRUE,
  output_prefix = "~/Desktop/importABCatlas_9007_safe_validation"
)
```

## Focused two-route regression

```r
source("TEST_FAILED_DATASETS_9007.R")

result <- test_failed_datasets_9007(
  cache_root = "~/importABCatlas_9004_large_fetch_cache",
  release = "20260711",
  download_missing = FALSE,
  cleanup_expression_after_success = TRUE,
  output_prefix = "~/Desktop/importABCatlas_9007_failed_routes"
)

result$all_selected_routes_pass
```

The focused script reuses the retained large files, discovers the measured axis with `genes = NULL`, asserts the published feature counts, performs strict explicit simple/complex queries using returned identifiers, explicitly tests `.iabdupN` Macaque identifiers, preserves failed files, and deletes the exact H5AD only after success when requested.
