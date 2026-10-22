# importABCatlas 0.2.0.9007 (automatic Python and wet-lab quick start)

- Make a private Python 3.11 environment the default. It is resolved and installed automatically with reticulate/uv in a clean helper R process, so users no longer need to set a Python path, edit `.Renviron`, manage Conda, or depend on Homebrew Python.
- Ignore ambient Python variables in normal managed mode, preventing stale interpreter settings from silently selecting an incomplete environment.
- Add `setup_importABCatlas()` for explicit one-time setup and self-repair. The first backend call also triggers setup automatically.
- Keep deliberate external-Python use as an advanced opt-in through `setup_environment(python=...)`.
- Add health checks, requirement-version tracking, a setup lock for concurrent R sessions, persistent setup logs, and automatic recreation when the managed cache is removed or the package requirements change.
- Rewrite the README and beginner setup instructions around a short R-only workflow suitable for wet-lab users.

# importABCatlas 0.2.0.9006 (feature-axis one-to-one reconciliation)

- Fix native HMBA Macaque feature matching when duplicated NCBI identifiers and duplicated gene symbols coexist. Metadata now preserves duplicate occurrence/group-size information, while matching prefers Ensembl, NCBI, transcript and other stable identifiers before symbols. `ENSEMBL:` and `NCBIGene:` prefixes are reconciled conservatively.
- Add the native HMBA identifier columns `ensembl_gene_identifier`, `ncbi_gene_identifier`, and `ncbi_gene_symbol` to the supported H5AD/metadata alias vocabulary.
- Resolve an ambiguous alias only after trying all stronger stable identifiers. If both metadata and H5AD retain the same duplicated stable identifier, rows are aligned one-to-one by deterministic occurrence only when duplicate-group sizes agree; rows are never aggregated.
- Fix dense imputed-MERFISH extraction by canonicalising metadata supersets to one metadata row per authoritative H5AD feature before indexing. The selected dense HDF5 column index is now strictly increasing and duplicate-free, as required by h5py.
- Reject explicit requests in which two metadata identifiers resolve to the same measured H5AD feature, with guidance to run `genes = NULL` and reuse the returned measured identifiers.
- Extend provenance with one-to-one feature-axis reconciliation counts, collapsed metadata-alias rows, match methods, and examples.
- Add focused live regression coverage for explicit disambiguated HMBA Macaque identifiers and for imputed-MERFISH authoritative-axis simple/complex queries.
- Keep the two routes marked `LIVE_FETCH_RETEST_PENDING_9006` until the supplied live test passes against the retained 7.53 GB and 50.2 GB H5AD files.

# importABCatlas 0.2.0.9004 (clean-install and final large-fetch validation candidate)

- Repackage the 0.2.0.9003 code as a clean source tree with no installed lazy-load databases.
- Add `INSTALL_CLEAN_9004.R`, which installs into a dedicated versioned R library from a fresh `Rscript --vanilla` process and verifies the namespace, package version, 23-route catalogue, and optional Python dependencies.
- Add `TEST_REMAINING_LARGE_FETCHES_9004.R`, which refuses to use a different loaded namespace, performs a synthetic R/Python smoke test before large downloads, validates checkpoint/package compatibility, and tests the seven size-gated expression routes with resumable downloads and exact cleanup.
- Keep all data-access behavior, exact-file planning, disk-space checks, retry/resume logic, and release-20260711 route capabilities from 0.2.0.9003 unchanged.
- This version exists to prevent stale or overwritten installed-package state from being mistaken for a source-code defect.

# importABCatlas 0.2.0.9003 (capability and validation clarity)

- All 23 bundled routes now declare a modality, available expression representations, a default representation, and release-scoped live-validation status.
- `list_datasets()` exposes `modality`, `available_data_types`, `validation_release`, and `validation`.
- `MERFISH_imputed` and `HMBA_PatchSeq` are explicitly log2-only because no raw expression asset is published in the pinned Allen manifest; unsupported raw requests fail before the backend or any expression download.
- Other bundled routes advertise both raw and log2 only where both exact manifest assets were observed.
- Improve invalid `data_type` errors so they identify the argument and allowed values.
- Repair the comprehensive validation harness version check, unsupported-representation classification, and invalid-data-type assertion.
- Add a consolidated release-20260711 validation table: every route passed live metadata loading and exact-file planning; 16 routes completed live simple and complex fetches, while 7 large routes were not downloaded because of the configured size cap.
- No dataset route was removed: every route maps to a published Allen asset. Unpublished representations are omitted rather than presented as working options.

# importABCatlas 0.2.0.9002 (live-validation outcome repair)

- Avoid the fragile bucket-wide `releases/` listing when a YYYYMMDD release is pinned; the exact manifest object is loaded directly.
- Add package-level retry/backoff and configurable connect/read timeouts for transient Allen/S3 failures.
- Add resumable `.part` downloads with S3 Range requests, exact byte-count checks, checksum verification, visible progress/rate/ETA, and retained partial bytes after an interrupted transfer.
- `plan_fetch()` continues to report exact selected-object sizes; it now also reports resumable partial bytes and the remaining additional download. Directory totals are never substituted.
- Treat genuinely empty optional manifest categories (`expression_matrices`, `image_volumes`, or `mapmycells`) as empty in `list_files()` while propagating network and malformed-manifest errors.
- Disambiguate duplicated HMBA native gene identifiers without aggregating rows. The original identifier is retained in `source_gene_identifier`; ambiguous original-ID queries fail with the generated unique identifiers.
- Declare route-specific expression defaults. `MERFISH_imputed` and `HMBA_PatchSeq` default to `log2`; all other bundled routes default to `raw`. Explicit unavailable representations fail before downloading.
- Pin `setup_environment()` by default to the official `abc_atlas_access` commit used in the user live-validation environment, and probe/install the explicit S3 downloader dependencies.
- Expand regression coverage for optional empty categories, route defaults, exact-size retries, duplicate genes, partial-download resume, and completed-partial finalization.

# importABCatlas 0.2.0.9001 (comprehensive validation candidate)

- `plan_fetch()` now obtains exact selected-object sizes from verified local files, manifest file attributes, or S3 object metadata; directory totals are never substituted.
- `fetch_data()` performs an exact pre-download disk-space check with configurable safety factor and free-space reserve, and rejects unknown sizes by default.
- Verified cached expression files can be read without renewed download consent; uncached files still require explicit consent.
- Native HMBA Human, Macaque and Marmoset routes use shared cell metadata, donor-based species filtering before previewing, and exact species gene tables. Example-expression tables cannot be mistaken for gene metadata.
- Patch-seq and all expression routing match raw/log2 suffixes case-insensitively and retain the single-file-within-dataset safety rule.
- Added synthetic simple, multi-file complex, log2, cache, exact-size, insufficient-space, unknown-size and pre-download error regression tests.
- Added `scripts/test_all_functions.R` and `TEST_ALL_FUNCTIONS.R`, which checkpoint results and delete each isolated dataset cache after testing.

# importABCatlas 0.2.0.9000 (source release candidate)

- Correct package identity and separate R, Python, cache and documentation responsibilities.
- Replace repeated dataset implementations with a manifest-driven registry and shared engine.
- Add intended routes for consensus, developing mouse, MERFISH measured/imputed, Zhuang1-4, HMBA native/spatial/Patch-seq and SEA-AD Multiregion/CaH.
- Repair native gene-source selection and taxonomy alias handling.
- Add explicit download consent, exact identifier preservation, safe joins, selective dense/CSR/CSC reads, size checks and provenance.
- Add Seurat/SCE/SpatialExperiment/matrix choices, with distinct transformed and imputed semantics.
- Add synthetic fixtures, regression tests, supplied R/cross-platform CI/live validators, safe branch migration helper, and paired benchmark scripts.
- Validation remains local-fixture-only; R, live atlas and scientific benchmarks are release blockers.

## Live manifest/fetch routing fix

- Treat Allen SDK's explicit empty optional sub-directory response as an empty file category without hiding unrelated manifest errors.
- Restrict expression planning to the selected dataset's expression directories; unrelated taxonomy/metadata directories can no longer break MERFISH or other fetch plans.
- Validate explicit `file_map` entries directly against their named expression directory.

## 0.2.0.9000 interpreter hotfix

- Resolve the Python executable exactly once per backend call and pass it directly to `system2()`.
- Guarantee that an explicit `python=` argument has highest precedence, followed by `IMPORTABCATLAS_PYTHON`, the R option, reticulate, the legacy variable, and PATH.
- Backend errors now report the exact Python executable used.
- Align beginner documentation on `IMPORTABCATLAS_PYTHON`.

## Focused PMDBS / Aging Mouse live-schema repair

- Preserve configured virtual-environment Python executable paths; do not resolve venv symlinks to the base Homebrew interpreter.
- Treat Allen SDK empty optional sub-directories as empty file categories while preserving genuine manifest errors.
- PMDBS uses `ASAP-PMDBS-10X/gene` and native PMDBS cluster labels. The very large SEA-AD and Siletti/WHB MapMyCells files remain optional reference mappings and are not forced into every load.
- Aging Mouse joins native `cell_cluster_annotations` and cell-level WMB `cell_cross_mapping_annotations` by exact `cell_label`, suffixing overlapping columns rather than overwriting native annotations.
- Aging Mouse retains the WMB gene universe because the observed `Zeng-Aging-Mouse-10Xv3` metadata manifest does not expose a `gene` metadata table.

## Candidate update: fetch planning and progress
- Restores visible Allen download progress/ETA in interactive R sessions.
- Preserves virtualenv Python symlinks in both runtime resolution and setup_environment().
- Keeps both raw and log2 expression modes; raw remains the default and only the selected representation is routed.
- plan_fetch() now prints the exact expression manifest entries that would be used, requested output dimensions, and Allen directory-size context without downloading expression.

## v8 candidate (2026-10-09)
- Fetch plans no longer show directory totals as exact file sizes. Individual
  source file sizes are reported when independently known; otherwise unknown.
- Added `scripts/test_all_datasets.R`: static-catalogue generic dispatch,
  load/plan-only testing, per-dataset runtime, checkpoints and safe cleanup of
  newly downloaded metadata.
- Imputed MERFISH is planned with log2 in the validation harness.
- HMBA PatchSeq unresolved expression routing remains an explicit live
  validation issue; no unverified filename is fabricated.
