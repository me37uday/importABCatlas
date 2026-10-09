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
