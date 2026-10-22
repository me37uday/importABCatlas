# Audit and validation record

Initial audit date: 2026-10-08; comprehensive-candidate update: 2026-10-16. Inputs: the uploaded `importABCatlas-main.zip` and `importABCatlas_status.xlsx`. The Google manuscript was not accessible. The code ZIP, not an independently cloned current GitHub commit, is the repair baseline. Primary documentation was checked online; live S3 file access from the code execution environment failed.

## Workbook interpretation

The status table is in the uploaded workbook's Sheet2; the unfilled timing layout is in Sheet1. The emails use different sheet numbering. Ten family rows are listed. Existing-function colouring marks five families: WMB, Aging Mouse, WHB, PMDBS and HMBA. Of those, only WHB and PMDBS have working load status in the supplied sheet. The fetch column is mostly blank, and one entry says yet to be reviewed; this is not evidence that their expression fetches passed. No measured runtimes are present.

The replacement audit workbook translates the status information into explicit words. The supplied workbook is retained unchanged in the parent bundle. The expanded 23-route count is an importer design count separating components and modalities, not a correction claiming that Allen publishes exactly 23 biological studies.

## Confirmed source defects and changes

The original DESCRIPTION named the package `importAllen`, contrary to the requested/importABCatlas interface, and contained unfinished author/licence metadata. The replacement uses the intended package name and declares dependencies, but leaves real authorship and licensing as explicit owner decisions.

Aging Mouse read WHB human gene metadata. It now uses the documented WMB mouse annotations. PMDBS likewise read WHB genes; it now uses its own `ASAP-PMDBS-10X/gene` table. Wrong gene universes can cause plausible-looking but incorrect selections even when a metadata loader appears to run.

The WMB/WHB taxonomy pivot discarded original cluster aliases and reconstructed them from row positions. The replacement preserves the source aliases and checks join cardinality. General joins no longer silently multiply cells or overwrite colliding annotations. Exact cell and gene order is checked against H5AD indices.

HMBA directory capitalisation is inconsistent across official release notes and notebooks. It was too strong to attribute the reported HMBA failure solely to one guessed spelling. The fix resolves a unique actual manifest directory, records it and refuses absent or ambiguous matches; the true live failure cause still requires the live log.

The original setup depended on a hard-coded Python 3.10 installation and mixed environment assumptions. The replacement provisions a private Python 3.11 environment automatically with reticulate/uv, verifies every dependency, and runs the backend out of process. Advanced external-Python use is opt-in. The Allen SDK is pinned to the tested commit rather than a moving default branch.

The expression implementation was rewritten to plan selected source shards, preserve string IDs, validate sparse/dense H5AD selection, retain native feature universes, distinguish transformed/imputed expression, and construct optional R frameworks. This is a significant implementation change and therefore increases the importance of independent live reference comparisons; it is not a small patch presumed safe without testing.

## Executed evidence

The comprehensive candidate test log records **108 passed, 3 skipped**. The passing tests include all 23 controlled route fixtures, exact selected-file sizing through local/manifest/S3 metadata paths, insufficient-space and unknown-size safeguards, offline-cache errors, simple and multi-source complex fetches, order preservation, raw/log2 semantics, cache reuse and corrupt-H5AD handling. See the machine-readable environment and scope files rather than relying on this prose alone.

Synthetic HDF5 fixtures exercised dense, CSR and CSC layouts; order permutations; huge numeric-looking IDs; duplicate/missing genes and cells; correct taxonomy aliases; many-to-one joins; explicit file mappings; imputed counts rejection; native gene resolution; request/response handling; full-file-download consent; and all 23 route configurations. Benchmark safety tests reject synthetic data as manuscript timings. A separate test of the candidate-installation helper uses a disposable local Git repository and verifies dry-run, backup, branch and history safeguards.

The three explicit skips document R integration, live Allen access, and independent AnnData comparison. The independent comparison code is supplied, but AnnData was unavailable locally. Successful fixture routing means a controlled mock manifest with the expected schema works; it does not independently verify that the upstream manifest currently has that schema.

## Not executed / no claim made

R parsing, R package installation, roxygen2, R tests, Seurat/SCE/SpatialExperiment constructors, R CMD check, GitHub Actions, actual Allen metadata/expression downloads, independent real-source comparison, all-shard validation and biological cold/warm benchmarks were not executed here. The execution container had no R and no outbound DNS. Attempts to install R failed because repositories could not be resolved. Python fixture execution was feasible without network. Live runtime cells remain NOT_MEASURED, never zero.

The current candidate implements native HMBA species access by filtering the shared cell table with donor `species_genus` before preview selection and selecting exact `human_gene`, `macaque_gene`, or `marmoset_gene` metadata. Patch-seq uses only its declared dataset directory and an unambiguous single-file fallback. These repairs pass controlled tests but still require the supplied live harness before being reported as verified.

## What constitutes a release pass

Run the local R tests and R CMD check with optional dependencies installed. Then exercise each applicable dataset/representation through metadata loading, plan generation, authorised fetch, independent AnnData source comparison, and every promised R object constructor. Check exact IDs, axes, values, donor joins, taxonomy aliases and coordinates. Cover more than one file group and both biological cohorts when present. Explicitly classify unavailable representations and non-applicable cohorts rather than silently skipping them.

A 404, changed field, non-unique taxonomy join, or unavailable raw/log2 product requires investigation. `FAIL_OR_UNAVAILABLE` is intentionally not a pass. Record the exact manifest, SDK commit/version, package commit and environment. Re-run after changing the registry or reader. Archive genuine output and reviewer-approved exceptions before switching status to validated.
