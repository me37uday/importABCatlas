# Native R release and validation checklist

This checklist defines the minimum validation expected before the
native-R implementation of importABCatlas is treated as release-ready or
merged into a release branch.

The package has two complementary validation layers: portable package
validation using local fixtures and CI, and real-data acceptance against
the Allen Brain Cell Atlas release. Passing one layer does not replace
the other.

## 1. Repository state

-   [ ] Working tree is clean.
-   [ ] Intended branch and commit hash are recorded.
-   [ ] No generated build/check artifacts are tracked.
-   [ ] No obsolete Python or reticulate runtime dependency has been
    introduced.
-   [ ] README installation instructions point to the intended branch or
    release.
-   [ ] `git diff --check` is clean.

Record the branch, commit, validation date, and validator.

## 2. Unit and fixture tests

Run the complete package test suite with `testthat::test_local(".")`.

-   [ ] Zero failed tests.
-   [ ] Zero unexpected warnings.
-   [ ] Any skips are understood and documented.
-   [ ] Matrix materialization passes.
-   [ ] Seurat materialization passes.
-   [ ] SingleCellExperiment materialization passes.
-   [ ] SpatialExperiment materialization passes.

The fixture suite should remain independent of large Atlas expression
downloads.

## 3. Package build and check

Run `R CMD build .` followed by `R CMD check --no-manual` on the built
tarball.

-   [ ] Package builds successfully.
-   [ ] Package installs successfully.
-   [ ] Zero `R CMD check` errors.
-   [ ] Zero `R CMD check` warnings.
-   [ ] Any notes are reviewed and documented.
-   [ ] Vignette builds successfully.
-   [ ] Examples pass.

## 4. Cross-platform CI

The native-R GitHub Actions workflow must complete on:

-   [ ] Ubuntu
-   [ ] macOS
-   [ ] Windows

Review failures rather than weakening tests solely to make CI green.
Platform-specific filesystem paths, temporary directories, line endings,
and system dependencies must be treated as portability concerns.

## 5. Real-data route acceptance

Validate every supported route against the intended Allen release.

For each route verify:

-   [ ] `load_data()` succeeds.
-   [ ] Cell identifiers are preserved exactly as character data.
-   [ ] Gene selection is valid for the selected expression
    representation.
-   [ ] `plan_fetch()` selects the expected representation and files.
-   [ ] Required files can be resolved from the validation cache.
-   [ ] `fetch_data()` returns the requested cells and genes.
-   [ ] Row order is correct.
-   [ ] Column order is correct.
-   [ ] Expression provenance is correct.
-   [ ] H5AD layout is supported.
-   [ ] Applicable object materializers succeed.

Record the aggregate result and commit hash.

## 6. Offline / zero-download acceptance

After the required validation files are cached, rerun acceptance without
permitting additional downloads.

-   [ ] Every supported route passes in no-download mode.
-   [ ] Planned additional download bytes are zero.
-   [ ] No route silently accesses a different representation or file.
-   [ ] The validation result is tied to the exact package commit.

## 7. Feature-universe audit

For each expression route compare metadata gene identifiers with the
H5AD feature universe and classify the relationship (`EXACT`,
`H5AD_SUBSET`, or `DIFFERENT` as appropriate).

-   [ ] Every route is classified.
-   [ ] H5AD duplicate identifiers are checked.
-   [ ] Metadata-only identifiers are quantified.
-   [ ] H5AD-only identifiers are quantified.
-   [ ] Known upstream identifier malformations are documented rather
    than silently rewritten.
-   [ ] Gene requests unavailable in a selected expression matrix fail
    clearly rather than being silently dropped.

## 8. Known special cases

Explicitly recheck routes representing previously discovered edge cases:

-   [ ] AgingMouse gene-identifier irregularities.
-   [ ] HMBA_Macaque accession-qualified / duplicate gene identifiers.
-   [ ] Zhuang spatial cell identifiers.
-   [ ] MERFISH long cell identifiers.
-   [ ] MERFISH_imputed expression feature subset.
-   [ ] Dense H5AD expression matrices.
-   [ ] Raw versus upstream log2/imputed expression semantics.

Extend this list whenever a new regression class is discovered.

## 9. Documentation

-   [ ] README matches the current API.
-   [ ] README contains no obsolete Python setup instructions.
-   [ ] Installation instructions were tested.
-   [ ] `list_datasets()`, `load_data()`, `plan_fetch()`, and
    `fetch_data()` examples reflect the current interface.
-   [ ] Raw/log2 semantics are documented.
-   [ ] Cache/offline behavior is documented.
-   [ ] Optional output-package requirements are documented.
-   [ ] Collaborator testing instructions are current.

## 10. Release evidence

Record the package commit, Allen release, unit-test result,
`R CMD check` result, Ubuntu/macOS/Windows CI status, real-data route
results, offline acceptance, feature-universe audit, and validation
report location.

A release decision should be based on this evidence, not solely on
whether the package installs successfully.
