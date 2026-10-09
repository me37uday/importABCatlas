# importABCatlas

## Release candidate 0.2.0.9000 - 8 October 2026

Release-aware access to Allen Brain Cell Atlas metadata and selected expression from R. The package delegates manifests, downloads and cache integrity to the official Allen Python SDK, then validates identifiers, joins and on-disk H5AD selection before creating an R analysis object.

**Validation status: 75 Python fixture tests passed, 3 explicitly skipped. No live Allen downloads, R execution, R CMD check, or biological runtime benchmarks were possible in the delivery environment. This is a source release candidate, not a certified all-dataset production release.** See `validation/python-tests.txt`, `validation/coverage.csv`, and `docs/AUDIT_AND_VALIDATION.md`. Native-species HMBA file discovery particularly needs live confirmation. No benchmark values have been invented.

The registry contains **23 intended expression routes** covering the documented dataset families, with measured and imputed MERFISH, four Zhuang brains, two consensus components, HMBA aligned/native/spatial/Patch-seq routes, and two SEA-AD releases kept distinct. WMB-10X combines its v2/v3/Multi file groups through the expression resolver. A registry entry is not a successful live test.

## Install the source candidate

Unzip the bundle. In R, use the actual path to its `importABCatlas` folder:

```r
install.packages("remotes")
remotes::install_local("/path/to/importABCatlas", dependencies = NA,
                       upgrade = "never")
library(importABCatlas)
setup_environment(install = TRUE)  # explicit internet-dependent installation
example_data("matrix")            # tiny synthetic fixture; no atlas download
```

The package requires R >=4.2 and an explicitly selected Python >=3.10. Python 3.11 is configured in the supplied, as-yet-unexecuted cross-platform CI. `setup_environment()` uses the official SDK's main branch by default; pin a tested commit with `upstream_ref` before publishing. It selects that Python for the current R session. In subsequent sessions run `setup_environment()` or configure `ABCATLAS_PYTHON`.

For optional analysis objects:

```r
install.packages("SeuratObject")
install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SpatialExperiment"))
```

## Metadata first, download plan second, expression last

```r
library(importABCatlas)
setup_environment()  # selects the previously installed environment
list_datasets()     # bundled routes; no network

cache <- "~/abc_atlas_cache"  # outside the Git repository
release <- "20260711"         # latest release documented when audited

# A preview is for schema inspection, NOT neuronal/non-neuronal benchmarking.
d <- load_data_MERFISH(download_base = cache, release = release,
                       preview_n = 20L)
names(d$cell_metadata)
head(d$gene_data)

ids <- head(d$cell_metadata$cell_label, 5L)
genes <- head(d$gene_data$gene_identifier, 3L)
plan_fetch(data = d, cell_ids = ids, genes = genes)

# Run only after reviewing source files, directory size context and disk space.
x <- fetch_data(data = d, cell_ids = ids, genes = genes,
                return_type = "sce", allow_downloads = TRUE)

# Reuses the local cache; missing files cause an error rather than a download.
x2 <- fetch_data(data = d, cell_ids = ids, genes = genes,
                 return_type = "sce", offline = TRUE)
```

`return_type` accepts `seurat` (default), `sce`, `spatial` (SpatialExperiment), or `matrix` (a list containing a sparse matrix, cell metadata, gene metadata and provenance). For imputed Allen MERFISH use `load_data_MERFISH(imputed=TRUE)` and `data_type="log2"`; it is predicted expression, never counts. Gene row names are source identifiers, not forcibly unique symbols.

**A small output does not imply a small download.** Full source H5AD files are cached before selected rows/columns are read. Metadata can itself be large. The output guard is not a total-memory or network-transfer guarantee. The importer does not claim to eliminate HPC requirements for whole-atlas analyses.

## Discover everything in the selected manifest

`list_files()` lists metadata, expression matrices, image volumes and MapMyCells assets. `get_metadata()` reads arbitrary CSV metadata; `download_atlas_file()` retrieves a declared asset with explicit consent. Image volumes, polygons, ATAC files and electrophysiology are not automatically converted to RNA assays. Resources outside the selected ABC manifest are not implicitly supported.

## Documentation and validation

- `docs/BEGINNER_GUIDE.md`: installation, Git concepts, safe branch workflow, editing, debugging, tests and release steps.
- `docs/DATASET_GUIDE.md`: complete route registry, SEA-AD and MERFISH considerations, schema overrides and non-expression assets.
- `docs/BENCHMARK_PROTOCOL.md`: reproducible cold/warm measurements and one-figure generation without fabricated values.
- `docs/AUDIT_AND_VALIDATION.md`: original defects, workbook interpretation, actual test evidence and remaining gaps.
- `docs/ARCHITECTURE.md`: boundaries and data-flow decisions.
- `docs/RELEASE_CHECKLIST.md`: mandatory release gates, authorship and licensing decisions.
- `docs/SOURCES.md`: primary sources, data release notes and journal policies.

Run Python tests with `python -m pytest tests/python -q`. Run R tests with `devtools::test()` after installing development and optional dependencies. `scripts/validate_live.R` starts in catalogue-only mode; no automatic multi-GB downloads are hidden in the unit tests. The CI workflow is supplied, not already green.

## Maintenance and attribution

Edit code in `R/` and `inst/python/`, and routing in `inst/python/registry.json`. Roxygen comments beside functions are the source of help text. This ZIP includes deliberately labelled bootstrap Rd/NAMESPACE files because R/roxygen2 could not run here; generate canonical documentation once with `Rscript scripts/regenerate_docs.R`, then use `devtools::document()` normally.

`DESCRIPTION` still requires the real maintainer and authors. `LICENSE` explicitly records that the supplied original had no completed code licence. The rights holders must resolve this before public release; no licence was invented. Atlas data retain their own source terms and citation requirements. Do not commit atlas files or clinical data to this repository.
