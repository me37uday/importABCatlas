# importABCatlas

`importABCatlas` lets you work with Allen Brain Cell Atlas data from R without
manually downloading files or configuring Python.

The usual workflow is:

1. choose a dataset;
2. load its cell and gene information;
3. check which expression file is needed and how large it is;
4. download only after you confirm;
5. receive an R matrix, Seurat object, SingleCellExperiment, or SpatialExperiment.

## Install

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "me37uday/importABCatlas",
  dependencies = NA,
  upgrade = "never"
)
```

Restart R or RStudio, then run these three lines:

```r
library(importABCatlas)
setup_importABCatlas()
diagnose()
```

That is the complete software setup. The package creates and maintains its own
isolated Python environment. You do **not** need to find a Python path, install
Conda, edit `.Renviron`, or use Homebrew. First-time setup requires internet
access and can take a few minutes. Later sessions reuse the same environment.

Older development versions asked users to set Python environment variables.
Normal managed mode now ignores those old paths. Calling
`setup_importABCatlas()` also resets any stale package-specific external mode in
the current R session.

If you skip `setup_importABCatlas()`, the package runs it automatically the first
time Python is needed.

Check the small built-in example:

```r
example <- example_data("matrix")
dim(example$expression)
```

The result should be 4 genes by 3 cells.

## A small first example

The example below uses `Zhuang4`, one of the smaller spatial datasets.

```r
library(importABCatlas)

# 1. Load metadata. Expression is not downloaded yet.
d <- load_data(
  dataset = "Zhuang4",
  download_base = "~/abc_atlas_cache",
  release = "20260711",
  preview_n = 1000L
)

# 2. Choose a few cells.
cells <- head(d$cell_metadata$cell_label, 5L)

# 3. See the exact file, its size, cache status, and required disk space.
plan <- plan_fetch(
  data = d,
  cell_ids = cells,
  genes = NULL,
  data_type = "raw"
)
print(plan)

# 4. Download and return an R sparse-matrix result.
x <- fetch_data(
  data = d,
  cell_ids = cells,
  genes = NULL,
  data_type = "raw",
  return_type = "matrix",
  allow_downloads = TRUE
)

x$expression
x$cell_metadata
x$gene_data
```

`genes = NULL` is the safest first query. It asks the H5AD file which genes are
actually measured. You can then make a smaller strict query:

```r
measured_genes <- head(x$gene_data$gene_identifier, 20L)

small <- fetch_data(
  data = d,
  cell_ids = cells,
  genes = measured_genes,
  data_type = "raw",
  return_type = "matrix",
  allow_downloads = FALSE,
  offline = TRUE
)
```

## The three main functions

### `load_data()`

Loads cell metadata, gene metadata, donor information, taxonomy, and spatial
coordinates where available. It does not load expression.

```r
d <- load_data("MERFISH", preview_n = 1000L)
```

See all supported routes:

```r
list_datasets()
```

Dataset-specific helpers are also available, for example:

```r
load_data_WMB()
load_data_WHB()
load_data_MERFISH()
load_data_PMDBS()
load_data_AgingMouse()
load_data_SEAAD_Multiregion()
```

### `plan_fetch()`

Shows exactly what would be downloaded without downloading expression. It
reports:

- the selected H5AD file or files;
- exact file size, not folder size;
- raw or log2 representation;
- cached and resumable bytes;
- additional transfer required;
- free disk space and whether it is sufficient.

Always inspect the plan before a large download.

### `fetch_data()`

Downloads or reuses the planned H5AD file and returns only the requested cells
and genes. Downloads show progress, transfer speed, and ETA. Interrupted
transfers are saved as `.part` files and resume when the same cache is used.

Expression downloads require explicit consent:

```r
allow_downloads = TRUE
```

## R output types

Sparse matrix plus metadata:

```r
return_type = "matrix"
```

SingleCellExperiment:

```r
return_type = "sce"
```

Seurat:

```r
return_type = "seurat"
```

SpatialExperiment, for routes with anatomical coordinates:

```r
return_type = "spatial"
```

Optional object types require their corresponding R packages. Matrix output has
the fewest dependencies and is the best first test.

## Raw and log2 data

Most routes offer both `raw` and `log2`:

```r
data_type = "raw"
data_type = "log2"
```

See what each route publishes:

```r
list_datasets()[, c(
  "dataset",
  "available_data_types",
  "default_data_type",
  "validation"
)]
```

`MERFISH_imputed` and `HMBA_PatchSeq` are log2-only in release `20260711`.
The package rejects unavailable representations before downloading.

## Large files and disk safety

A request for only a few cells can still require a complete multi-gigabyte H5AD
file. `plan_fetch()` therefore checks the exact source object before transfer.
By default, importABCatlas:

- blocks downloads when the exact size is unknown;
- checks free disk space;
- keeps a 1 GiB reserve;
- adds a 5% safety margin;
- retries temporary network failures;
- resumes partial downloads.

The cache is stored under the `download_base` that you choose, normally:

```text
~/abc_atlas_cache
```

Delete that directory only when you intentionally want to remove all downloaded
atlas files.

## If something goes wrong

Start with:

```r
diagnose()
```

Repair the private Python environment without choosing a path:

```r
setup_importABCatlas(repair = TRUE)
```

Then restart R and retry the same command. Keep the diagnostic directory shown
in any backend error message when reporting a problem.

For a bug report, include:

- the complete R command;
- the complete error message;
- `packageVersion("importABCatlas")`;
- `diagnose()`;
- operating system;
- dataset and Allen release;
- whether the file was new, cached, or partially downloaded.

## Advanced: use your own Python

Most users should not do this. To deliberately use an existing environment:

```r
setup_environment(python = "/absolute/path/to/python")
```

The environment is accepted only when all required Python packages are present.
Run `setup_importABCatlas()` to return to automatic managed mode.

## Testing and development

Contributor instructions are in:

- `docs/COLLABORATOR_TESTING.md`
- `docs/GITHUB_AND_COLLABORATOR_STEPS_9007.md`
- `docs/BEGINNER_GUIDE.md`
- `docs/DATASET_GUIDE.md`
- `docs/BENCHMARK_PROTOCOL.md`

Run package tests before opening a pull request:

```r
devtools::test()
devtools::document()
```

```bash
python -m pytest -q tests/python
R CMD build .
R CMD check --no-manual importABCatlas_0.2.0.9007.tar.gz
```

## Data source

Allen Brain Cell Atlas access documentation:
https://alleninstitute.github.io/abc_atlas_access/

The automatic Python environment uses the managed-environment mechanism
recommended by `reticulate`:
https://rstudio.github.io/reticulate/reference/py_require.html

Allen datasets retain their own terms and licences. The package repository must
also contain the final maintainer details and chosen software licence before a
formal public release.
