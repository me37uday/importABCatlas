# Collaborator testing guide

This guide is for collaborators reviewing the native-R implementation of
importABCatlas. The native-R package reads supported ABC Atlas data
directly from R and does not require Python, reticulate, Conda, or a
separately configured Python environment.

## 1. Install the branch

Start a fresh R session and run:

``` r
install.packages("BiocManager")

BiocManager::install(
  c("rhdf5", "HDF5Array"),
  ask = FALSE,
  update = FALSE
)

install.packages("remotes")

remotes::install_github(
  "me37uday/importABCatlas",
  ref = "r-native",
  upgrade = "never"
)

library(importABCatlas)
packageVersion("importABCatlas")
sessionInfo()
```

`rhdf5` and `HDF5Array` are Bioconductor packages used to read ABC Atlas H5AD expression files directly in R. Installing them explicitly prevents dependency-resolution failures on fresh R installations. No Python installation is required.

If installation fails, report the complete installation output before
changing your R environment.

## 2. Inspect available datasets

Run `list_datasets()` and confirm that the expected routes and
expression representations are present.

## 3. Metadata-only smoke test

``` r
abc <- load_data("AgingMouse", n_cells = 5, seed = 1)
abc$cell_metadata
abc$gene_data
```

This should not require downloading a large expression matrix. Check
that cell labels and gene identifiers look intact.

## 4. Plan before fetching expression

``` r
genes <- head(unique(abc$gene_data$gene_identifier), 5)
plan <- plan_fetch(abc, genes = genes)
plan
```

Inspect the selected representation, required files, file sizes, cache
state, and additional download required. A small cell request can still
require a large source H5AD file.

## 5. Matrix output

When the required file size is acceptable:

``` r
x <- fetch_data(abc, genes = genes, return_type = "matrix")
dim(x$expression)
head(x$cell_metadata)
head(x$gene_data)
```

Expression rows and columns should correspond to the returned gene and
cell metadata.

## 6. Optional output formats

The package also supports `return_type = "seurat"`, `"sce"`, and
`"spatial"` when the corresponding optional R packages are installed.
For SpatialExperiment, use a spatial route such as `Zhuang4` and confirm
that spatial coordinates are present.

## 7. Raw and transformed expression

Use `list_datasets()` to determine which representations are available.
Do not assume every route has both `raw` and `log2`.

The package preserves the upstream representation and does not convert
transformed or imputed expression into raw counts. Some expression
matrices contain only a subset of the genes in broader metadata;
unavailable requested genes should produce a clear error rather than
being silently dropped.

## 8. Cache and offline behavior

Inspect the active cache with `abc_download_dir()`. After required files
are cached, test `offline = TRUE`. Use offline mode only when required
files are already available locally.

## 9. What to report

For successful tests report:

-   operating system;
-   R version;
-   importABCatlas version / commit;
-   dataset;
-   data type;
-   return type;
-   result.

For failures additionally report:

-   exact R command;
-   complete error/warning;
-   whether the required file was cached;
-   whether a partial download existed;
-   whether offline mode was enabled;
-   expected behavior;
-   observed behavior;
-   `sessionInfo()`.

Do not include passwords, access tokens, credentials, or sensitive data.

## 10. Scope

Collaborator testing is intended to detect installation, portability,
usability, routing, cache, and object-construction problems on
independent systems. It complements, but does not replace, automated
tests, cross-platform CI, full real-data route acceptance, offline
zero-download acceptance, and feature-universe auditing.
