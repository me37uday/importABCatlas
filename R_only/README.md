# importABCatlas (R-only)

A version of `importABCatlas` that runs entirely in R. It has the same functions
and returns the same results as the reticulate-based package in the repository
root, without Python.

It lives in `R_only/` while it is being developed. It is a complete package with
its own `DESCRIPTION`, `R/`, `man/` and `tests/`, and the root package's build
ignores it.

## Usage

```r
devtools::load_all("R_only")        # or: devtools::install("R_only")

options(importABCatlas.download_base = "/data/abc_atlas")   # optional

WHB <- load_data("WHB")             # also "WMB", "AgingMouse", "HMBA", "PMDBS"
names(WHB$unique_values)

obj <- fetch_data(
  metadata  = WHB$cell_metadata,
  filters   = list(feature_matrix_label = "WHB-10Xv3-Nonneurons",
                   region_of_interest_label = "Human CEN"),
  gene_data = WHB$gene_data
)
```

`abc_cache_connect()` returns a handle for direct file access:

```r
cache <- abc_cache_connect()
cache$get_file_path("WHB-10Xv3", "WHB-10Xv3-Neurons/raw")       # downloads if needed
cache$get_metadata_dataframe("ASAP-PMDBS-10X", "donor")
```

For a full walkthrough (choosing cells, fetching their expression, working
offline), see the vignette: `vignette("importABCatlas", package = "importABCatlas")`
after installing with `build_vignettes = TRUE`, or read
[vignettes/importABCatlas.Rmd](vignettes/importABCatlas.Rmd). The vignette is
precomputed because its code downloads several GB. To regenerate it after a
change, edit `importABCatlas.Rmd.orig` and run `vignettes/precompile.R`.

## How it works

| Step | Python package (`abc_atlas_access`) | R-only |
|---|---|---|
| List releases, pick a manifest | boto3 | S3 listing over HTTPS (`curl`, `xml2`) |
| Download and verify files | boto3 + md5 | `curl::multi_download` (resumable `.part` file) + `tools::md5sum` |
| Read metadata CSVs | `pandas.read_csv` | `data.table::fread` with pandas' NA strings |
| Read expression (h5ad) | `anndata` (backed) | `HDF5Array::H5SparseMatrix`: only the selected cells and genes are read, in blocks, and nothing is densified |

The download directory layout is the same as the Python package's
(`releases/<version>/manifest.json`, `_manifest_last_used.txt`,
`_downloaded_data.json`, and files at their manifest `relative_path`). A
directory filled by either package works with the other one, and nothing is
downloaded twice.

## Differences from the reticulate version

- `setup_environment()` is gone because there is no Python environment to
  set up.
- Metadata columns get clean R types. Integer columns are `integer` instead of
  `double`, missing values are `NA` instead of `NaN`, and columns that pandas
  read with mixed types (which reticulate turned into list columns, e.g.
  `alignment_job_id`) are plain vectors. The values are identical.
- Genes are matched to the expression file by `gene_identifier`, not by
  position. With a gene table from another dataset, `fetch_data()` stops with a
  clear error instead of failing inside Python.
- Interrupted downloads resume instead of starting over.

## Tests

```r
devtools::test("R_only")
```

The offline tests build a small fixture cache (manifest, CSVs and h5ad files).
The online tests only list releases and download one small file. They are
skipped on CRAN and when the network is unavailable.
