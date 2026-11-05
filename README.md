# importABCatlas

`importABCatlas` provides an R interface to the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).

It lets you explore available datasets, select cells and genes, check how much data needs to be downloaded, and load expression data directly into R.

**No Python installation is required.**

## Installation

Install the package from the `r-native` branch:

```r
install.packages("remotes")

remotes::install_github(
  "me37uday/importABCatlas",
  ref = "r-native"
)
```

Then load it:

```r
library(importABCatlas)
```

## Quick start

The basic workflow is: choose a dataset, select cells, check the download, and fetch expression.

### 1. See the available datasets

```r
list_datasets()
```

This shows the supported datasets, including which are spatial and which expression types are available.

### 2. Choose cells

```r
abc <- load_data(
  "AgingMouse",
  filters = list(anatomical_division_label = "STR"),
  n_cells = 500,
  seed = 1
)
```

`load_data()` loads metadata and selects cells. It does **not** download the large expression matrices.

For exact cells:

```r
abc <- load_data(
  "AgingMouse",
  cell_ids = my_cell_labels
)
```

### 3. Check the download

```r
genes <- c("P2ry12", "Cx3cr1", "Tmem119")

plan_fetch(
  abc,
  genes = genes
)
```

`plan_fetch()` shows which expression files are needed, their sizes, whether they are already cached, and how much additional data must be downloaded. It does **not** start the download.

### 4. Fetch expression

```r
x <- fetch_data(
  abc,
  genes = genes,
  return_type = "matrix"
)
```

Required expression files are downloaded if necessary and reused from the local cache on future requests.

## Choosing an output format

### Sparse matrix

```r
x <- fetch_data(abc, genes = genes, return_type = "matrix")
```

This is the simplest option and does not require an additional analysis framework.

### Seurat

```r
x <- fetch_data(abc, genes = genes, return_type = "seurat")
```

Requires `SeuratObject`.

### SingleCellExperiment

```r
x <- fetch_data(abc, genes = genes, return_type = "sce")
```

Requires `SingleCellExperiment`.

### SpatialExperiment

```r
spatial <- load_data(
  "Zhuang4",
  n_cells = 500,
  seed = 1
)

x <- fetch_data(
  spatial,
  genes = genes,
  return_type = "spatial"
)
```

Requires `SpatialExperiment`. Spatial coordinates provided by the dataset are included in the returned object.

## Raw and transformed expression

Use `list_datasets()` to see which representations are available.

- `"raw"` — raw expression with count semantics.
- `"log2"` — expression transformed upstream by the data provider.

If `data_type` is not specified, the package uses the dataset's default representation.

```r
x <- fetch_data(
  abc,
  genes = genes,
  data_type = "raw",
  return_type = "matrix"
)
```

The package preserves the representation supplied by the ABC Atlas and does not silently convert transformed expression into counts.

Some representations may contain expression for only a subset of the genes in the broader dataset metadata. If a requested gene is not available in that expression matrix, `fetch_data()` reports it rather than silently dropping it.

## Downloads and local cache

Downloaded files are cached so they do not need to be downloaded again.

Set the cache location with:

```r
options(importABCatlas.download_base = "/path/to/abc_atlas")
```

Check the active location:

```r
abc_download_dir()
```

If required files are already cached, network access can be disabled:

```r
abc <- load_data("AgingMouse", offline = TRUE)

x <- fetch_data(
  abc,
  genes = genes,
  offline = TRUE
)
```

## Complete example

```r
library(importABCatlas)

list_datasets()

abc <- load_data(
  "AgingMouse",
  filters = list(anatomical_division_label = "STR"),
  n_cells = 500,
  seed = 1
)

genes <- c("P2ry12", "Cx3cr1", "Tmem119")

plan_fetch(abc, genes = genes)

x <- fetch_data(
  abc,
  genes = genes,
  return_type = "matrix"
)

x
```

## Getting help

```r
?list_datasets
?load_data
?plan_fetch
?fetch_data
```

For a more detailed walkthrough:

```r
vignette("importABCatlas", package = "importABCatlas")
```
