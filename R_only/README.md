# importABCatlas

`importABCatlas` is an R-native, release-aware interface to the
[Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).

It discovers supported dataset routes, loads and joins cell and gene metadata,
selects exact or sampled cells, plans expression downloads before transfer, and
reads selected H5AD expression directly in R.

No Python installation is required.

## Discover datasets

    library(importABCatlas)
    list_datasets()

`list_datasets()` reports the curated routes supported by the package together
with their modality, spatial and imputed status, and available expression
representations. It performs no network access.

## Load metadata and select cells

`load_data()` assembles cell and gene metadata. It does not download the large
expression matrices.

    abc <- load_data(
      "AgingMouse",
      filters = list(
        anatomical_division_label = "STR"
      ),
      n_cells = 500,
      seed = 1
    )

The result is an `abc_data` object containing cell metadata, gene metadata,
useful distinct metadata values, and dataset/cache/release provenance.

For an exact set of cells, use `cell_ids`. Requested order is preserved:

    abc <- load_data(
      "AgingMouse",
      cell_ids = my_cell_labels
    )

## Plan before downloading expression

ABC Atlas expression files can be large. Inspect the request before transfer:

    genes <- c("P2ry12", "Cx3cr1", "Tmem119")

    plan <- plan_fetch(
      abc,
      genes = genes
    )

    plan

`plan_fetch()` reports the selected representation, requested dimensions,
expression files, cache state, resumable partial downloads and additional
download required. It does not download expression matrices.

## Fetch expression

Pass the `abc_data` object directly:

    x <- fetch_data(
      abc,
      genes = genes,
      return_type = "matrix"
    )

Four output representations are supported:

    x_matrix <- fetch_data(abc, genes = genes, return_type = "matrix")
    x_seurat <- fetch_data(abc, genes = genes, return_type = "seurat")
    x_sce <- fetch_data(abc, genes = genes, return_type = "sce")
    x_spatial <- fetch_data(spatial_abc, genes = genes, return_type = "spatial")

`return_type = "matrix"` returns an `abc_matrix` bundle containing sparse
expression, selected cell metadata, selected gene metadata, representation and
provenance.

Seurat, SingleCellExperiment and SpatialExperiment outputs require their
corresponding optional packages.

## Expression semantics

Dataset routes declare the representations actually published upstream.

- `data_type = "raw"` retains count semantics.
- `data_type = "log2"` remains upstream transformed expression and is never
  presented as counts.
- `MERFISH_imputed` represents upstream published imputed spatial expression.
  `importABCatlas` does not perform the imputation.

Use `list_datasets()` to inspect each route's default and available
representations.

## Gene identity

Expression rows use canonical `gene_identifier` values. `gene_symbol` remains
available as annotation and for gene lookup.

Gene symbols can be supplied to `genes` when they resolve uniquely. If a symbol
is ambiguous, the package asks for a canonical identifier rather than silently
choosing one.

## Spatial data

Spatial routes are identified by `list_datasets()`.

For `return_type = "spatial"`, verified spatial coordinates are carried into a
SpatialExperiment. Automatic coordinate selection is restricted to routes
declared spatial. Otherwise, verified `coordinate_columns` must be supplied
explicitly.

This prevents dimensionality-reduction coordinates from being mistaken for
physical spatial positions.

## Cache and releases

    options(importABCatlas.download_base = "/data/abc_atlas")
    abc_download_dir()

Cache resolution order is:

1. explicit `download_base`;
2. `options(importABCatlas.download_base = ...)`;
3. `ABC_ATLAS_CACHE`;
4. the package's per-user cache directory.

Downloads are cached and resumable.

A release can be pinned explicitly for reproducibility:

    abc <- load_data(
      "AgingMouse",
      release = "20260711"
    )

After required files are cached, `offline = TRUE` prohibits network access.

## Direct cache access

`abc_cache_connect()` provides lower-level access to files and metadata in the
release manifest:

    cache <- abc_cache_connect()

    donors <- cache$get_metadata_dataframe(
      directory = "Zeng-Aging-Mouse-10Xv3",
      file_name = "donor"
    )

## Development

The native package currently lives in `R_only/` in this repository.

    devtools::load_all("R_only")
    devtools::test("R_only")

For a complete walkthrough:

    vignette("importABCatlas", package = "importABCatlas")
