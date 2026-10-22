# Independent collaborator testing — 0.2.0.9007

Use a fresh R session and a separate atlas cache. Do not install over an
`importABCatlas` namespace that is currently loaded.

## Install from a GitHub branch

```r
install.packages("remotes")
remotes::install_github(
  "me37uday/importABCatlas@<branch-name>",
  dependencies = NA,
  upgrade = "never",
  force = TRUE
)
```

Restart R, then run the one-time automatic setup:

```r
library(importABCatlas)
setup_importABCatlas()
diagnose()
```

No Python path should be configured. The package creates and reuses its own
private Python 3.11 environment.

## Minimum smoke test

```r
packageVersion("importABCatlas")
list_datasets()

example <- example_data("matrix")
stopifnot(identical(dim(example$expression), c(4L, 3L)))
```

## Recommended live test

Choose one route, pin release `20260711`, and use a new cache:

```r
d <- load_data(
  dataset = "Zhuang4",
  download_base = "~/abc_atlas_collaborator_cache",
  release = "20260711",
  preview_n = 1000L
)

cells <- head(d$cell_metadata$cell_label, 3L)

p <- plan_fetch(data = d, cell_ids = cells, genes = NULL)
print(p)

x <- fetch_data(
  data = d,
  cell_ids = cells,
  genes = NULL,
  return_type = "matrix",
  allow_downloads = TRUE
)

stopifnot(identical(colnames(x$expression), cells))
```

After this fetch, use identifiers from `x$gene_data$gene_identifier` for
explicit-gene tests.

## Focused feature-axis regression

Run `TEST_FAILED_DATASETS_9007.R` for `HMBA_Macaque` and
`MERFISH_imputed`. The script uses the managed Python environment and does not
need a Python argument. It can reuse an existing cache and will not download
absent large files unless `download_missing = TRUE` is set deliberately.

## What to report

Open one GitHub issue per distinct problem and include:

- package version and Git commit;
- operating system and architecture;
- R version and `.libPaths()`;
- output of `diagnose()`;
- Allen release;
- dataset route and raw/log2 representation;
- exact command;
- complete error message and retained diagnostic directory;
- whether the source file was cached, partial, or newly downloaded;
- free disk space before the operation;
- whether the error is reproducible in a fresh R session.

Do not upload Allen expression matrices or sensitive donor/cell-level data to
GitHub issues.
