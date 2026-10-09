# Candidate v8: metadata/plan validation

After installing the unzipped source (folder containing `DESCRIPTION`) and restarting R:

```r
library(importABCatlas)
Sys.setenv(IMPORTABCATLAS_PYTHON = '/Users/IEO7976/.virtualenvs/r-reticulate-env/bin/python3')
source('scripts/test_all_datasets.R')
results <- test_all_datasets(
  download_base = '~/abc_atlas_validation_isolated',
  release = '20260711', preview_n = 20L,
  plan_cells = 5L, plan_genes = 3L,
  cleanup = TRUE,
  output_prefix = 'importABCatlas_dataset_validation_v8')
metadata_runtime <- setNames(results$metadata_elapsed_seconds, results$dataset)
subset(results, load_status != 'PASS' | plan_status != 'PASS')
```

The script never calls `fetch_data`. It uses the static catalogue and
`load_data(dataset)` so every registered route is attempted. Imputed MERFISH
plans log2; other routes default to raw. It deletes only **new** metadata files
created during the current iteration, leaving pre-existing files and expression
files untouched. Manifest files are retained; to measure truly cold runs use a
fresh cache and record that fact in the manuscript.

Selected file sizes are reported only when the SDK manifest contains explicit
per-file bytes or a verified cached file is discoverable. Unknown is not zero.

HMBA PatchSeq may still need a live verified file_map if its expression asset
names differ from the metadata feature_matrix_label. The tool does not guess.
