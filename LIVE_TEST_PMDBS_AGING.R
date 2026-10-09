library(importABCatlas)

# Use the virtualenv path without resolving its symlink.
good_python <- "/Users/IEO7976/.virtualenvs/r-reticulate-env/bin/python3"
Sys.setenv(IMPORTABCATLAS_PYTHON = good_python)
options(importABCatlas.python = good_python)
stopifnot(identical(importABCatlas:::.abc_python(), good_python))

# Live catalogue should tolerate empty optional Allen file categories.
print(head(list_datasets(live = TRUE), 20))

# PMDBS: native cells + own human gene table. Reference mappings are optional.
pmdbs <- load_data_PMDBS(download_base = "~/abc_atlas_cache", release = "20260711", preview_n = 20L)
stopifnot(nrow(pmdbs$cell_metadata) == 20L)
stopifnot(all(c("cell_label", "cluster_label") %in% names(pmdbs$cell_metadata)))
stopifnot(all(c("gene_identifier", "gene_symbol", "molecular_type", "description") %in% names(pmdbs$gene_data)))
print(dim(pmdbs$cell_metadata)); print(dim(pmdbs$gene_data)); print(head(pmdbs$cell_metadata)); print(head(pmdbs$gene_data))

# Optional PMDBS reference mappings (large files; these may already be cached).
pmdbs_seaad <- get_metadata("ASAP-PMDBS-taxonomy", "mmc_results_seaad", download_base = "~/abc_atlas_cache", release = "20260711", preview_n = 3L)
pmdbs_whb <- get_metadata("ASAP-PMDBS-taxonomy", "mmc_results_siletti_whb", download_base = "~/abc_atlas_cache", release = "20260711", preview_n = 3L)
stopifnot("cell_label" %in% names(pmdbs_seaad), "cell_label" %in% names(pmdbs_whb))

# Aging Mouse: native annotations plus WMB cross-mapping.
aging <- load_data_AgingMouse(download_base = "~/abc_atlas_cache", release = "20260711", preview_n = 20L)
stopifnot(nrow(aging$cell_metadata) == 20L)
stopifnot(all(c("cell_label", "cluster_alias", "cluster_name", "wmb_cluster_alias", "wmb_class_name") %in% names(aging$cell_metadata)))
stopifnot(all(c("gene_identifier", "gene_symbol") %in% names(aging$gene_data)))
print(dim(aging$cell_metadata)); print(dim(aging$gene_data)); print(head(aging$cell_metadata)); print(head(aging$gene_data))

# Plan only: do not download expression matrices yet.
pids <- head(pmdbs$cell_metadata$cell_label, 5L); pgenes <- head(pmdbs$gene_data$gene_identifier, 3L)
aids <- head(aging$cell_metadata$cell_label, 5L); agenes <- head(aging$gene_data$gene_identifier, 3L)
print(plan_fetch(pmdbs, cell_ids = pids, genes = pgenes, data_type = "raw"))
print(plan_fetch(aging, cell_ids = aids, genes = agenes, data_type = "raw"))
