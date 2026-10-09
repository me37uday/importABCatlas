# Dataset scope, scientific handling and extension

## Scope and evidence

The current official ABC documentation covers more distinct expression components than the ten high-level rows in the supplied workbook. The 23 entries below are importer routes, not 23 independent biological studies. They are configured for the documented release 20260711. Every route passed a synthetic expected-schema fixture; none has a completed live test in this delivery environment. Exact live filenames, schema details and optional object conversion remain release gates. [S1]

| Route | Primary metadata directory | Special treatment |
|---|---|---|
| `WMB` | `WMB-10X` | Routes v2/v3/Multi expression shards from WMB metadata. |
| `ConsensusMouseAIBS` | `Consensus-WMB-AIBS-10X` | Consensus taxonomy; reuses WMB expression and gene source. |
| `ConsensusMouseMacosko` | `Consensus-WMB-Macosko-10X` | Separate expression files and native feature universe. |
| `AgingMouse` | `Zeng-Aging-Mouse-10Xv3` | WMB genes; separate native and WMB mapping annotations. |
| `DevelopingMouse` | `Developing-Mouse-Vis-Cortex-10X` | Visual cortex developmental taxonomy; preserve actual donor age. |
| `MERFISH` | `MERFISH-C57BL6J-638850` | Measured Allen panel; physical coordinates. |
| `MERFISH_imputed` | `MERFISH-C57BL6J-638850` | Predicted log2 only; actual H5AD var defines panel. |
| `Zhuang1` | `Zhuang-ABCA-1` | Native measured spatial panel and coordinates. |
| `Zhuang2` | `Zhuang-ABCA-2` | Native measured spatial panel and coordinates. |
| `Zhuang3` | `Zhuang-ABCA-3` | Native measured spatial panel and coordinates. |
| `Zhuang4` | `Zhuang-ABCA-4` | Native measured spatial panel and coordinates. |
| `WHB` | `WHB-10Xv3` | Native human feature universe and taxonomy. |
| `PMDBS` | `ASAP-PMDBS-10X` | Native PMDBS genes; optional MapMyCells mapping is separate. |
| `HMBA` | `HMBA-10XMultiome-BG-Aligned` | Aligned cross-species feature universe. |
| `HMBA_Human` | `HMBA-10XMultiome-BG` | Species-specific files require live discovery or verified overrides. |
| `HMBA_Macaque` | `HMBA-10XMultiome-BG` | Species-specific files require live discovery or verified overrides. |
| `HMBA_Marmoset` | `HMBA-10XMultiome-BG` | Species-specific files require live discovery or verified overrides. |
| `HMBA_spatial_Human` | `HMBA-MERSCOPE-H22.30.001-BG` | Native measured spatial panel and coordinates. |
| `HMBA_spatial_Macaque` | `HMBA-MERSCOPE-QM23.50.001-BG` | Native measured spatial panel and coordinates. |
| `HMBA_spatial_Marmoset` | `HMBA-Xenium-CJ23.56.004-BG` | Native measured spatial panel and coordinates. |
| `HMBA_PatchSeq` | `HMBA-Macaque-PatchSeq` | RNA subset; NCBI-style gene IDs; non-neuronal cohort N/A. |
| `SEAAD_Multiregion` | `SEA-AD-Multiregion-10X` | Regional expression shards, donor/disease joins, own taxonomy. |
| `SEAAD_CaH` | `SEA-AD-CaH-10X` | Separate Caudate release and taxonomy; no implicit merge. |


The broad generic asset API is separate from these expression routes. `list_files()` enumerates the selected manifest; `get_metadata()` reads named metadata CSVs; `download_atlas_file()` retrieves named assets such as volumes, geometry and MapMyCells models. The importer does not transform every downloadable file into Seurat. It does not claim access to controlled/external SEA-AD resources, to an arbitrary future manifest schema, or automatic ATAC/electrophysiology integration.

## SEA-AD: why two dedicated routes are necessary

The multiregion source uses `SEA-AD-Multiregion-10X` and `SEA-AD-Multiregion-taxonomy`; CaH uses its own corresponding directories. Both obtain genes from their own expression metadata, not WHB merely because the species is human. The multiregion release includes nuclei beyond the quality-filtered reference used to define the taxonomy; retaining source cells is not equivalent to asserting that all are suitable for analysis. [S6-S8]

Cell identity has two documented forms: `cell_label` combines a barcode with a sample identity, while `exp_component_name` can index expression hosted in another layout. A barcode alone is not a global cell identifier. The reader first tries exact full cell labels against H5AD observation names, then the complete alternate identifiers; it fails when neither gives an exact match. It does not trim suffixes or guess sample joins. An explicit `obs_column` specifies a metadata column to match against the H5AD index, not an instruction to round or rewrite IDs.

Expression is separated into regional files, selected using `feature_matrix_label`. The multiregion example lists region/file prefixes such as AnG-10X, DFC-10X and HIP-10X. Requesting a single region reduces the number of H5AD files required, not necessarily the size of a single file. The original metadata still has to be downloaded. Reference-only annotations, low-quality cells, missing mappings and intentionally excluded cells must be distinguished in the analysis rather than converted into a default biological class. [S6-S7]

Donor, library and disease tables represent different entities. Joins use source keys and require a unique right-hand record; colliding columns are retained with suffixes. This avoids multiplying a cell through an ambiguous disease or library join. Clinical column names with spaces are kept unchanged. Missing disease fields remain missing, not 'control'. The supplied loader makes no clinical interpretation and performs no case-control analysis.

For downstream study design, retain donor and region identifiers. Thousands of nuclei from one donor do not constitute thousands of independent donors. Check donor overlap before treating CaH and multiregion selections as independent cohorts. Distinguish technical sampling for a runtime benchmark from a donor-aware biological comparison. These are analysis-design considerations, not statistical features supplied by the importer.

A first schema check is:

```r
s <- load_data_SEAAD_Multiregion(download_base = cache,
                                release = "20260711", preview_n = 100L)
names(s$cell_metadata)
unique(s$cell_metadata$feature_matrix_label)
head(s$cell_metadata[c("cell_label", "exp_component_name")])
```

After inspecting the actual metadata, a multiregion cohort can use the documented `Class` values `Neuronal: GABAergic`, `Neuronal: Glutamatergic`, or `Non-neuronal and Non-neural`. Verify the values and capitalisation in the pinned release; do not transfer them to other taxonomies by assumption. This code is a live workflow supplied for execution, not a result run here. [S6]

```r
neurons <- load_data_SEAAD_Multiregion(
  download_base = cache, release = "20260711",
  filters = list(feature_matrix_label = "AnG-10X",
                 Class = c("Neuronal: GABAergic", "Neuronal: Glutamatergic")),
  n_cells = 100L, seed = 1729L)
plan_fetch(data = neurons,
           genes = head(neurons$gene_data$gene_identifier, 3L))
```

No unreviewed cell is silently reassigned to make that filter yield 100 cells. Broader selections may require considerable metadata RAM because full annotation precedes `n_cells` sampling. Save reviewed exact IDs and use `cell_ids` for subsequent runs; the main cell CSV is then filtered in chunks. Supplemental metadata and indices can still consume memory.

The access path is RNA expression even where a library was generated by multiome profiling. A corresponding ATAC matrix, fragment file, spatial assay or neuropathology image is not implied to be downloaded or combined. UMAP x/y coordinates in SEA-AD are visual embeddings, so this registry is not marked as spatial. [S6-S8]

## MERFISH: measured and imputed are different products

Measured Allen MERFISH and Zhuang1-4 use their own measured feature panels and cell metadata. The imputed Allen route uses the same cells but a different expression directory. It has predicted log2 expression for a broader marker panel, not additional independently measured transcripts. The source notebook identifies 8,460 imputed markers and a large dense H5AD source; these are upstream descriptions, not package benchmark measurements. [S3-S4]

```r
m <- load_data_MERFISH(source = "Allen", download_base = cache,
                       release = "20260711", preview_n = 20L)
z <- load_data_MERFISH(source = "Zhuang4", download_base = cache,
                       release = "20260711", preview_n = 20L)
i <- load_data_MERFISH(imputed = TRUE, download_base = cache,
                       release = "20260711", preview_n = 20L)
plan_fetch(data = i, genes = c("Calb2", "Baiap3", "Lypd1"),
           data_type = "log2")
```

A gene appearing in the WMB candidate annotation is not proof that it appears in the imputed matrix. The imputed reader validates the real H5AD feature index; when all genes are requested it restricts to that authoritative universe after the authorised download. A conservative pre-download output guard can use the larger candidate set. Explicit small gene selections are the simplest way to test this route.

Cell identifiers remain strings; some Zhuang identifiers exceed ordinary integer precision by many digits. Coordinates retain source units and section/specimen fields. A two-dimensional plane is not automatically a common three-dimensional reference. `return_type="spatial"` uses native slab coordinates when supplied by HMBA, otherwise declared anatomical x/y for a spatial route, or explicitly provided verified coordinate columns. Additional CCF coordinates are a separate annotation, not a replacement for native positions. [S4, S9]

```r
# Inspect names first; only use these columns when the source defines them.
xy <- fetch_data(data = m, cell_ids = head(m$cell_metadata$cell_label, 5L),
                 genes = head(m$gene_data$gene_identifier, 3L),
                 return_type = "spatial", allow_downloads = TRUE)
```

Do not call an imputed output a counts assay or use it as independent evidence that a predicted gene was measured. A merged measured/imputed object would need separate assays and explicit interpretation; this importer intentionally returns them separately.

## HMBA: aligned genes, native species, spatial assays and Patch-seq

The aligned route and native-species routes have different feature universes. Cross-species alignment is not equivalent to renaming all genes to human symbols. Native Human, Macaque and Marmoset file discovery is constrained to the actual manifest and requires unique metadata matches. The descriptive source did not provide an independently downloadable complete native file listing in this execution environment. This is a known validation gap, not a hidden claim of success. [S1, S10]

```r
h <- load_data_HMBA(download_base = cache, preview_n = 20L)
hn <- load_data_HMBA_species("Human", download_base = cache,
                              preview_n = 20L)
hs <- load_data_HMBA_spatial("Marmoset", download_base = cache,
                              preview_n = 20L)
```

Official notebooks and release notes use differing x/X capitalisation in some directory names. Resolution accepts a unique manifest-backed case match, not an invented path. HMBA also has an upstream normalization correction dated 20260331; pinning a recent reviewed manifest is scientifically meaningful, not merely convenient. [S1]

MERSCOPE Human/Macaque and Xenium Marmoset keep their own panels and slab annotations. Patch-seq preserves NCBI-style source gene identifiers rather than assuming every atlas uses Ensembl. Its electrophysiology sweeps and other one-to-many measurements should be accessed separately, not left-joined into one row per cell without a defined aggregation. A non-neuronal Patch-seq benchmark is not applicable to the documented neuronal collection. [S9, S11]

## PMDBS, consensus and aging mappings

PMDBS uses its own gene table. Its MapMyCells output maps cells to another reference; a mapped class is not a new native taxonomy and has mapping confidence semantics. Read the actual mapping file and select a scientifically justified confidence rule, keeping unmapped cells explicit. MapMyCells CSVs may start with provenance comments; the reader skips only the leading comment block, preserving hex colour strings elsewhere. [S5]

```r
p <- load_data_PMDBS(download_base = cache, preview_n = 100L)
a <- get_metadata("ASAP-PMDBS-taxonomy", "mmc_results_siletti_whb",
                  download_base = cache, release = "20260711")
p <- join_metadata(p, a, by = "cell_label", prefix = "whb_map_")
```

Inspect `list_files("ASAP-PMDBS-taxonomy", ...)` before using an assumed mapping name on another release. `join_metadata` rejects non-unique keys and a known manifest mismatch. A manually constructed annotation without a manifest attribute remains the caller's provenance responsibility.

Consensus AIBS reuses WMB expression but applies a consensus taxonomy; Macosko has a distinct feature universe. Do not pad an unmeasured gene with zero to create a merged matrix. Aging Mouse retains its source cluster annotations and the separate WMB correspondence rather than replacing source labels with row numbers. Developmental ages are preserved as provided; avoid recoding one postnatal day to another merely to imitate a plotting notebook. [S1, S12]

## Verified overrides and new manifest assets

When a resolver fails, query the actual selected manifest and inspect the relevant metadata before editing the registry. For native HMBA, for example:

```r
f <- list_files("HMBA-10XMultiome-BG", download_base = cache,
                release = "20260711")
f[f$kind == "metadata", ]
# Replace the two values below with exact names that you have verified above.
h <- load_data_HMBA_species("Human", download_base = cache,
      overrides = list(cell_file = "EXACT_VERIFIED_CELL_FILE",
                       gene_file = "EXACT_VERIFIED_GENE_FILE"))
```

Placeholders are intentionally not guessed real filenames. Supported override keys are directory, cell_file, gene_directory, gene_file, taxonomy, expression_directories, spatial, file_selector and imputed. A file_map passed to plan_fetch/fetch_data maps an exact feature_matrix_label to a named list containing directory and file_name, both validated against the manifest. Extend fixtures and execute a live source comparison before making a temporary override the default registry rule.

For non-expression files, discover the exact names and use `download_atlas_file(directory, file_name, allow_downloads=TRUE)`. Image/geometry decoding and modality-specific analysis remain separate. All sources are listed in `SOURCES.md`.
