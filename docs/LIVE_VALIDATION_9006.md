# Live validation status for candidate 0.2.0.9006

Pinned Allen release: `20260711`.

## Evidence established before this candidate

- 23/23 bundled routes completed live metadata loading.
- 23/23 routes completed exact selected-file planning with object-level sizes.
- 21/23 routes completed live simple and complex expression retrieval, warm-cache equality, axis/order validation, and available object-construction checks.
- The complete HMBA Macaque raw H5AD (7,528,695,907 bytes) and imputed MERFISH log2 H5AD (50,206,825,799 bytes) were downloaded and retained locally.

## Defects exposed by the 0.2.0.9005 focused run

### Native HMBA Macaque

`genes = NULL` reached the real H5AD but failed on a duplicated source NCBI identifier whose gene symbol was also duplicated. Native Macaque gene metadata carries multiple identifier systems, including Ensembl and NCBI fields. Candidate 0.2.0.9006 now:

1. reads Ensembl, NCBI, transcript, mapped-NCBI and symbol aliases from metadata and H5AD `var`;
2. reconciles `ENSEMBL:` and `NCBIGene:` prefix variants conservatively;
3. tries stable identifiers before symbols;
4. preserves duplicate source-identifier occurrence and group size;
5. uses occurrence alignment only when both duplicate-group sizes agree;
6. never aggregates duplicate rows.

The focused test additionally performs an explicit fetch of generated `.iabdupN` identifiers.

### Imputed MERFISH

The real log2 H5AD completed a simple fetch and warm-cache repeat. The 0.2.0.9005 authoritative-axis scan then allowed more than one WMB metadata alias to resolve to the same imputed H5AD feature, which produced a duplicate/non-increasing dense HDF5 index. Candidate 0.2.0.9006 now:

1. chooses at most one canonical metadata row per H5AD feature;
2. prefers stable identifier matches to symbol-only aliases;
3. records collapsed metadata aliases in provenance;
4. verifies one-to-one mapping in every selected source;
5. supplies h5py a strictly increasing, duplicate-free dense column index;
6. rejects explicit requests that map two identifiers to one feature.

## Release gate

Run `TEST_FAILED_DATASETS_9007.R` against the retained H5AD files. Do not mark the two routes as live-fetch validated until:

```r
result$all_selected_routes_pass
# TRUE
```

The focused script requires 35,219 measured native Macaque genes and 8,460 imputed MERFISH genes, then tests explicit simple and complex queries, warm-cache equality, overlap equality, duplicate Macaque identifiers, the previous problem genes, optional AnnData comparison, optional SCE/Seurat construction, and exact-file cleanup after success.
