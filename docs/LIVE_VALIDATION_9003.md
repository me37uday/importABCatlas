# Live validation and upstream accessibility (release 20260711)

This table consolidates the original all-route live run and the focused 0.2.0.9002 rerun. It distinguishes a completed expression fetch from a route that was deliberately not downloaded because its exact selected H5AD exceeded the configured test cap.

## Conclusion

- Every one of the 23 bundled routes completed live metadata loading and exact-file `plan_fetch()` routing against Allen release `20260711`.
- Sixteen routes completed live simple and complex `fetch_data()` tests, warm-cache equality, and axis checks.
- Seven routes were not expression-fetched in the safe rerun because the selected source files ranged from approximately 7 to 47 GiB. These are `MERFISH_imputed`, `WHB`, `PMDBS`, `HMBA`, `HMBA_Human`, `HMBA_Macaque`, and `SEAAD_CaH`. Exact S3 object existence and size were verified; this is `NOT_RUN_SIZE_CAP`, not a functional failure.
- No dataset route is removed. All map to published Allen resources and exact manifest objects in release `20260711`.
- Two representations are absent from the validated manifest and are therefore removed from the advertised choices: raw expression for `MERFISH_imputed`, and raw expression for `HMBA_PatchSeq`. Both routes remain accessible as log2-only expression routes.

## Route table

| Route | Modality | Supported expression | Default | Load | Plan | Live fetch | Primary source GiB |
|---|---|---|---|---|---|---|---:|
| `WMB` | `single_cell_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 6.90 |
| `ConsensusMouseAIBS` | `single_cell_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 5.68 |
| `ConsensusMouseMacosko` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 3.51 |
| `AgingMouse` | `single_cell_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 12.90 |
| `DevelopingMouse` | `single_cell_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 6.83 |
| `MERFISH` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 7.10 |
| `MERFISH_imputed` | `imputed_spatial_transcriptomics` | `log2` | `log2` | PASS | PASS | NOT_RUN_SIZE_CAP | 46.76 |
| `Zhuang1` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 1.11 |
| `Zhuang2` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.48 |
| `Zhuang3` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.61 |
| `Zhuang4` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.06 |
| `WHB` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | NOT_RUN_SIZE_CAP | 30.61 |
| `PMDBS` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | NOT_RUN_SIZE_CAP | 14.68 |
| `HMBA` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | NOT_RUN_SIZE_CAP | 14.34 |
| `HMBA_Human` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | NOT_RUN_SIZE_CAP | 14.87 |
| `HMBA_Macaque` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | NOT_RUN_SIZE_CAP | 7.01 |
| `HMBA_Marmoset` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 3.61 |
| `HMBA_spatial_Human` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.63 |
| `HMBA_spatial_Macaque` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.67 |
| `HMBA_spatial_Marmoset` | `spatial_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.29 |
| `HMBA_PatchSeq` | `patch_seq_transcriptomics` | `log2` | `log2` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 0.06 |
| `SEAAD_Multiregion` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | PASS_LIVE_SIMPLE_COMPLEX | 2.40 |
| `SEAAD_CaH` | `single_nucleus_transcriptomics` | `raw,log2` | `raw` | PASS | PASS | NOT_RUN_SIZE_CAP | 10.73 |

## Upstream accessibility evidence

Official Allen documentation describes the release as a public AWS dataset and states that expression H5AD files must be downloaded locally. The package therefore cannot prove a 30–50 GiB route by returning a tiny subset without first transferring the complete selected H5AD.

- ABC Atlas data access and release notes: https://alleninstitute.github.io/abc_atlas_access/
- Imputed MERFISH: https://alleninstitute.github.io/abc_atlas_access/descriptions/MERFISH-C57BL6J-638850-imputed.html
- HMBA individual species: https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-10XMultiome-BG.html
- HMBA spatial: https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-Spatial-BG.html
- HMBA Macaque Patch-seq: https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-Macaque-PatchSeq-BG.html
- PMDBS: https://alleninstitute.github.io/abc_atlas_access/descriptions/ASAP-PMDBS-10X.html
- WHB: https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html

## Interpretation of status values

- `PASS_LIVE_SIMPLE_COMPLEX`: simple and complex expression requests completed, axes matched, and the warm-cache result was identical.
- `NOT_RUN_SIZE_CAP`: metadata and exact file planning passed, but expression transfer was intentionally skipped by the validation cap. It must not be described as a passed fetch or as a failed route.

