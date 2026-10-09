# Architecture and invariants

## Data flow

R API -> JSON request in a temporary directory -> explicit external Python interpreter -> official AbcProjectCache -> validated metadata / selected H5AD reads -> CSV plus column schema and MatrixMarket -> R sparse matrix / optional analysis object.

The Allen SDK owns release manifests, cache paths, downloads and integrity handling. importABCatlas does not maintain guessed S3 URLs, a competing downloader, or hard-coded installation of an interpreter. Python is an implementation dependency; the public workflow is R. The bridge starts a Python process for each call, so startup, CSV transfer and sparse conversion overhead are real costs included in function timing.

## Scientific invariants

Cell and gene identifiers remain exact strings. Requested cell/gene order is restored after indexed reads. Each annotation join must be many-to-one and preserve the number and order of cells. Missing and ambiguous genes cause errors, not zero imputation. Taxonomy aliases are source keys, never generated row numbers. Colliding metadata columns receive source-specific suffixes rather than overwriting values. Manifest mismatches between metadata, fetches and supplemental joins fail.

Feature symbols are annotations, not guaranteed unique matrix keys. Native species genes are not interchangeable with an aligned orthologue universe. Each source H5AD supplies the authoritative expression index. Default reads use X, not an automatically guessed raw/layer slot. The importer refuses ambiguous routing and has explicit overrides for reviewed schema changes.

Dense, CSR and CSC encodings are handled without converting the entire expression matrix to dense. Sparse output is genes by cells. Dense-equivalent selected-output size is checked before conversion; full observation/feature index arrays and some metadata are still resident. This is selective expression access after whole-file download, not remote cell-level streaming and not a proof of constant memory.

## Assay and spatial semantics

Only integer, nonnegative raw values can default to a counts label, and users can explicitly prevent that label with counts_semantics=FALSE. Numerical integrality is necessary but not sufficient evidence for biological count semantics; review upstream methods. Imputed or log2 products cannot be forced into counts. Log2 normalization is retained and named rather than silently reinterpreted as natural-log normalisation.

Anatomical coordinates are stored with original metadata and their selected column names in provenance. A UMAP is not a physical coordinate system. Coordinate units, coordinate reference, specimen and section identifiers remain the source's responsibility; no automatic coordinate registration is performed. Related but non-cell-level resources are retrieved separately.

## Configuration and extension

The registry is data, not separate duplicate code for each atlas. A generic directory route can use verified overrides. Dataset-specific exceptional rules currently include Aging mouse mappings, SEA-AD alternate expression identifiers and native-species HMBA metadata discovery. Files must first exist in the selected manifest. A unique case-insensitive directory match is allowed only against that manifest, because upstream release notes and notebooks use inconsistent capitalisation in places.

Temporary bridge files are removed on success and retained with error details on failure. Atlas data remain in the configured persistent cache. R failures after a successful Python computation also leave bridge output for debugging. Logs can include paths or selected metadata and must be reviewed before public sharing.

## Deliberate limitations

No biological reclassification, automatic neuron/non-neuron inference, statistical testing, batch correction, imputation, cross-species orthology inference or image reconstruction is performed. No automatic ATAC/SCE multimodal construction is implied by a source name containing Multiome. No promise covers resources outside the selected ABC manifest. No full-atlas or R verification was executed in this delivery environment.
