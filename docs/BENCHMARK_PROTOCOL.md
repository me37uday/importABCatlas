# Genuine cold/warm benchmarks and Figure 1

## Current evidence

No live biological runtime, download volume or resident-memory measurement was obtained in this environment. The original timing sheet is blank. The bundle therefore contains a runnable protocol and plotting code, not an invented figure. Unit-test duration is a software-test result, not an atlas-import benchmark.

## Scientific question and units

Measure two public functions separately: `load_data` (manifest/cache resolution, metadata retrieval, validated joins and R metadata object construction) and `fetch_data` (source-file retrieval, selected expression reads, transfer and the requested R analysis-object construction). Each call starts its own Python subprocess. That startup and bridge overhead belongs in the measured function time; R process startup and the worker's RDS save/read do not.

The natural experimental unit is a fixed dataset/representation, cohort definition, selected cell list, selected gene list and output format. A biological cohort such as non-neuronal cells must be defined from the actual source taxonomy. Do not match a substring like 'NN' across all atlases; WMB's NN-IMN-GC grouping is not synonymous with non-neuronal cells. Unknown assignments and developmental precursors need reviewed treatment. Patch-seq non-neuronal is N/A, not a zero-second result.

For a compact primary experiment, use 100 cells and a fixed small panel such as three available genes per dataset. The same exact identifiers must be used for cold and warm runs. A different native panel across datasets is acceptable for an access benchmark if disclosed; it does not establish a harmonized biological comparison. Use a second scaling experiment with larger selections if claiming memory or size scaling. Do not infer whole-atlas feasibility from a 100-cell output.

## Define cold and warm precisely

A cold replicate gets its own new empty application cache. Time load_data first, then fetch_data; at fetch time metadata and the manifest are already cached, but the selected expression is not. A warm replicate uses that same populated cache but fresh R processes and reruns the same functions and selections. There is no substitution of a serialized finished object for fetch_data.

The operating-system page cache, network routes and server caches are not cleared. Therefore describe 'cold application cache', not a physically cold machine. Warm calls remain in the same online mode as cold calls so changing offline semantics does not confound the comparison. The runner suppresses warm measurements when the preceding cold workflow failed rather than mislabelling a partial download as a warm cache.

Each dataset/cohort/technical replicate receives a distinct cache, preventing a previous non-neuronal run from warming the neuronal run. The cost is repeated large downloads and disk use. The script never deletes any cache. Review the entire plan and available storage before enabling many routes. A sequential run order is supplied; record that order and consider balanced randomisation for a stronger study because network conditions can drift.

## Prepare reviewed selections outside timing

Edit `scripts/benchmark_config.json`. All 46 potential dataset/cohort profiles initially have enabled=false and cohort_reviewed=false. Keep absent cohorts disabled and document N/A in the coverage record. Set exact source-column filters, choose a representation, genes and output type, and have a collaborator review the biological definition before enabling a profile. Except for the known neuronal Patch-seq route, an empty filter is refused.

For SEA-AD Multiregion the source documents a Class field with neuronal GABAergic/Glutamatergic and non-neuronal/non-neural values. Those exact values are shown in DATASET_GUIDE.md. Do not copy them into WMB, PMDBS or developing-mouse profiles without inspecting those source taxonomies. PMDBS mapping confidence and unassigned cells need an explicit policy. Freeze classifications and select cells once; do not resample for each timing replicate.

```sh
Rscript scripts/prepare_benchmark.R scripts/benchmark_config.json
```

This can download metadata into a separate selection cache and can require substantial RAM. It writes `benchmark_profiles.json` containing exact string identifiers, seed, filters, release, representation and output format. These preparation times are not benchmark observations. Inspect each frozen selection's source file plan with plan_fetch before authorising expression downloads. In particular, a random 100-cell whole-mouse sample may span many expression shards, while another atlas's sample fits in one; source layout is part of the result and must be reported.

## Run real measurements

Ensure the selected Python interpreter has psutil, pandas and matplotlib, and the R package plus chosen R framework is installed. Ensure `ABCATLAS_PYTHON` is exported for child R processes; an option set only in an interactive R session is not inherited by a new process.

```sh
python -m pip install -r scripts/requirements-benchmark.txt
python scripts/benchmark.py --profiles benchmark_profiles.json --cache-root /large-disk/abc_benchmark_caches --out benchmark_results --replicates 3 --allow-downloads
```

Use the appropriate Python executable path on your machine. The runner creates a unique run directory rather than overwriting earlier results. No numeric result is written until a real worker returns it. Failed stages retain logs and a failure/NOT_RUN status, not a runtime of zero. The filename printed at the end identifies the actual results CSV.

The recorded `seconds` is elapsed wall time inside the R function. `sampled_peak_process_tree_rss_mb` is the largest summed R-plus-descendant resident-memory reading observed at roughly 50 ms sampling while the stage flag is running. It can miss a transient peak; it is not a theoretical memory bound. `net_cache_growth_bytes` sums changes in cached file lengths, not on-the-wire network bytes and not physical allocated disk blocks. Do not relabel these metrics as exact traffic or exact peak RSS.

Keep full worker result JSON files, which record function provenance and R session information, together with source plans, package/SDK commit identifiers, selection hashes, machine CPU/RAM/storage details and network description. Record disk free space before/after and any interruptions. Three technical repetitions quantify machine/workflow variability, not biological uncertainty across donors.

## Generate one figure only after validation

```sh
python scripts/plot_benchmarks.py benchmark_results/ACTUAL_RUN/benchmark_results.csv --output figure1.png
```

The script produces a four-column heatmap for load-cold, load-warm, fetch-cold and fetch-warm, plus a vector SVG and summary CSV. Its logarithmic seconds scale accommodates different source-file layouts. It rejects synthetic-origin data, nonpositive/nonfinite values, inconsistent selections, mixed formats/transformations/manifests within a cohort, duplicate conditions and a lack of complete cold/warm pairs. It summarizes complete successful pairs; missing conditions remain blank. It cannot detect a human deliberately falsifying a LIVE_ATLAS label, so retain and review raw worker evidence.

The plotted output must be inspected at the journal's final width. With many dataset/cohort rows a four-column heatmap is more economical than dozens of separate plots, but the figure may still be too tall. Use compact readable route names, explain abbreviations, and put detailed per-repeat values in supplementary data. Do not silently omit difficult datasets; state exclusions and completeness explicitly. The supplied figure generator has not been run on biological timings because there are none here.

A suitable author-written legend should identify the machine, software/data versions, exact cell/gene selection sizes, cohort definitions, output object, three technical repetitions, median summary, log scale, cold-cache definition, and N/A/missing conditions. Do not present the draft legend as evidence that those measurements already exist.

## Baseline, interpretation and submission threshold

Before claiming improvement, compare equivalent outputs with the official SDK plus explicit backed AnnData selection and R conversion, or another justified existing workflow. Match cells, genes, source manifests, transformed/raw values and output format; include the download and cache state symmetrically. Validate exact matrix equivalence before comparing speed. The independent AnnData checks in validate_live.R establish a correctness comparison for small subsets, not a complete performance baseline.

A within-tool cold/warm comparison demonstrates caching behaviour, not superiority to other software. Differences between neuronal and non-neuronal samples may reflect shard layout, sparsity, panel size or annotation joins rather than neuronal biology. The hypothesis that whole-mouse neurons will be slowest must remain a hypothesis until the measurements support it. First-run source transfer can dominate small-output import; that finding may be more informative than a claimed tiny object-construction speedup.

Only replace the manuscript's pending-results statements after all biological runs and exclusions are reviewed. Keep the local fixture test report separate from the scientific benchmark table. The reporting workbook intentionally has blank timing cells and explicit NOT_MEASURED states until genuine measurements are supplied.
