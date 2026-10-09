# Building, testing and maintaining importABCatlas

## 1. What you have, and what you do not yet have

The ZIP contains a replacement source package, synthetic test data, test scripts, a benchmark runner, and documentation. It does not contain the atlas itself. The original ZIP and workbook were audited, not modified in place. The Google manuscript was not accessible because the available Google Drive integration is disabled by the administrator; the publication material is therefore not an edit of that document.

The executable Python subset passed 75 fixture tests, with 3 explicit skips. R and live downloads could not run in the delivery environment. Start with the included synthetic example, then run the live validator. Do not advertise every route as validated until the dataset-specific results pass on a real machine. The data release default is 20260711, based on Allen's current documentation at the audit date, not a freshly downloaded manifest certificate.

There are three distinct things to keep separate: the Git repository holds small source files and their history; your R and Python installations run that source; the atlas cache holds potentially very large data files. Keep these in different folders. Deleting a source checkout should not delete your downloaded atlas. Never place the cache inside the repository.

## 2. Install the prerequisites

Install a current supported R release from CRAN, and optionally RStudio Desktop from Posit. R is the language runtime; RStudio is an editor and interface around it. Install Python >=3.10, Git, and enough free disk space for the source files listed by the manifest. Python 3.11 is the intended CI test version, but that CI has not yet been run. Use a Python/R combination supported by the packages available for your operating system rather than copying an old hard-coded Python 3.10 path.

On Windows, use the Python installer option that exposes Python on PATH, then restart RStudio. Rtools is needed when R packages must be compiled from source; select the Rtools version matching your R installation. On macOS, compiler tools may be needed for source builds. Linux package names vary by distribution. Prefer binary R packages where available. Do not run a package installer as administrator merely to bypass an unexplained error.

In an operating-system terminal, check the installations:

```sh
git --version
python --version
Rscript --version
```

On macOS/Linux the executable may be named `python3`. The R Console accepts R commands such as `library(...)`; the Terminal accepts shell commands such as `git status`. Pasting a shell command into the R Console is a common beginner error. In RStudio these are separate tabs.

## 3. Install and run the supplied candidate locally

Unzip the bundle to a normal writable folder. The package root is the folder containing `DESCRIPTION`, `NAMESPACE`, `R`, `inst`, and `tests`. It is not the parent bundle folder. Replace the example path below with that package root. On Windows use forward slashes, such as `C:/Users/you/Documents/importABCatlas_bundle/importABCatlas`.

```r
install.packages("remotes")
remotes::install_local("/path/to/importABCatlas", dependencies = NA,
                       upgrade = "never")
library(importABCatlas)
setup_environment(install = TRUE)
x <- example_data("matrix")
dim(x$expression)
# Expected: 4 genes by 3 synthetic cells
x$cell_metadata$cell_label
```

The initial `remotes` command installs the R package and its required R dependencies. `dependencies=NA` avoids automatically installing every optional analysis framework. `setup_environment(install=TRUE)` explicitly creates a dedicated Python virtual environment and installs numerical dependencies plus the official Allen SDK. A virtual environment is an isolated package directory: it prevents this project from accidentally changing Python packages used by unrelated work. It is not a data cache and should not be committed to Git.

The setup function does not silently install Python itself. Where discovery is ambiguous, provide the actual executable path:

```r
setup_environment(install = TRUE,
                  python = "/absolute/path/to/python")
```

Use the path returned by `Sys.which("python")`, or the known path of the intended installation. After successful setup, `getOption("importABCatlas.python")` shows the dedicated interpreter selected for this R session. A new R session forgets ordinary options, so call `setup_environment()` again with its default `install=FALSE`; this selects an existing environment rather than reinstalling it.

For noninteractive scripts and CI, set the `IMPORTABCATLAS_PYTHON` environment variable to that exact interpreter. `usethis::edit_r_environ()` can open your personal `.Renviron`; add a line of the form `IMPORTABCATLAS_PYTHON=/absolute/path/to/venv/bin/python`, or the Windows `Scripts/python.exe` equivalent, then restart R. Keep personal environment files out of the repository. A path in `options(importABCatlas.python=...)` takes precedence.

Before a paper release, replace setup's default `upstream_ref="main"` with a tested official SDK commit hash. Main is moving software, not a reproducibility lock. `diagnose()` reports the interpreter and installed versions; save it alongside results. Record `pip freeze` from that same interpreter and an R dependency snapshot only after validation.

## 4. Select an output format

The least demanding output is `return_type="matrix"`. It returns a list, not a naked matrix: use `x$expression`, `x$cell_metadata`, `x$gene_data`, and `x$provenance`. The expression matrix is genes by cells and sparse. Its column names are exact cell labels and its row names are source gene identifiers.

Install optional R frameworks explicitly:

```r
install.packages("SeuratObject")
install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SpatialExperiment"))
example_data("seurat")
example_data("sce")
example_data("spatial")
```

Seurat raw count data enter the counts assay. Transformed data enter a Seurat v5 data layer without fabricated counts. Upstream log2 data do not necessarily follow Seurat's usual natural-log convention: do not blindly call count normalizers, invert with `expm1`, or treat imputed expression as observations. SCE assay names deliberately distinguish `counts`, `raw`, `log2`, and `imputed_log2`. `SpatialExperiment` requires valid anatomical coordinates, not an embedding named x/y that happens to be a UMAP. The package does not attach artificial Visium image objects to MERFISH data.

## 5. The normal atlas workflow

First inspect the bundled routes. Then query the real manifest; this is a network operation, unlike listing the bundled registry.

```r
library(importABCatlas)
setup_environment()
cache <- "~/abc_atlas_cache"
release <- "20260711"
list_datasets()
live_files <- list_files(download_base = cache, release = release)
head(live_files)
```

A manifest is an index of named, versioned assets. Its date is not necessarily each constituent file's publication date: a newer release can reference older unchanged files. Pin the manifest rather than constructing guessed S3 URLs. Case is resolved only when a unique corresponding directory exists in that actual manifest.

For a first experiment, read a metadata preview. This can still download a large entire CSV; the preview limits rows parsed/returned, not bytes transferred.

```r
d <- load_data_MERFISH(download_base = cache, release = release,
                       preview_n = 20L)
names(d$cell_metadata)
head(d$gene_data)
ids <- head(d$cell_metadata$cell_label, 5L)
genes <- head(d$gene_data$gene_identifier, 3L)
plan_fetch(data = d, cell_ids = ids, genes = genes)
```

Read the plan before approving expression downloads. It reports selected source files, selected output dimensions and directory-level size context. Directory sizes are not exact sizes of the selected files and not a promise about network traffic. Five selected cells may still require a multi-gigabyte H5AD file. Confirm disk capacity and an appropriate network connection.

```r
x <- fetch_data(data = d, cell_ids = ids, genes = genes,
                return_type = "matrix", allow_downloads = TRUE,
                max_output_mb = 128)
stopifnot(identical(colnames(x$expression), ids))
x_again <- fetch_data(data = d, cell_ids = ids, genes = genes,
                      return_type = "matrix", offline = TRUE)
```

Offline means missing local files are an error. It does not mean an empty cache will work. `max_output_mb` limits the dense-equivalent selected output; it is not a hard cap on total memory. H5AD index arrays, metadata, sparse intermediates, CSV bridging and R objects also consume memory.

For biological subsets, inspect exact metadata column names and values. Filters are named lists of exact matches, not expressions and not fuzzy searches. For instance, after confirming the source actually contains those values, a filter might select one `feature_matrix_label` and particular cell-type annotations. Missing categories must not be silently invented. A preview and `head()` are schema checks, not representative sampling. For a reproducible selection use reviewed filters and `n_cells=100L, seed=1729L`, save the exact cell IDs, then reload by `cell_ids`.

Never convert cell labels to numbers. Some spatial IDs are far larger than a floating-point integer can exactly represent. Excel can irreversibly round them. Read identifiers as strings directly from the source and keep them as strings through R, Python, JSON and CSV.

## 6. Git and GitHub: the concepts

Git is a local version-control system. GitHub is a website hosting a copy of a Git repository plus collaboration features. You can use Git without GitHub and can run this package without publishing anything.

A repository is a folder with source files and a hidden `.git` directory containing history. A clone is a local copy including that history. `origin` is the conventional name of the remote repository from which you cloned. A branch is a named line of development; create one for this repair so the main branch remains unchanged until review.

Your working tree contains the files you are editing. `git diff` shows changes relative to the last commit. The staging area contains the changes selected for the next commit; `git add` stages them. A commit is a named snapshot recorded locally with an explanatory message. A commit is not automatically on GitHub. `git push` sends local commits to a remote branch.

`git fetch` downloads remote history without merging it into your working branch. `git pull` fetches and integrates the remote branch; this can create a merge or a rebase depending on settings. This guide uses `git pull --ff-only` on a clean main branch: it updates only when no conflicting local history must be reconciled. A pull request is a GitHub review page proposing that one branch be merged into another. It is not the same operation as `git pull`.

A fork is your own GitHub-hosted copy of someone else's repository. Use a fork when you do not have permission to push a branch to the original. A tag is a stable name for a particular commit, such as `v0.2.0`. A GitHub Release adds human-readable notes and downloadable assets to a tag. A DOI archive preserves a specific software version independently of a moving main branch.

## 7. Put the update onto GitHub safely

Configure the author identity for your local commits. These are your real values, not package placeholders.

```sh
git config --global user.name "Your Name"
git config --global user.email "your-address@example.org"
git clone https://github.com/me37uday/importABCatlas.git
cd importABCatlas
git status
git switch main
git pull --ff-only
git switch -c fix/atlas-release-20260711
```

If the repository uses a default branch named something other than `main`, substitute its actual name. Stop when Git reports uncommitted changes or a refused fast-forward; do not force it. Copy personal work elsewhere or commit it on a separate branch before proceeding. Authentication for pushing can use GitHub Desktop's sign-in, an SSH key, or a properly scoped token through a credential manager; never put tokens in scripts, remotes, issue reports or this ZIP.

Keep the unpacked candidate in a separate folder, not inside the clone. The supplied `scripts/apply_candidate.py` has a dry-run mode, requires a clean feature branch, preserves `.git`, and creates a backup before replacement. Run it from the candidate folder, passing absolute paths:

```sh
python scripts/apply_candidate.py --candidate /path/to/candidate/importABCatlas --target /path/to/cloned/importABCatlas
```

Read the proposed replacement list, then add `--apply` to execute it. It replaces the package's owned directories and files, removes named obsolete loaders/docs/Python source from the original layout, and leaves unrelated files alone. The script is supplied for local use and was tested against a disposable local Git repository; it does not contact GitHub, commit, push, or delete your Git history. The original code remains recoverable from both the backup archive and prior commits.

Why not simply paste new files over old files? Old `R/load_data_WMB.R` and similar files can continue defining functions and override the new wrappers depending on load order. Old `man/` help can disagree with the API. Old `README.Rmd` can regenerate an obsolete README. Replacing the named package-owned trees avoids these stale duplicates.

Open `importABCatlas.Rproj` in RStudio. Follow Sections 8-10, then inspect and stage deliberately:

```sh
git status
git diff --stat
git diff
git add R inst man tests scripts docs .github DESCRIPTION NAMESPACE README.md NEWS.md LICENSE .gitignore .Rbuildignore importABCatlas.Rproj
git add -u
git diff --cached --stat
git diff --cached
git commit -m "Add release-aware ABC atlas import and validation"
git push -u origin fix/atlas-release-20260711
```

`git add -u` stages tracked deletions, including removed legacy files. Review that no downloaded data, passwords, logs with personal paths, or virtual environments are staged. The command does not make a licence decision for you. Before a public release, the actual rights holders must replace the licence placeholder and the real maintainers must fill in `DESCRIPTION`.

On GitHub, open the pushed branch and choose Compare & pull request. Explain the observed defects, changes, actual fixture results, and outstanding live/R validation. Do not write 'all datasets tested' when only synthetic tests have passed. Ask a collaborator to review the scientific joins and run live tests on a second machine. Merge only when the required checks and data-specific release gates pass. After merging:

```sh
git switch main
git pull --ff-only
```

For a reviewed release, change the development version to the intended stable version, commit that change, then create and push an annotated tag. Do this only after the release checklist is complete, not as an immediate step on this candidate.

```sh
git tag -a v0.2.0 -m "Validated importABCatlas 0.2.0"
git push origin v0.2.0
```

## 8. Which files to edit

Edit `R/load_data.R` for the public metadata API, `R/fetch_data.R` for R object construction, `R/wrappers.R` for convenience names, `R/setup.R` for environment setup, and `R/backend.R` for the subprocess bridge. Edit `inst/python/registry.json` for route configuration and `inst/python/abc_backend.py` for cache resolution, table joins and H5AD selection. `tests/python/` checks numerical/data logic; `tests/testthat/` checks the R interface and objects.

Lines beginning `#'` immediately above an R function are roxygen2 source documentation. They describe parameters, return values, exports and examples. Roxygen converts them into `man/*.Rd` help files and `NAMESPACE`. A generated-file warning means edit the source comments and regenerate, not that you are forbidden to change the package.

This candidate includes a clearly labelled bootstrap help file and namespace because R was unavailable. Regenerate them once, from the package root:

```r
install.packages(c("devtools", "roxygen2", "testthat", "rcmdcheck"))
source("scripts/regenerate_docs.R")
devtools::load_all()
devtools::test()
devtools::check(args = "--no-manual")
```

On later edits use `devtools::document()`. Check and commit the resulting `man/` and `NAMESPACE` changes together with the source. `devtools::load_all()` loads the checkout for development; it does not permanently install it. A new plain R session with `library(importABCatlas)` may still load an older installed copy. Use `packageVersion()` and `find.package()` to confirm what is running, or reinstall the local source after changes.

## 9. Test at several levels

Use the dedicated Python interpreter, not a different system Python, for the Python tests. In a terminal, run its full path followed by the arguments below, or activate that environment first. Install the additional test/reference libraries once:

```sh
python -m pip install pytest anndata -r scripts/requirements-benchmark.txt
python -m pytest tests/python -q
```

The bundled suite tests exact huge IDs, missing/duplicate keys, safe joins, taxonomy pivots, file disambiguation, dense/CSR/CSC H5AD layouts, order preservation, gene selection, imputed assay semantics, download consent and JSON/MatrixMarket bridging. It includes synthetic routing checks for all 23 registry entries. Synthetic files cannot establish that every real upstream schema matches the assumptions.

Run `devtools::test()` for R and optional object construction. Skipped optional package tests are not passes. Install SeuratObject and the Bioconductor packages, then rerun. `R CMD check` also examines package metadata, documentation and dependencies; a passing Python suite does not imply R package correctness. The supplied CI checks Linux, macOS and Windows but has not been executed on your GitHub repository.

For live validation, set `IMPORTABCATLAS_PYTHON`, then run from the package root:

```sh
Rscript scripts/validate_live.R
Rscript scripts/validate_live.R --metadata
Rscript scripts/validate_live.R --expression --allow-downloads
```

The first command records a live catalogue only. The second may download large metadata files. The third can download large raw and transformed expression files and requires `anndata` plus all optional R frameworks. Limit a first run by setting `ABCATLAS_DATASETS=MERFISH` in the process environment. In R, a portable alternative is:

```r
Sys.setenv(ABCATLAS_DATASETS = "MERFISH",
           ABCATLAS_CACHE = path.expand("~/abc_atlas_cache"),
           ABCATLAS_RELEASE = "20260711",
           ABCATLAS_PYTHON = getOption("importABCatlas.python"))
system2("Rscript", c("scripts/validate_live.R", "--metadata"))
```

The expression validator compares a small selected matrix against an independent AnnData reader, then checks the R object axes, values and cell annotations. It records failures per route rather than replacing them with success. A preview of one source file is not coverage of every anatomical shard, donor, taxonomy class, or neuronal/non-neuronal cohort. Extend reviewed selections across multiple file groups before a release. Benchmarking is separate and must follow `docs/BENCHMARK_PROTOCOL.md`.

## 10. Debug without guessing

Start with the smallest reproducible example: one dataset, pinned release, one source-file group, five exact cells, three native genes, and matrix output. This isolates upstream/network problems from Seurat or Bioconductor problems. Keep the original error message and traceback. The bridge retains a temporary diagnostic directory after failure containing the request and Python log; copy it to a private working folder before closing the session. Inspect for personal paths or sensitive metadata before sharing.

```r
diagnose()
sessionInfo()
packageVersion("importABCatlas")
find.package("importABCatlas")
traceback()
```

A 'Python not found' error means the interpreter path is missing or points to a removed environment. A missing module usually means the wrong Python was selected or dependencies were installed elsewhere. Do not reinstall R as the first response. Compare `diagnose()` with the interpreter receiving `python -m pip install`.

A missing manifest directory or ambiguous file means the schema resolver lacks a verified match. Query `list_files()` for the same release, inspect the source documentation, and use an explicit reviewed override. Do not change the default release to latest just to silence the error. A missing gene is not a zero-expression gene: confirm whether the identifier is in the measured panel, imputed panel, species-native universe or aligned orthologue set.

A duplicate join-key error prevents row multiplication. Determine whether the annotation is genuinely one record per cell. Electrophysiology sweeps and polygon vertices often are not; keep them in separate tables instead of removing duplicates arbitrarily. An unknown taxonomy label should remain unknown until the mapping is verified. Missing spatial coordinates require deliberate cell exclusion or a different output format, not fabricated coordinates.

For a download/hash error, check connectivity, disk space, cache permissions, release and the exact named source file. Preserve logs. The official SDK owns integrity checking; do not disable checks or substitute an unrelated similarly named H5AD. A interrupted transfer may require repairing only that asset using documented SDK behaviour. Never delete an entire shared cache as the first troubleshooting step.

For R debugging, `debugonce(load_data)` enters a step-through browser on the next call; `n` advances, `c` continues and `Q` quits. Most numerical logic is in Python, so inspect the retained request/log or add a small pytest case there rather than trying to step inside Python from the R browser. After changing code, add a regression test that fails on the old behaviour and passes on the intended fix.

## 11. Resolve Git problems without losing work

A rejected push usually means the remote changed or you lack permission. Do not use force push as a beginner workaround. Fetch and inspect the history, then discuss the change with collaborators. For a feature branch, a reviewed merge from the updated main branch is a straightforward option. A merge conflict marks the exact overlapping text with `<<<<<<<`, `=======`, and `>>>>>>>`; edit the file to the desired combined result, remove the markers, test, stage and commit it.

`git restore path/to/file` discards uncommitted edits to that file; use it only when you intend that loss. To undo an already shared faulty commit, prefer `git revert COMMIT_ID`, which creates a new correcting commit and preserves history. `git reset --hard` and forced pushes can destroy or rewrite work and are not needed for this workflow. Before any uncertain operation, make a backup and inspect `git status`.

## 12. A durable release and reproducibility record

Fill in real authors, roles, maintainer email, code licence and change log. Run documentation generation, Python tests, R tests, R CMD check and live comparisons. Archive all source-file selections, manifests, software versions, optional framework versions, platform details and reviewed cohort definitions used in benchmarks. Commit small reproducibility materials only when source terms allow them; keep bulk data in its authorised repository.

Create a stable tag and an archival DOI only after the code is validated. Test installation from a clean machine following the README, not merely from your configured developer session. Keep the package available and maintained for at least the journal's required period. The manuscript must distinguish access convenience from algorithmic novelty, and measured results from hoped-for improvements. Source and submission-policy references are in `docs/SOURCES.md` and the separate publication notes.
