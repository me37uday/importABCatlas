# Building, testing and maintaining importABCatlas

## 1. The three folders to keep separate

There are three different things:

1. the Git repository, which contains source code;
2. the private Python environment, which importABCatlas manages automatically;
3. the atlas cache, which can contain very large downloaded data files.

Never place the atlas cache inside the Git repository. Deleting a source checkout
should not delete downloaded atlas data, and deleting a cache should not delete
source history.

## 2. Install R and the package

Install a supported R release and, optionally, RStudio. Git is needed only if you
want to contribute code. You do not need to install Python, Conda, or Homebrew
Python for normal use.

Install importABCatlas from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "me37uday/importABCatlas",
  dependencies = NA,
  upgrade = "never"
)
```

Restart R or RStudio, then run:

```r
library(importABCatlas)
setup_importABCatlas()
```

`setup_importABCatlas()` creates a private Python 3.11 environment using
reticulate and uv. It does not require administrator access and does not modify
system Python. The setup happens once and the environment is reused. If its
cache is removed, the package recreates it automatically.

The first function that needs Python also runs setup automatically, so the
explicit setup command is optional. It is recommended because the user can see
setup progress before starting a data download.

Check the installation:

```r
diagnose()
x <- example_data("matrix")
dim(x$expression)
# Expected: 4 genes by 3 cells
```

If setup fails because the internet connection was interrupted:

```r
setup_importABCatlas(repair = TRUE)
```

No Python path should be added to `.Renviron` or `.Rprofile` for ordinary use.
Old `IMPORTABCATLAS_PYTHON` or `RETICULATE_PYTHON` values are ignored in
managed mode, and `setup_importABCatlas()` resets stale package-specific
external mode for the current R session.

## 3. Choose an output format

The simplest output is `return_type="matrix"`. It returns a list with:

```r
x$expression
x$cell_metadata
x$gene_data
x$provenance
```

Optional frameworks can be installed separately:

```r
install.packages("SeuratObject")
install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SpatialExperiment"))
```

Then use `return_type="seurat"`, `"sce"`, or `"spatial"`. Start with matrix
output when checking a new machine.

## 4. The normal atlas workflow

List the package's supported routes:

```r
list_datasets()
```

Load metadata first. This does not download expression:

```r
cache <- "~/abc_atlas_cache"

d <- load_data(
  dataset = "Zhuang4",
  download_base = cache,
  release = "20260711",
  preview_n = 1000L
)
```

Select cells and inspect the exact expression download:

```r
cells <- head(d$cell_metadata$cell_label, 5L)

plan <- plan_fetch(
  data = d,
  cell_ids = cells,
  genes = NULL,
  data_type = "raw"
)

print(plan)
```

The plan reports the exact H5AD file, exact object size, cache state, retained
partial bytes, additional transfer, and disk-space result. It does not download
expression and does not substitute a folder total for a file size.

After reviewing the plan:

```r
x <- fetch_data(
  data = d,
  cell_ids = cells,
  genes = NULL,
  data_type = "raw",
  return_type = "matrix",
  allow_downloads = TRUE
)
```

A request for only a few cells can still require a complete multi-gigabyte H5AD.
Downloads show progress, speed, and ETA. Interrupted files are retained with a
`.part` suffix and resume when the same cache is used again.

After the first `genes=NULL` query, use the returned measured identifiers for
strict smaller queries:

```r
genes <- head(x$gene_data$gene_identifier, 20L)

y <- fetch_data(
  data = d,
  cell_ids = cells,
  genes = genes,
  data_type = "raw",
  return_type = "matrix",
  allow_downloads = FALSE,
  offline = TRUE
)
```

Keep cell and gene identifiers as character strings. Some cell labels are too
large to represent safely as numbers, and spreadsheet software can round them.

## 5. Troubleshooting

Run:

```r
diagnose()
```

If Python is not ready:

```r
setup_importABCatlas(repair = TRUE)
```

If a data download stops, rerun the same command with the same `download_base`.
Do not delete the `.part` file; it is used for resumption.

When reporting a problem, include the complete command, complete error,
`packageVersion("importABCatlas")`, `diagnose()`, operating system, dataset,
release, cache state, and free disk space.

Advanced developers who intentionally need an existing Python can opt out of the
managed environment for the current R session:

```r
setup_environment(python = "/absolute/path/to/python")
```

Normal users should not need this.

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

The package itself manages Python automatically. The standalone Python fixture tests are developer tests and can run in any compatible development Python. Install their additional test/reference libraries once:

```sh
python -m pip install pytest anndata -r scripts/requirements-benchmark.txt
python -m pytest tests/python -q
```

The bundled suite tests exact huge IDs, missing/duplicate keys, safe joins, taxonomy pivots, file disambiguation, dense/CSR/CSC H5AD layouts, order preservation, gene selection, imputed assay semantics, download consent and JSON/MatrixMarket bridging. It includes synthetic routing checks for all 23 registry entries. Synthetic files cannot establish that every real upstream schema matches the assumptions.

Run `devtools::test()` for R and optional object construction. Skipped optional package tests are not passes. Install SeuratObject and the Bioconductor packages, then rerun. `R CMD check` also examines package metadata, documentation and dependencies; a passing Python suite does not imply R package correctness. The supplied CI checks Linux, macOS and Windows but has not been executed on your GitHub repository.

For live validation, run `setup_importABCatlas()` once, then use the package scripts from the repository root:

```sh
Rscript scripts/validate_live.R
Rscript scripts/validate_live.R --metadata
Rscript scripts/validate_live.R --expression --allow-downloads
```

The first command records a live catalogue only. The second may download large metadata files. The third can download large raw and transformed expression files and requires `anndata` plus all optional R frameworks. Limit a first run by setting `ABCATLAS_DATASETS=MERFISH` in the process environment. In R, a portable alternative is:

```r
Sys.setenv(ABCATLAS_DATASETS = "MERFISH",
           ABCATLAS_CACHE = path.expand("~/abc_atlas_cache"),
           ABCATLAS_RELEASE = "20260711")
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
