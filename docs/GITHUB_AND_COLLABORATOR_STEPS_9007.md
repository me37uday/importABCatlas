# Put 0.2.0.9007 on GitHub and ask collaborators to test

## 1. Test the candidate locally

Install from the extracted source folder in a fresh R session:

```r
install.packages("remotes")
remotes::install_local(
  "/absolute/path/to/importABCatlas-0.2.0.9007",
  dependencies = NA,
  upgrade = "never",
  force = TRUE
)
```

Restart R or RStudio, then run:

```r
library(importABCatlas)
setup_importABCatlas()
diagnose()

example <- example_data("matrix")
stopifnot(identical(dim(example$expression), c(4L, 3L)))
```

No Python path should be set. The package creates and reuses its own managed
Python environment.

Run the automatic-environment regression:

```r
source("/absolute/path/to/importABCatlas-0.2.0.9007/TEST_AUTOMATIC_SETUP_9007.R")
test_automatic_setup_9007(run_live_plan = TRUE)
```

Run the focused live feature-axis regression when the large source files are
already cached, or set `download_missing = TRUE` only after reviewing disk and
network requirements:

```r
source("/absolute/path/to/importABCatlas-0.2.0.9007/TEST_FAILED_DATASETS_9007.R")

result <- test_failed_datasets_9007(
  cache_root = "~/importABCatlas_9004_large_fetch_cache",
  download_missing = FALSE,
  cleanup_expression_after_success = TRUE,
  output_prefix = "~/Desktop/importABCatlas_9007_failed_routes"
)

result$all_selected_routes_pass
```

## 2. Apply the candidate to a Git branch

From a clean clone:

```sh
git switch main
git pull --ff-only
git switch -c automatic-python-9007
```

Extract the candidate outside the clone. Copy it into the clone while
preserving `.git`:

```sh
rsync -av --delete \
  --exclude='.git' \
  --exclude='*.Rcheck' \
  --exclude='*.tar.gz' \
  /absolute/path/to/importABCatlas-0.2.0.9007/ \
  /absolute/path/to/importABCatlas-clone/
```

Review the changes and verify that no downloaded Allen files are present:

```sh
git status --short
git diff --stat
git diff
find . -type f -size +50M -not -path './.git/*' -print
```

## 3. Run release gates

From the repository root, in R:

```r
devtools::test()
devtools::document()
```

In Terminal:

```sh
python -m pytest -q tests/python
R CMD build .
R CMD check --no-manual importABCatlas_0.2.0.9007.tar.gz
```

The bundled Python suite should pass. Review every R check error or warning
before pushing.

## 4. Commit and open a pull request

```sh
git add -A
git status
git diff --cached --stat
git commit -m "Automate Python setup and simplify installation"
git push -u origin automatic-python-9007
```

Create a GitHub pull request from `automatic-python-9007` into `main`. Wait for
all automated checks, inspect **Files changed**, and merge only after checks
pass.

## 5. Collaborator installation

Before merge, collaborators can install the branch directly:

```r
install.packages("remotes")
remotes::install_github(
  "me37uday/importABCatlas@automatic-python-9007",
  dependencies = NA,
  upgrade = "never",
  force = TRUE
)
```

After restarting R:

```r
library(importABCatlas)
setup_importABCatlas()
diagnose()

example <- example_data("matrix")
stopifnot(identical(dim(example$expression), c(4L, 3L)))
```

Collaborators do not need to install Python, Conda, or Homebrew, and should not
set a Python path. Ask each collaborator to use a separate atlas cache and test
at least one live route.

## 6. What collaborators should report

For every issue, include:

- Git commit and package version;
- operating system and architecture;
- R version;
- complete `diagnose()` output;
- Allen release and dataset route;
- raw or log2 representation;
- exact command and complete error;
- whether the source file was new, partial, or cached;
- available disk space;
- retained diagnostic-directory path.

Do not upload Allen H5AD files or sensitive donor/cell-level data to GitHub.
