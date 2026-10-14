# Regenerate the precomputed vignette from its source, importABCatlas.Rmd.orig.
#
# The vignette code downloads several GB from the ABC Atlas, so it is run here,
# by hand, rather than at package build time. Run from the R_only/ folder,
# pointing ABC_ATLAS_CACHE at a download folder with enough space:
#
#   ABC_ATLAS_CACHE=/data/abc_atlas Rscript vignettes/precompile.R

pkgload::load_all(".", helpers = FALSE, attach_testthat = FALSE, quiet = TRUE)
knitr::knit("vignettes/importABCatlas.Rmd.orig", "vignettes/importABCatlas.Rmd")
