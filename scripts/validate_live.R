# Run from the package root after installation.
# Default = manifest catalogue only. --metadata may download GB of metadata.
# --expression permits full H5AD downloads and requires --allow-downloads.
library(importABCatlas)
args <- commandArgs(trailingOnly = TRUE)
metadata_test <- "--metadata" %in% args || "--expression" %in% args
expression_test <- "--expression" %in% args
allow <- "--allow-downloads" %in% args
if (expression_test && !allow) stop("Expression tests require --allow-downloads; first review plan_fetch().")
cache <- Sys.getenv("ABCATLAS_CACHE", "~/abc_atlas_cache")
release <- Sys.getenv("ABCATLAS_RELEASE", "20260711")
out <- Sys.getenv("ABCATLAS_VALIDATION", "validation-live")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(diagnose(), file.path(out,"environment.json"), auto_unbox=TRUE, pretty=TRUE)
files <- list_files(download_base=cache,release=release)
write.csv(files,file.path(out,"live-catalog.csv"),row.names=FALSE)
registry <- list_datasets()
selected <- Sys.getenv("ABCATLAS_DATASETS", "")
if (nzchar(selected)) registry <- registry[registry$dataset %in% strsplit(selected,",",fixed=TRUE)[[1]],]
rows <- list()
record <- function(id,stage,status,message="",elapsed=NA_real_) {
  rows[[length(rows)+1L]] <<- data.frame(dataset=id,stage=stage,status=status,message=message,
     elapsed_seconds=elapsed,manifest=release,data_origin="LIVE_ATLAS",stringsAsFactors=FALSE)
  write.csv(do.call(rbind,rows),file.path(out,"results.csv"),row.names=FALSE,na="")
}
for (id in registry$dataset) {
  if (!metadata_test) {record(id,"load_data","NOT_RUN","Catalogue-only mode"); next}
  start <- proc.time()[["elapsed"]]
  d <- tryCatch(load_data(id,download_base=cache,release=release,preview_n=100L),error=function(e)e)
  if (inherits(d,"error")) {record(id,"load_data","FAIL",conditionMessage(d));next}
  stopifnot(!anyDuplicated(d$cell_metadata$cell_label),is.character(d$cell_metadata$cell_label))
  record(id,"load_data","PASS_METADATA_PREVIEW",elapsed=proc.time()[["elapsed"]]-start)
  if (!expression_test) {record(id,"fetch_data","NOT_RUN","No expression consent");next}
  ids <- head(d$cell_metadata$cell_label,5L)
  genes <- if (id=="MERFISH_imputed") c("Calb2","Baiap3","Lypd1") else head(d$gene_data$gene_identifier,3L)
  types <- if (id=="MERFISH_imputed") "log2" else c("raw","log2")
  for (kind in types) {
    start <- proc.time()[["elapsed"]]
    attempt <- tryCatch({
      p <- plan_fetch(data=d,cell_ids=ids,genes=genes,data_type=kind)
      jsonlite::write_json(p,file.path(out,paste0(id,"-",kind,"-plan.json")),auto_unbox=TRUE,pretty=TRUE)
      x <- fetch_data(data=d,cell_ids=ids,genes=genes,data_type=kind,return_type="matrix",
                      allow_downloads=TRUE,verify_reference=TRUE,max_output_mb=32)
      stopifnot(identical(colnames(x$expression),ids),isTRUE(x$provenance$reference_comparison))
      # Independent source comparison above, then exercise all object constructors.
      formats <- c("seurat","sce",if (isTRUE(d$provenance$spec$spatial)) "spatial")
      for (format in formats) {
        y <- fetch_data(data=d,cell_ids=ids,genes=genes,data_type=kind,return_type=format,
                        offline=TRUE,max_output_mb=32)
        stopifnot(identical(dim(y),dim(x$expression)),identical(colnames(y),colnames(x$expression)),
                  identical(rownames(y),rownames(x$expression)))
        ym <- if (format=="seurat") SeuratObject::LayerData(y, assay="RNA",
                 layer=if (x$provenance$assay_kind=="counts") "counts" else "data") else
                 SummarizedExperiment::assay(y, 1L)
        stopifnot(isTRUE(all.equal(as.matrix(ym), as.matrix(x$expression), tolerance=1e-8)))
        yc <- if (format=="seurat") y[[]] else as.data.frame(SummarizedExperiment::colData(y))
        for (column in names(x$cell_metadata))
          stopifnot(identical(as.character(yc[[column]]), as.character(x$cell_metadata[[column]])))
      }
      jsonlite::write_json(x$provenance,file.path(out,paste0(id,"-",kind,"-provenance.json")),auto_unbox=TRUE,pretty=TRUE)
      TRUE
    },error=function(e)e)
    if (inherits(attempt,"error")) record(id,paste0("fetch_data/",kind),"FAIL_OR_UNAVAILABLE",conditionMessage(attempt)) else
      record(id,paste0("fetch_data/",kind),"PASS_SOURCE_COMPARISON_AND_R_OBJECTS",elapsed=proc.time()[["elapsed"]]-start)
  }
}
cat("Results are in",normalizePath(out),"\nA small preview test is not a full-atlas or all-cohort validation.\n")
