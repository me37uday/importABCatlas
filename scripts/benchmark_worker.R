# Internal fresh-process worker; run through benchmark.py.
library(importABCatlas)
args <- commandArgs(trailingOnly=TRUE)
req <- jsonlite::read_json(args[[1]],simplifyVector=FALSE)
p <- req$profile
stage <- req$stage
out <- req$out_dir
dir.create(out,recursive=TRUE,showWarnings=FALSE)
flag <- file.path(out,"stage.txt")
error <- NULL
tryCatch({
  if(stage=="fetch_data") d <- readRDS(req$metadata_rds)
  gc()
  writeLines("running",flag)
  start <- proc.time()[["elapsed"]]
  if(stage=="load_data") {
    d <- load_data(p$dataset,download_base=req$cache,release=p$release,
                   cell_ids=unlist(p$cell_ids,use.names=FALSE),annotations=TRUE)
    seconds <- proc.time()[["elapsed"]]-start
    writeLines("done",flag)
    saveRDS(d,req$metadata_rds) # serialization is outside timed function
    dimensions <- c(nrow(d$gene_data),nrow(d$cell_metadata))
    provenance <- d$provenance
  } else {
    value <- fetch_data(data=d,genes=unlist(p$genes,use.names=FALSE),data_type=p$data_type,
                        return_type=p$return_type,allow_downloads=TRUE,max_output_mb=256)
    seconds <- proc.time()[["elapsed"]]-start
    writeLines("done",flag)
    dimensions <- if(p$return_type=="matrix") dim(value$expression) else dim(value)
    provenance <- if(p$return_type=="matrix") value$provenance else if(p$return_type=="seurat") value@misc$importABCatlas else
      S4Vectors::metadata(value)$importABCatlas
  }
  jsonlite::write_json(list(status="PASS",seconds=seconds,n_genes=dimensions[[1]],n_cells=dimensions[[2]],
      provenance=provenance,session=utils::capture.output(sessionInfo())),file.path(out,"result.json"),auto_unbox=TRUE,pretty=TRUE)
},error=function(e) {
  writeLines("failed",flag)
  jsonlite::write_json(list(status="FAIL",error=conditionMessage(e)),file.path(out,"result.json"),auto_unbox=TRUE,pretty=TRUE)
  error <<- e
})
if(!is.null(error)) quit(status=1L)
