# Review scripts/benchmark_config.json first. Nothing runs while enabled=false.
# Called once outside timed runs; metadata downloads here are not benchmark results.
library(importABCatlas)
args <- commandArgs(trailingOnly=TRUE)
config <- if(length(args)) args[[1]] else "scripts/benchmark_config.json"
x <- jsonlite::read_json(config,simplifyVector=FALSE)
profiles <- list()
for (p in x$profiles) {
  if (!isTRUE(p$enabled)) next
  if (!isTRUE(p$cohort_reviewed)) stop("Review the exact neuronal/non-neuronal definition for ",p$dataset)
  if (!length(p$filters) && !(p$dataset=="HMBA_PatchSeq" && p$cohort=="neuronal"))
    stop("An exact biological cohort filter is required for ", p$dataset, "/", p$cohort)
  d <- load_data(p$dataset, download_base=x$selection_cache, release=x$release,
                 filters=p$filters, n_cells=as.integer(p$n_cells), seed=as.integer(x$seed))
  genes <- unlist(p$genes,use.names=FALSE)
  if (!length(genes)) {
    genes <- if(p$dataset=="MERFISH_imputed") c("Calb2","Baiap3","Lypd1") else head(d$gene_data$gene_identifier,3L)
  }
  profiles[[length(profiles)+1L]] <- list(dataset=p$dataset,cohort=p$cohort,
      cell_ids=as.list(d$cell_metadata$cell_label),genes=as.list(genes),release=d$release,
      data_type=p$data_type,return_type=p$return_type,cohort_definition=p$filters,
      selection_seed=x$seed,cohort_reviewed=TRUE,annotations=TRUE)
}
if (!length(profiles)) stop("No enabled, reviewed profiles. Read docs/BENCHMARK_PROTOCOL.md.")
jsonlite::write_json(list(profiles=profiles),"benchmark_profiles.json",auto_unbox=TRUE,pretty=TRUE,null="null")
cat("Frozen selections written to benchmark_profiles.json. These contain no timings.\n")
