# Dataset convenience wrappers retain the original names; all share load_data().

#' Load WMB metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_WMB <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("WMB", download_base = download_base, ...)
}

#' Load WHB metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_WHB <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("WHB", download_base = download_base, ...)
}

#' Load AgingMouse metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_AgingMouse <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("AgingMouse", download_base = download_base, ...)
}

#' Load PMDBS metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_PMDBS <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("PMDBS", download_base = download_base, ...)
}

#' Load HMBA metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_HMBA <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("HMBA", download_base = download_base, ...)
}

#' Load DevelopingMouse metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_DevelopingMouse <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("DevelopingMouse", download_base = download_base, ...)
}

#' Load SEAAD_Multiregion metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_SEAAD_Multiregion <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("SEAAD_Multiregion", download_base = download_base, ...)
}

#' Load SEAAD_CaH metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_SEAAD_CaH <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("SEAAD_CaH", download_base = download_base, ...)
}

#' Load HMBA_PatchSeq metadata
#' @param download_base Cache directory outside the repository.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_HMBA_PatchSeq <- function(download_base = "~/abc_atlas_cache", ...) {
  load_data("HMBA_PatchSeq", download_base = download_base, ...)
}

#' Load measured or imputed Allen or Zhuang MERFISH metadata
#' @param source Allen or Zhuang1 through Zhuang4.
#' @param imputed Request the Allen imputed product, not measured counts.
#' @param ... Further arguments to load_data.
#' @return An abc_data list. Fetch imputed expression with data_type='log2'.
#' @export
load_data_MERFISH <- function(source = c("Allen", "Zhuang1", "Zhuang2", "Zhuang3", "Zhuang4"),
                              imputed = FALSE, ...) {
  source <- match.arg(source)
  if (imputed && source != "Allen") stop("The imputed route is defined only for the Allen product.", call. = FALSE)
  id <- if (source != "Allen") source else if (imputed) "MERFISH_imputed" else "MERFISH"
  load_data(id, ...)
}

#' Load a consensus mouse component without conflating gene universes
#' @param source AIBS or Macosko.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_ConsensusMouse <- function(source = c("AIBS", "Macosko"), ...) {
  load_data(paste0("ConsensusMouse", match.arg(source)), ...)
}

#' Load an HMBA native-species expression route
#' @param species Human, Macaque or Marmoset.
#' @param ... Further arguments to load_data; exact file overrides may be needed.
#' @return An abc_data list. Species-specific files are discovered, not invented.
#' @export
load_data_HMBA_species <- function(species = c("Human", "Macaque", "Marmoset"), ...) {
  load_data(paste0("HMBA_", match.arg(species)), ...)
}

#' Load HMBA MERSCOPE or Xenium spatial metadata
#' @param species Human, Macaque or Marmoset.
#' @param ... Further arguments to load_data.
#' @return An abc_data list.
#' @export
load_data_HMBA_spatial <- function(species = c("Human", "Macaque", "Marmoset"), ...) {
  load_data(paste0("HMBA_spatial_", match.arg(species)), ...)
}
