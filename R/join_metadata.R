#' Attach a separately retrieved cell-level annotation without changing row count
#'
#' Useful for PMDBS MapMyCells assignments or MERFISH CCF coordinates. Long-form
#' region maps and Patch-seq sweep features are not cell-level and must remain separate.
#' @param data An abc_data object returned by load_data.
#' @param annotation A data frame from get_metadata.
#' @param by The common, exact join key.
#' @param prefix Prefix applied to all annotation columns except the key.
#' @return An abc_data object with a validated many-to-one left join.
#' @export
join_metadata <- function(data, annotation, by = "cell_label", prefix = "annotation_") {
  stopifnot(inherits(data, "abc_data"), is.data.frame(annotation))
  if (!by %in% names(data$cell_metadata) || !by %in% names(annotation)) stop("Join key missing.", call. = FALSE)
  if (!is.character(data$cell_metadata[[by]]) || !is.character(annotation[[by]])) stop("Join keys must be character, not rounded numeric IDs.", call. = FALSE)
  if (anyNA(annotation[[by]]) || anyDuplicated(annotation[[by]])) stop("Annotation key is missing or duplicated; a cell join would be unsafe.", call. = FALSE)
  annotation_manifest <- attr(annotation, "manifest")
  normalize_manifest <- function(x) sub("^releases/([0-9]+)/manifest.json$", "\\1", x)
  if (!is.null(annotation_manifest) && !is.null(data$release) &&
      !identical(normalize_manifest(annotation_manifest), normalize_manifest(data$release)))
    stop("Annotation and cell metadata use different manifests. Reload the annotation for the same release.", call. = FALSE)
  cols <- setdiff(names(annotation), by)
  new <- paste0(prefix, cols)
  if (any(new %in% names(data$cell_metadata)) || anyDuplicated(new)) stop("Annotation names collide; choose another prefix.", call. = FALSE)
  idx <- match(data$cell_metadata[[by]], annotation[[by]])
  for (i in seq_along(cols)) data$cell_metadata[[new[[i]]]] <- annotation[[cols[[i]]]][idx]
  data$provenance$supplemental_joins <- c(data$provenance$supplemental_joins,
      list(list(by = by, prefix = prefix, manifest = attr(annotation, "manifest"),
                matched_cells = sum(!is.na(idx)), unmatched_cells = sum(is.na(idx)))))
  data
}
