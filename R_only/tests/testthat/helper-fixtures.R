# Offline fixtures: a miniature ABC Atlas cache (manifest + CSVs + h5ad) laid
# out exactly like a real download directory, so the loaders can run end to end
# with offline = TRUE.

# Write a minimal h5ad holding `counts` (cells x genes) as CSR, in the on-disk
# layout used by the ABC Atlas files (anndata >= 0.8, `_index` attributes).
write_h5ad <- function(path, counts, cell_labels, gene_ids, gene_symbols) {
  csr <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "RsparseMatrix")
  rhdf5::h5createFile(path)
  fid <- rhdf5::H5Fopen(path)
  on.exit(rhdf5::H5Fclose(fid))

  set_attrs <- function(group, attrs) {
    gid <- rhdf5::H5Gopen(fid, group)
    on.exit(rhdf5::H5Gclose(gid))
    for (n in names(attrs)) rhdf5::h5writeAttribute(attrs[[n]], gid, n)
  }

  rhdf5::h5createGroup(fid, "X")
  rhdf5::h5write(as.numeric(csr@x), fid, "X/data")
  rhdf5::h5write(as.integer(csr@j), fid, "X/indices")
  rhdf5::h5write(as.integer(csr@p), fid, "X/indptr")
  set_attrs("X", list(`encoding-type` = "csr_matrix", `encoding-version` = "0.1.0",
                      shape = as.integer(dim(counts))))

  rhdf5::h5createGroup(fid, "obs")
  rhdf5::h5write(cell_labels, fid, "obs/cell_label")
  set_attrs("obs", list(`_index` = "cell_label", `encoding-type` = "dataframe"))

  rhdf5::h5createGroup(fid, "var")
  rhdf5::h5write(gene_ids, fid, "var/gene_identifier")
  rhdf5::h5write(gene_symbols, fid, "var/gene_symbol")
  set_attrs("var", list(`_index` = "gene_identifier", `encoding-type` = "dataframe"))
  invisible(path)
}

# Build a fixture cache in `dir` for a WHB-like dataset "TEST-10X" with two
# expression files, a taxonomy, a region-of-interest map and a gene table.
# `var_ids` overrides the gene identifiers stored in the h5ad files.
# Returns a list with the directory and the true dense counts (genes x cells).
make_fixture_cache <- function(dir = withr::local_tempdir(.local_envir = parent.frame()),
                               var_ids = NULL) {
  version <- "20250101"
  listing <- list()
  add_file <- function(directory, kind, name, rel, sub = NULL) {
    path <- file.path(dir, rel)
    entry <- list(version = version, relative_path = rel,
                  url = paste0("https://example.invalid/", rel),
                  size = file.size(path), file_hash = unname(tools::md5sum(path)))
    files <- list(files = stats::setNames(list(entry), tools::file_ext(rel)))
    if (is.null(sub)) {
      listing[[directory]][[kind]][[name]] <<- files
    } else {
      listing[[directory]][[kind]][[name]][[sub]] <<- files
    }
  }
  write_csv <- function(df, directory, name) {
    rel <- file.path("metadata", directory, version, paste0(name, ".csv"))
    dir.create(file.path(dir, dirname(rel)), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, file.path(dir, rel), row.names = FALSE, na = "")
    add_file(directory, "metadata", name, rel)
  }

  cells <- data.frame(
    cell_label = c("c1", "c2", "c3", "c4", "c5", "c6"),
    feature_matrix_label = c("TEST-A", "TEST-A", "TEST-A", "TEST-B", "TEST-B", "TEST-B"),
    dataset_label = "TEST-10X",
    donor_sex = c("M", "F", "M", "F", NA, "M"),
    region_of_interest_label = c("R1", "R2", "R1", "R2", "R1", "R1"),
    cluster_alias = c(1L, 2L, 1L, 2L, 2L, 1L)
  )
  write_csv(cells, "TEST-10X", "cell_metadata")

  genes <- data.frame(
    gene_identifier = c("G1", "G2", "G3", "G4"),
    gene_symbol = c("Aaa", "Bbb", "Aaa", "Ddd"),   # duplicated symbol on purpose
    name = c("alpha ", "beta", "gamma", "delta")   # trailing space must survive
  )
  write_csv(genes, "TEST-10X", "gene")

  write_csv(data.frame(region_of_interest_label = c("R1", "R2"),
                       color_hex_triplet = c("#111111", "#222222")),
            "TEST-10X", "region_of_interest_structure_map")

  write_csv(data.frame(
    cluster_alias = c(1L, 1L, 2L, 2L),
    cluster_annotation_term_set_name = c("cluster", "supercluster", "cluster", "supercluster"),
    cluster_annotation_term_name = c("cl1", "sup1", "cl2", NA),
    color_hex_triplet = c("#aa0000", "#bb0000", "#cc0000", "#dd0000")
  ), "TEST-taxonomy", "cluster_to_cluster_annotation_membership")
  write_csv(data.frame(label = c("L1", "L2"), name = c("supercluster", "cluster")),
            "TEST-taxonomy", "cluster_annotation_term_set")

  # Expression: genes x cells truth, split across two files by cell.
  truth <- matrix(c(1, 0, 3, 0,
                    0, 2, 0, 0,
                    5, 0, 0, 1,
                    0, 0, 7, 0,
                    0, 4, 0, 0,
                    0, 0, 0, 0), nrow = 4,
                  dimnames = list(genes$gene_identifier, cells$cell_label))
  for (f in c("TEST-A", "TEST-B")) {
    in_file <- cells$cell_label[cells$feature_matrix_label == f]
    in_file <- rev(in_file)  # on-disk order differs from metadata order
    rel <- file.path("expression_matrices", "TEST-10X", version, paste0(f, "-raw.h5ad"))
    dir.create(file.path(dir, dirname(rel)), recursive = TRUE, showWarnings = FALSE)
    write_h5ad(file.path(dir, rel), t(truth[, in_file]), in_file,
               var_ids %||% genes$gene_identifier, genes$gene_symbol)
    add_file("TEST-10X", "expression_matrices", f, rel, sub = "raw")
  }

  manifest <- list(version = version, resource_uri = "s3://example.invalid/",
                   directory_listing = list(), file_listing = listing)
  dir.create(file.path(dir, "releases", version), recursive = TRUE)
  jsonlite::write_json(manifest, file.path(dir, "releases", version, "manifest.json"),
                       auto_unbox = TRUE)

  list(dir = dir, truth = truth, cells = cells, genes = genes)
}

# Dataset spec for the fixture, in the format of the .abc_datasets registry.
fixture_spec <- list(
  cell = list(dir = "TEST-10X", col_classes = c(cell_label = "character")),
  taxonomy = "TEST-taxonomy",
  roi = list(dir = "TEST-10X", file = "region_of_interest_structure_map"),
  gene = list(dir = "TEST-10X"),
  cols_of_interest = c("donor_sex", "cluster", "supercluster", "not_a_column")
)
