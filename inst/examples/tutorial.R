# Executable companion to docs/TUTORIAL.md. Run after installing AnnR:
# Rscript -e 'source(system.file("examples", "tutorial.R", package="AnnR"))'
library(AnnR)
file <- system.file("extdata", "annr_toy.h5ad", package = "AnnR", mustWork = TRUE)
ad <- read_h5ad(file)
stopifnot(identical(dim(ad), c(120L, 40L)),
          identical(obs_names(ad), rownames(ad$obs)),
          identical(var_names(ad), rownames(ad$var)))
print(ad)
print(head(h5ad_info(file)))
pca <- read_h5ad_elem(file, "obsm/X_pca")
stopifnot(identical(dim(pca), c(120L, 2L)))
partial <- read_h5ad(file, X = FALSE, layers = "counts", obsm = "X_pca",
                     varm = FALSE, obsp = FALSE, varp = FALSE, raw = FALSE, uns = FALSE)
stopifnot(is.null(partial$X), identical(names(partial$layers), "counts"))
# Subsetting cells also subsets layers, embeddings and graphs; raw keeps genes.
small <- ad[ad$obs$group == "groupA", seq_len(10L)]
stopifnot(identical(dim(small), c(40L, 10L)),
          identical(dim(small$obsp$connectivities), c(40L, 40L)),
          identical(dim(small$raw$X), c(40L, 40L)))
# Use a temporary directory so this example does not overwrite user files.
out_dir <- tempfile("annr_tutorial_")
dir.create(out_dir)
out <- file.path(out_dir, "annr_roundtrip.h5ad")
write_h5ad(ad, out, compression = "gzip", overwrite = FALSE, threads = 2L)
back <- read_h5ad(out, threads = 1L)
stopifnot(isTRUE(all.equal(as.matrix(back$X), as.matrix(ad$X))),
          identical(back$obs, ad$obs))
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  sce <- as_SingleCellExperiment(ad)
  stopifnot(identical(dim(sce), rev(dim(ad))),
            isTRUE(all.equal(unname(as.matrix(SummarizedExperiment::assay(sce, "X"))),
                             unname(t(as.matrix(ad$X))))))
  write_h5ad(sce, file.path(out_dir, "sce_roundtrip.h5ad"), compression = "gzip")
}
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  seu <- as_Seurat(ad, counts = "counts", data = "X")
  stopifnot(ncol(seu) == nrow(ad), nrow(seu) == ncol(ad))
  write_h5ad(seu, file.path(out_dir, "seurat_roundtrip.h5ad"), compression = "gzip")
}
print(annr_build_info())
cat("Tutorial completed. Output files are in:", out_dir, "\n")
