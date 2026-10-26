test_that("SingleCellExperiment conversion is lossless", {
  skip_if_not_installed("SingleCellExperiment")
  ad <- read_fixture(fixture("modern.h5ad"))
  sce <- read_fixture(fixture("modern.h5ad"), as = "SingleCellExperiment")
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), rev(dim(ad)))
  expect_s4_class(SummarizedExperiment::assay(sce, "X"), "dgCMatrix")
  expect_equal(as.matrix(SummarizedExperiment::assay(sce, "X")), t(as.matrix(ad$X)), ignore_attr = TRUE)
  expect_identical(colnames(sce), obs_names(ad))
  expect_true("X_pca" %in% SingleCellExperiment::reducedDimNames(sce))
  expect_true("raw" %in% SingleCellExperiment::altExpNames(sce))
  f <- tmp_h5ad()
  write_h5ad(sce, f)
  back <- read_fixture(f)
  expect_identical(back$obs, ad$obs)
  expect_equal(as.matrix(back$X), as.matrix(ad$X))
  expect_identical(class(back$X), class(ad$X))
  expect_identical(class(back$layers$counts), class(ad$layers$counts))
})

test_that("Seurat conversion is lossless", {
  skip_if_not_installed("SeuratObject")
  ad <- read_fixture(fixture("modern.h5ad"))
  seu <- as_Seurat(ad)
  expect_s4_class(seu, "Seurat")
  expect_equal(ncol(seu), nrow(ad))
  f <- tmp_h5ad()
  write_h5ad(seu, f)
  back <- read_fixture(f)
  expect_identical(obs_names(back), obs_names(ad))
  expect_identical(var_names(back), var_names(ad))
  expect_identical(back$obs, ad$obs)
  expect_equal(as.matrix(back$X), as.matrix(ad$X))
  expect_equal(as.matrix(back$layers$counts), as.matrix(ad$layers$counts))
  expect_equal(back$obsm$X_pca, ad$obsm$X_pca, ignore_attr = TRUE)
})


test_that("changed Seurat identities cannot silently reuse original metadata", {
  skip_if_not_installed("SeuratObject")
  ad <- read_fixture(fixture("modern.h5ad"))
  seu <- as_Seurat(ad)
  renamed <- SeuratObject::RenameCells(seu, new.names = paste0("renamed", colnames(seu)))
  expect_error(as_AnnR(renamed), "identities have changed")
  subsetted <- seu[, seq_len(ncol(seu) - 1L)]
  expect_error(as_AnnR(subsetted), "identities have changed")
  reordered <- seu[, rev(seq_len(ncol(seu)))]
  # Seurat preserves its original cell order for this subsetting operation.
  expect_identical(colnames(reordered), colnames(seu))
  fewer_features <- seu[seq_len(nrow(seu) - 1L), ]
  expect_error(as_AnnR(fewer_features), "identities have changed")
  expect_equal(as.matrix(as_AnnR(seu)$X), as.matrix(ad$X))
})
