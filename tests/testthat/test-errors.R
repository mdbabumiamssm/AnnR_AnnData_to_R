test_that("bad inputs give clear errors", {
  expect_error(read_h5ad(tempfile()), "does not exist")
  bad <- tempfile(fileext = ".h5ad")
  writeBin(as.raw(1:200), bad)
  expect_error(read_h5ad(bad), "HDF5")
  expect_error(read_h5ad(fixture("modern.h5ad"), layers = 3), "TRUE, FALSE")
  expect_error(write_h5ad(AnnR(X = matrix(1, 1, 1)), file.path(tempfile(), "x.h5ad")), "directory")
  expect_error(AnnR(X = matrix(1, 2, 2), obs = data.frame(a = 1:3)), "obs/var")
})
