test_that("malformed pairwise matrices and raw annotations are rejected", {
  m <- matrix(1, 2, 3)
  expect_error(AnnR(X = m, obsp = list(graph = 1:2)), "must be 2 x 2")
  expect_error(AnnR(X = m, varp = list(graph = array(1, c(3, 3, 3)))), "must be 3 x 3")
  expect_error(AnnR(X = m, raw = list(X = 1:2)), "two-dimensional")
  expect_error(AnnR(X = m, raw = list(X = m, var = data.frame(a = 1:2))), "raw variable")
})

test_that("failed replacement preserves the existing destination", {
  f <- tmp_h5ad()
  on.exit(unlink(f))
  original <- AnnR(X = matrix(1, 2, 3))
  write_h5ad(original, f)
  before <- readBin(f, "raw", n = file.info(f)$size)
  expect_error(write_h5ad(original, f, overwrite = FALSE), "file exists")
  broken <- original
  broken$uns <- list("invalid/key" = 1)
  expect_error(write_h5ad(broken, f), "may not contain")
  expect_identical(readBin(f, "raw", n = file.info(f)$size), before)
  replacement <- AnnR(X = matrix(2, 2, 3))
  write_h5ad(replacement, f)
  expect_equal(as.matrix(read_h5ad(f)$X), replacement$X)
})

test_that("failed final rename preserves the original file", {
  skip_on_os("windows")
  f <- tmp_h5ad()
  on.exit(unlink(f))
  ad <- AnnR(X = matrix(1, 2, 3))
  write_h5ad(ad, f)
  before <- readBin(f, "raw", n = file.info(f)$size)
  fail_rename <- write_h5ad
  environment(fail_rename) <- list2env(list(file.rename = function(...) FALSE),
                                    parent = environment(write_h5ad))
  expect_error(fail_rename(ad, f), "existing file preserved")
  expect_identical(readBin(f, "raw", n = file.info(f)$size), before)
})

test_that("the bundled synthetic dataset is complete", {
  f <- system.file("extdata", "annr_toy.h5ad", package = "AnnR", mustWork = TRUE)
  ad <- read_h5ad(f)
  expect_identical(dim(ad), c(120L, 40L))
  expect_identical(dim(ad$obsm$X_pca), c(120L, 2L))
  expect_equal(as.numeric(Matrix::rowSums(ad$layers$counts)), as.numeric(ad$obs$total_counts))
  expect_equal(as.matrix(ad$X), log1p(as.matrix(ad$layers$counts) / ad$obs$total_counts * 1e4))
  expect_identical(as.character(ad$uns$author), "MD BABU MIA")
})
