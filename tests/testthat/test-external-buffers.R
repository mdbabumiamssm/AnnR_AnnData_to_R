test_that("external buffers preserve copy-on-modify, serialization and lifetime", {
  # Above the 2^16 element ALTREP threshold, but small enough for routine CI.
  m <- matrix(as.double(seq_len(80000)), nrow = 400L)
  f <- tmp_h5ad()
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(c(f, rds)))
  write_h5ad(AnnR(X = m), f, compression = "gzip")
  ad <- read_h5ad(f, threads = 2L)
  expect_equal(ad$X, m)
  copy <- ad$X
  copy[1, 1] <- -100
  expect_equal(ad$X[1, 1], 1)
  expect_equal(copy[1, 1], -100)
  saveRDS(ad, rds)
  restored <- readRDS(rds)
  expect_equal(restored$X, m)
  kept <- ad$X
  rm(ad)
  invisible(gc())
  expect_equal(kept, m)
  write_h5ad(restored, f, compression = "gzip", threads = 1L)
  expect_equal(read_h5ad(f, threads = 4L)$X, m)
})
