# Results must not depend on the number of threads (detects races / ordering bugs).
suppressPackageStartupMessages(library(AnnR))
files <- c(readLines("dev/validation/corpus_small.txt"), readLines("dev/work/syn/list.txt"),
           character())
bad <- 0
out1 <- tempfile(fileext = ".h5ad"); out2 <- tempfile(fileext = ".h5ad")
for (f in files) {
  ref <- suppressWarnings(read_h5ad(f, threads = 1))
  for (th in c(2, 3, 7, 16, 32)) {
    if (!identical(ref, suppressWarnings(read_h5ad(f, threads = th)))) { bad <- bad + 1; cat("READ MISMATCH", th, f, "\n") }
  }
  write_h5ad(ref, out1, compression = "gzip", threads = 1)
  write_h5ad(ref, out2, compression = "gzip", threads = 32)
  a <- suppressWarnings(read_h5ad(out1, threads = 4)); b <- suppressWarnings(read_h5ad(out2, threads = 4))
  attr(a, "h5ad_types") <- NULL; attr(b, "h5ad_types") <- NULL
  if (!identical(a, b)) { bad <- bad + 1; cat("WRITE MISMATCH", f, "\n") }
}
cat(sprintf("determinism: %d files x 6 thread counts, mismatches = %d\n", length(files), bad))
