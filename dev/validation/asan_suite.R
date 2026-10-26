# Exercise read/write/convert paths; run under AddressSanitizer/UBSan.
suppressPackageStartupMessages(library(AnnR))
files <- c(readLines("dev/validation/corpus_small.txt"), readLines("dev/work/syn/list.txt"))
out <- tempfile(fileext = ".h5ad")
n <- 0
for (f in files) {
  ad <- suppressWarnings(read_h5ad(f, threads = 8))
  for (comp in c("none", "gzip", "lzf", "zstd")) {
    write_h5ad(ad, out, compression = comp, threads = 8)
    back <- suppressWarnings(read_h5ad(out, threads = 8))
    n <- n + 1
  }
  sce <- suppressWarnings(read_h5ad(f, as = "SingleCellExperiment", threads = 8))
  write_h5ad(sce, out)
  seu <- suppressWarnings(read_h5ad(f, as = "Seurat"))
  write_h5ad(seu, out)
  sub <- ad[seq_len(min(5, nrow(ad))), ]
  write_h5ad(sub, out)
}
cat("ASAN SUITE OK:", length(files), "files,", n, "write/read cycles\n")
