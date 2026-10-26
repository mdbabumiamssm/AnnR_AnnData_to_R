# Read an .h5ad with AnnR and write it back.
# Usage: Rscript r_roundtrip.R <in.h5ad> <out.h5ad> [via=annr|sce|seurat] [compression=none]
suppressPackageStartupMessages(library(AnnR))
a <- commandArgs(trailingOnly = TRUE)
via <- if (length(a) >= 3) a[[3]] else "annr"
comp <- if (length(a) >= 4) a[[4]] else "none"
t0 <- proc.time()[["elapsed"]]
obj <- suppressWarnings(switch(via,
  annr = read_h5ad(a[[1]]),
  sce = read_h5ad(a[[1]], as = "SingleCellExperiment"),
  seurat = read_h5ad(a[[1]], as = "Seurat")
))
t1 <- proc.time()[["elapsed"]]
write_h5ad(obj, a[[2]], compression = comp)
t2 <- proc.time()[["elapsed"]]
cat(sprintf("read %.2fs write %.2fs\n", t1 - t0, t2 - t1))
