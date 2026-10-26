# Run from the repository root with AnnR installed:
# Rscript data-raw/make_toy.R
# Synthetic teaching data; no patient records or downloaded source data.
suppressPackageStartupMessages(library(AnnR))
set.seed(20261026)
n_cells <- 120L
n_genes <- 40L
cells <- sprintf("cell%03d", seq_len(n_cells))
genes <- sprintf("gene%03d", seq_len(n_genes))
group <- rep(c("groupA", "groupB", "groupC"), each = 40L)
batch <- rep(c("batch1", "batch2"), length.out = n_cells)
counts <- matrix(rpois(n_cells * n_genes, lambda = 0.6), nrow = n_cells,
                 dimnames = list(cells, genes))
for (g in seq_len(3L)) {
  rows <- which(group == c("groupA", "groupB", "groupC")[g])
  cols <- seq.int(1L + (g - 1L) * 5L, g * 5L)
  counts[rows, cols] <- counts[rows, cols] +
    matrix(rpois(length(rows) * length(cols), 4), nrow = length(rows))
}
normalized <- log1p(counts / rowSums(counts) * 1e4)
pc <- prcomp(normalized, center = TRUE, scale. = FALSE)
# Deterministic ring graph for illustrating obsp; this is not an inferred graph.
i <- seq_len(n_cells)
j <- c(seq.int(2L, n_cells), 1L)
graph <- Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1,
                             dims = c(n_cells, n_cells))
obs <- data.frame(group = factor(group), batch = factor(batch),
                  total_counts = as.integer(rowSums(counts)),
                  n_genes = as.integer(rowSums(counts > 0)), row.names = cells)
var <- data.frame(symbol = genes, is_marker = seq_len(n_genes) <= 15L,
                  row.names = genes)
ad <- AnnR(X = Matrix::Matrix(normalized, sparse = TRUE), obs = obs, var = var,
           layers = list(counts = Matrix::Matrix(counts, sparse = TRUE)),
           obsm = list(X_pca = unname(pc$x[, 1:2])),
           varm = list(PCs = unname(pc$rotation[, 1:2])),
           obsp = list(connectivities = graph),
           raw = list(X = Matrix::Matrix(counts, sparse = TRUE), var = var, varm = list()),
           uns = list(title = "AnnR synthetic tutorial dataset",
                      seed = 20261026L, license = "MIT",
                      author = "MD BABU MIA",
                      normalization = "log1p(counts / total_counts * 10000)",
                      graph = "undirected ring; demonstration only"))
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
write_h5ad(ad, "inst/extdata/annr_toy.h5ad", compression = "gzip")
write.csv(counts, "inst/extdata/toy_counts.csv")
write.csv(obs, "inst/extdata/toy_obs.csv")
write.csv(var, "inst/extdata/toy_var.csv")
cat("Created synthetic dataset:", n_cells, "cells x", n_genes, "genes\n")
