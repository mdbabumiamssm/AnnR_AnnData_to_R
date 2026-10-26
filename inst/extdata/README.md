# AnnR synthetic example data

Author: **MD BABU MIA**. License: MIT (same as AnnR).

`annr_toy.h5ad` contains 120 artificial cells and 40 artificial genes.
No human participants, patient data, or external datasets are involved.
It is suitable for learning file I/O, not for biological conclusions or
performance benchmarking.

* `X`: sparse log1p counts normalized to 10,000 counts per cell, cells x genes.
* `layers/counts`: sparse non-negative simulated counts.
* `obs`: group (three artificial groups), batch (two artificial batches),
  total_counts, and n_genes (number with nonzero counts).
* `var`: synthetic symbol and is_marker (first 15 genes).
* `obsm/X_pca`: first two principal component scores, cells x components.
* `varm/PCs`: matching feature loadings, genes x components.
* `obsp/connectivities`: symmetric ring graph for demonstrating graph storage.
* `raw`: count matrix and all 40 feature annotations.
* `uns`: provenance, normalization, license, seed, and graph description.

`toy_counts.csv`, `toy_obs.csv`, `toy_var.csv` provide human-readable inputs;
the first CSV column contains row names. Counts are cells x genes.

Regenerate from the checkout root with `Rscript data-raw/make_toy.R`.
Seed: 20261026; generated with R 4.6.1. PCA signs or compression bytes may
vary across numerical libraries. Check matrix values and metadata, rather
than expecting identical HDF5 bytes after regeneration.
