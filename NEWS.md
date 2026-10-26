# AnnR 0.1.0

* Initial release: native HDF5 reader and writer, selective element reads,
  sparse and dense matrices, legacy encodings, and dtype tracking.
* Conversion to and from SingleCellExperiment and Seurat v5.
* Parallel decoding and writing with gzip, LZF and optional Zstandard.
* Bundled synthetic single-cell example, reproducible generator, and tutorial.
* Reject malformed pairwise matrices and inconsistent raw annotations.
* Preserve existing output on failed file replacement.
* Validation drivers now exit unsuccessfully when comparisons fail.
* Reject changed Seurat identities to prevent reuse of stale original name maps.
