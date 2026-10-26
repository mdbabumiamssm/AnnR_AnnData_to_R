#' AnnR: fast and faithful AnnData (.h5ad) input/output for R
#'
#' AnnR reads and writes `.h5ad` files, the on-disk format of the Python
#' `anndata` package, with a native C++ engine. Numeric data are decoded in
#' parallel straight into R memory, sparse matrices map onto `Matrix` classes
#' without copies, and the on-disk dtypes of every element are remembered so
#' that [write_h5ad()] reproduces them.
#'
#' Main functions:
#' * [read_h5ad()] reads a file into an [AnnR] object, a
#'   `SingleCellExperiment` or a `Seurat` object.
#' * [write_h5ad()] writes an [AnnR], `SingleCellExperiment` or `Seurat`
#'   object.
#' * [h5ad_info()] lists the structure of a file without reading data.
#'
#' Options:
#' * `AnnR.threads`: number of threads (default: up to 16 available cores).
#'
#' @keywords internal
#' @useDynLib AnnR, .registration = TRUE
#' @importFrom methods new is as slot
#' @importClassesFrom Matrix dgCMatrix dgRMatrix lgCMatrix lgRMatrix
#' @importFrom Matrix t
"_PACKAGE"
