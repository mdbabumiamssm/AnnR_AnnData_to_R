#' Convert objects to AnnR
#'
#' Converts `SingleCellExperiment`, `SummarizedExperiment` and `Seurat`
#' objects to [AnnR] objects (cells x genes). Sparse `dgCMatrix` assays are
#' reinterpreted as `dgRMatrix` without copying.
#'
#' @param x Object to convert.
#' @param ... Method-specific arguments.
#' @return An [AnnR] object.
#' @export
as_AnnR <- function(x, ...) UseMethod("as_AnnR")

#' @export
as_AnnR.AnnR <- function(x, ...) x

#' @export
as_AnnR.default <- function(x, ...) {
  stop("cannot convert an object of class ", paste(class(x), collapse = "/"), " to AnnR", call. = FALSE)
}

.require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("package '", pkg, "' is required for this conversion", call. = FALSE)
  }
}

# Reinterpret a compressed sparse matrix as its transpose: CSR <-> CSC share
# the same three arrays, so this is O(1) and allocates nothing.
.t_free <- function(m) {
  cls <- class(m)[1L]
  to <- switch(cls,
    dgRMatrix = "dgCMatrix", dgCMatrix = "dgRMatrix",
    lgRMatrix = "lgCMatrix", lgCMatrix = "lgRMatrix",
    ngRMatrix = "ngCMatrix", ngCMatrix = "ngRMatrix",
    NULL
  )
  if (is.null(to)) return(Matrix::t(m))
  out <- methods::new(to)
  if (endsWith(cls, "RMatrix")) out@i <- m@j else out@j <- m@i
  out@p <- m@p
  if (!startsWith(cls, "n")) out@x <- m@x
  out@Dim <- rev(m@Dim)
  out@Dimnames <- rev(m@Dimnames)
  out
}

# AnnR (obs x var) -> features x cells for Bioconductor / Seurat.
.to_feature_major <- function(m, already = FALSE) {
  if (is.null(m)) return(NULL)
  if (methods::is(m, "sparseMatrix")) {
    out <- .t_free(m)
    if (methods::is(out, "RsparseMatrix")) out <- methods::as(out, "CsparseMatrix")
    return(out)
  }
  if (already) return(m)
  if (methods::is(m, "DelayedArray")) return(BiocGenerics::t(m))
  .t_dense(m)
}

# t() that also transposes the exact 64-bit integer patterns kept by the reader.
.t_dense <- function(m) {
  bits <- attr(m, "h5ad_int64")
  out <- t(m)
  if (!is.null(bits) && length(bits) == length(m)) {
    attr(out, "h5ad_int64") <- as.vector(t(matrix(bits, nrow(m), ncol(m))))
  }
  out
}

# features x cells -> AnnR (obs x var).
.to_obs_major <- function(m) {
  if (is.null(m)) return(NULL)
  if (methods::is(m, "DelayedArray")) m <- .realize_delayed(m)
  if (methods::is(m, "sparseMatrix")) {
    if (!methods::is(m, "CsparseMatrix") && !methods::is(m, "RsparseMatrix")) {
      m <- methods::as(m, "CsparseMatrix")
    }
    if (!methods::is(m, "generalMatrix")) m <- methods::as(m, "generalMatrix")
    return(.t_free(m))
  }
  if (methods::is(m, "Matrix")) m <- as.matrix(m)
  .t_dense(m)
}

.strip_dimnames <- function(m) {
  if (is.null(m)) return(m)
  if (methods::is(m, "Matrix")) {
    m@Dimnames <- list(NULL, NULL)
  } else if (is.matrix(m) || is.array(m)) {
    dimnames(m) <- NULL
  }
  m
}

.df_from_DF <- function(DF, names) {
  df <- as.data.frame(DF, optional = TRUE, stringsAsFactors = FALSE)
  if (ncol(df) == 0L) df <- data.frame(row.names = seq_len(nrow(DF)))
  rownames(df) <- names
  df
}

#' Convert an AnnR object to a SingleCellExperiment
#'
#' `X` becomes the assay named `X_name`, layers become further assays (all
#' features x cells), `obs`/`var` become `colData`/`rowData`, `obsm` becomes
#' `reducedDims`, `varm` is stored as `rowData(sce)$varm`, `obsp`/`varp` become
#' `colPairs`/`rowPairs`, `uns` becomes `metadata()` and `raw` becomes the
#' alternative experiment `"raw"`. Sparse matrices are reinterpreted without
#' copying. The h5ad dtype records are kept in `int_metadata(sce)` so that
#' [write_h5ad()] can restore them.
#'
#' @param x An [AnnR] object.
#' @param X_name Assay name for `X`.
#' @param ... Unused.
#' @return A `SingleCellExperiment`.
#' @export
as_SingleCellExperiment <- function(x, X_name = "X", ...) {
  .require("SingleCellExperiment")
  .require("SummarizedExperiment")
  .require("S4Vectors")
  fm <- isTRUE(attr(x, "annr_feature_major"))
  x <- unclass(x)
  cells <- rownames(x$obs)
  genes <- rownames(x$var)
  assays <- list()
  if (!is.null(x$X)) assays[[X_name]] <- .to_feature_major(x$X, fm)
  for (k in names(x$layers)) {
    nm <- if (k %in% names(assays)) paste0(k, "_layer") else k
    assays[[nm]] <- .to_feature_major(x$layers[[k]], fm)
  }
  assays <- lapply(assays, .strip_dimnames)
  col_data <- S4Vectors::DataFrame(x$obs, check.names = FALSE)
  row_data <- S4Vectors::DataFrame(x$var, check.names = FALSE)
  if (length(x$varm)) {
    vm <- S4Vectors::make_zero_col_DFrame(length(genes))
    for (k in names(x$varm)) vm[[k]] <- .strip_dimnames(x$varm[[k]])
    row_data$varm <- vm
  }
  red <- lapply(x$obsm, function(m) {
    if (is.data.frame(m)) m <- S4Vectors::DataFrame(m, check.names = FALSE)
    m <- .strip_dimnames(m)
    if (is.matrix(m)) rownames(m) <- cells
    m
  })
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays, colData = col_data, rowData = row_data,
    reducedDims = red, metadata = x$uns
  )
  rownames(sce) <- genes
  colnames(sce) <- cells
  for (k in names(x$obsp)) {
    SingleCellExperiment::colPair(sce, k) <- .as_dgC(x$obsp[[k]])
  }
  for (k in names(x$varp)) {
    SingleCellExperiment::rowPair(sce, k) <- .as_dgC(x$varp[[k]])
  }
  if (!is.null(x$raw) && !is.null(x$raw$X)) {
    rv <- x$raw$var %||% data.frame(row.names = .default_names(ncol(x$raw$X)))
    rd <- S4Vectors::DataFrame(rv, check.names = FALSE)
    if (length(x$raw$varm)) {
      vm <- S4Vectors::make_zero_col_DFrame(nrow(rv))
      for (k in names(x$raw$varm)) vm[[k]] <- .strip_dimnames(x$raw$varm[[k]])
      rd$varm <- vm
    }
    raw_sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(X = .strip_dimnames(.to_feature_major(x$raw$X, fm))), rowData = rd
    )
    rownames(raw_sce) <- rownames(rv)
    colnames(raw_sce) <- cells
    SingleCellExperiment::altExp(sce, "raw") <- raw_sce
  }
  types <- attr(structure(x, class = "AnnR"), "h5ad_types") %||% attr(x, "h5ad_types")
  im <- SingleCellExperiment::int_metadata(sce)
  im$AnnR_types <- types
  im$AnnR_X_name <- X_name
  SingleCellExperiment::int_metadata(sce) <- im
  sce
}

.as_dgC <- function(m) {
  if (methods::is(m, "sparseMatrix")) return(methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix"))
  methods::as(methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix"), "generalMatrix")
}

#' @rdname as_AnnR
#' @param X_name For `SingleCellExperiment`/`SummarizedExperiment`: the assay
#'   used as `X` (default: `"X"` when present, else the first assay). For
#'   `Seurat`: see [as_AnnR.Seurat()].
#' @export
as_AnnR.SummarizedExperiment <- function(x, X_name = NULL, ...) {
  .require("SummarizedExperiment")
  .require("S4Vectors")
  an <- SummarizedExperiment::assayNames(x)
  if (is.null(an) && length(SummarizedExperiment::assays(x))) {
    an <- paste0("assay", seq_along(SummarizedExperiment::assays(x)))
  }
  is_sce <- methods::is(x, "SingleCellExperiment")
  im <- if (is_sce) SingleCellExperiment::int_metadata(x) else list()
  if (is.null(X_name)) {
    X_name <- if (!is.null(im$AnnR_X_name) && im$AnnR_X_name %in% an) im$AnnR_X_name
              else if ("X" %in% an) "X" else if (length(an)) an[[1L]] else NULL
  }
  cells <- colnames(x) %||% .default_names(ncol(x))
  genes <- rownames(x) %||% .default_names(nrow(x))
  X <- if (!is.null(X_name)) .to_obs_major(SummarizedExperiment::assay(x, X_name, withDimnames = FALSE)) else NULL
  layers <- list()
  for (k in setdiff(an, X_name)) {
    layers[[k]] <- .to_obs_major(SummarizedExperiment::assay(x, k, withDimnames = FALSE))
  }
  rd <- SummarizedExperiment::rowData(x)
  varm <- list()
  if ("varm" %in% colnames(rd) && methods::is(rd$varm, "DataFrame")) {
    vm <- rd$varm
    for (k in colnames(vm)) varm[[k]] <- vm[[k]]
    rd$varm <- NULL
  }
  obs <- .df_from_DF(SummarizedExperiment::colData(x), cells)
  var <- .df_from_DF(rd, genes)
  obsm <- list()
  obsp <- list()
  varp <- list()
  raw <- NULL
  if (is_sce) {
    for (k in SingleCellExperiment::reducedDimNames(x)) {
      m <- SingleCellExperiment::reducedDim(x, k, withDimnames = FALSE)
      if (methods::is(m, "DataFrame")) {
        m <- as.data.frame(m, optional = TRUE)
        rownames(m) <- cells
      }
      obsm[[k]] <- .strip_dimnames(m)
    }
    for (k in SingleCellExperiment::colPairNames(x)) {
      obsp[[k]] <- SingleCellExperiment::colPair(x, k, asSparse = TRUE)
    }
    for (k in SingleCellExperiment::rowPairNames(x)) {
      varp[[k]] <- SingleCellExperiment::rowPair(x, k, asSparse = TRUE)
    }
    if ("raw" %in% SingleCellExperiment::altExpNames(x)) {
      ra <- SingleCellExperiment::altExp(x, "raw")
      rraw <- as_AnnR.SummarizedExperiment(ra)
      raw <- list(X = rraw$X, var = rraw$var, varm = rraw$varm)
    }
  }
  obsp <- lapply(obsp, .strip_dimnames)
  varp <- lapply(varp, .strip_dimnames)
  out <- AnnR(X = X, obs = obs, var = var, obsm = obsm, varm = varm, obsp = obsp, varp = varp,
              layers = layers, uns = S4Vectors::metadata(x), raw = raw)
  attr(out, "h5ad_types") <- im$AnnR_types
  out
}
