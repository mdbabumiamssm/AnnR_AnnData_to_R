#' AnnR objects
#'
#' An `AnnR` object is a plain named list mirroring Python's `AnnData`, with
#' the same orientation: `X` and every layer are `n_obs x n_vars`
#' (cells x genes). Sparse CSR data are held as `dgRMatrix` and CSC data as
#' `dgCMatrix`, both exactly as stored on disk, so no copy or transposition is
#' made when reading.
#'
#' Elements:
#' * `X`: matrix (dense, `dgRMatrix`/`dgCMatrix`) or `NULL`.
#' * `obs`, `var`: data frames whose row names are the obs/var names.
#' * `obsm`, `varm`, `obsp`, `varp`, `layers`, `uns`: named lists.
#' * `raw`: `NULL` or a list with `X`, `var` and `varm`.
#'
#' The attribute `"h5ad_types"` records the on-disk dtype and encoding of every
#' element read from a file; [write_h5ad()] uses it to reproduce the original
#' dtypes whenever this is lossless.
#'
#' @param X Main matrix (`n_obs x n_vars`) or `NULL`.
#' @param obs,var Data frames of observation / variable annotations.
#' @param obsm,varm,obsp,varp,layers,uns Named lists.
#' @param raw `NULL` or a list with elements `X`, `var` and `varm`.
#' @return An object of class `AnnR`.
#' @examples
#' m <- matrix(rpois(20, 2), nrow = 4,
#'             dimnames = list(paste0("cell", 1:4), paste0("gene", 1:5)))
#' ad <- AnnR(X = m, obs = data.frame(group = c("a", "a", "b", "b")))
#' ad
#' dim(ad)
#' @export
AnnR <- function(X = NULL, obs = NULL, var = NULL, obsm = list(), varm = list(),
                 obsp = list(), varp = list(), layers = list(), uns = list(), raw = NULL) {
  n_obs <- if (!is.null(obs)) nrow(obs) else if (!is.null(X)) nrow(X) else 0L
  n_var <- if (!is.null(var)) nrow(var) else if (!is.null(X)) ncol(X) else 0L
  obs_names <- if (!is.null(X) && !is.null(rownames(X))) rownames(X) else NULL
  var_names <- if (!is.null(X) && !is.null(colnames(X))) colnames(X) else NULL
  obs <- .as_annot_df(obs, n_obs, obs_names)
  var <- .as_annot_df(var, n_var, var_names)
  x <- structure(
    list(X = X, obs = obs, var = var, obsm = as.list(obsm), varm = as.list(varm),
         obsp = as.list(obsp), varp = as.list(varp), layers = as.list(layers),
         uns = as.list(uns), raw = raw),
    class = "AnnR"
  )
  .validate_annr(x)
  x
}

.as_annot_df <- function(df, n, names) {
  if (is.null(df)) {
    df <- data.frame(row.names = names %||% .default_names(n))
    if (nrow(df) == 0L && n > 0L) df <- data.frame(row.names = .default_names(n))
    return(df)
  }
  if (!is.data.frame(df)) df <- as.data.frame(df, optional = TRUE, stringsAsFactors = FALSE)
  if (.row_names_are_automatic(df)) {
    rownames(df) <- names %||% .default_names(nrow(df))
  }
  df
}

.row_names_are_automatic <- function(df) {
  rn <- .row_names_info(df, type = 1L)
  is.integer(rn) && length(rn) <= 2L && (length(rn) == 0L || is.na(rn[1L]) || rn[1L] < 0L) &&
    is.null(attr(df, "h5ad_index"))
}

.validate_annr <- function(x) {
  n_obs <- nrow(x$obs)
  n_var <- nrow(x$var)
  chk <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  if (!is.null(x$X)) {
    chk(length(dim(x$X)) == 2L, "`X` must be two-dimensional")
    chk(nrow(x$X) == n_obs && ncol(x$X) == n_var,
        sprintf("`X` is %d x %d but obs/var describe %d x %d", nrow(x$X), ncol(x$X), n_obs, n_var))
  }
  for (k in names(x$layers)) {
    l <- x$layers[[k]]
    chk(NROW(l) == n_obs && NCOL(l) == n_var, sprintf("layer '%s' has the wrong shape", k))
  }
  for (k in names(x$obsm)) chk(NROW(x$obsm[[k]]) == n_obs, sprintf("obsm '%s' must have %d rows", k, n_obs))
  for (k in names(x$varm)) chk(NROW(x$varm[[k]]) == n_var, sprintf("varm '%s' must have %d rows", k, n_var))
  for (k in names(x$obsp)) {
    chk(length(dim(x$obsp[[k]])) == 2L && all(dim(x$obsp[[k]]) == n_obs), sprintf("obsp '%s' must be %d x %d", k, n_obs, n_obs))
  }
  for (k in names(x$varp)) {
    chk(length(dim(x$varp[[k]])) == 2L && all(dim(x$varp[[k]]) == n_var), sprintf("varp '%s' must be %d x %d", k, n_var, n_var))
  }
  if (!is.null(x$raw) && !is.null(x$raw$X)) {
    chk(length(dim(x$raw$X)) == 2L, "`raw$X` must be two-dimensional")
    chk(nrow(x$raw$X) == n_obs, "`raw$X` must have one row per observation")
    if (!is.null(x$raw$var)) {
      chk(nrow(x$raw$var) == ncol(x$raw$X), "`raw$var` must have one row per raw variable")
    }
  }
  invisible(TRUE)
}

#' @export
dim.AnnR <- function(x) c(nrow(x$obs), nrow(x$var))

#' @export
dimnames.AnnR <- function(x) list(rownames(x$obs), rownames(x$var))

#' Observation and variable names
#'
#' @param x An [AnnR] object.
#' @return A character vector.
#' @export
obs_names <- function(x) rownames(x$obs)

#' @rdname obs_names
#' @export
var_names <- function(x) rownames(x$var)

.describe_matrix <- function(m) {
  if (is.null(m)) return("NULL")
  cls <- class(m)[1L]
  if (.is_sparse(m)) {
    sprintf("%s %d x %d, %s non-zeros", cls, nrow(m), ncol(m),
            format(length(m@x), big.mark = ","))
  } else {
    sprintf("%s %s", cls, paste(dim(m), collapse = " x "))
  }
}

#' @export
print.AnnR <- function(x, ...) {
  d <- dim(x)
  cat(sprintf("AnnR object with n_obs x n_vars = %d x %d\n", d[1L], d[2L]))
  show_keys <- function(label, keys) {
    if (length(keys)) {
      txt <- paste0("'", keys, "'", collapse = ", ")
      if (nchar(txt) > 200L) txt <- paste0(substr(txt, 1L, 197L), "...")
      cat(sprintf("    %s: %s\n", label, txt))
    }
  }
  show_keys("obs", names(x$obs))
  show_keys("var", names(x$var))
  show_keys("uns", names(x$uns))
  show_keys("obsm", names(x$obsm))
  show_keys("varm", names(x$varm))
  show_keys("layers", names(x$layers))
  show_keys("obsp", names(x$obsp))
  show_keys("varp", names(x$varp))
  if (!is.null(x$raw)) {
    cat(sprintf("    raw: %s vars\n", if (is.null(x$raw$var)) NCOL(x$raw$X) else nrow(x$raw$var)))
  }
  cat(sprintf("    X: %s", .describe_matrix(x$X)))
  types <- attr(x, "h5ad_types")
  if (!is.null(types)) {
    dt <- types$dtype[types$path %in% c("X", "X/data")]
    dt <- dt[nzchar(dt)]
    if (length(dt)) cat(sprintf(" (stored as %s)", dt[1L]))
  }
  cat("\n")
  invisible(x)
}

.subset_rows <- function(m, i) {
  if (is.null(m)) return(NULL)
  if (is.data.frame(m)) {
    out <- m[i, , drop = FALSE]
    hi <- attr(m, "h5ad_index")
    if (!is.null(hi)) attr(out, "h5ad_index") <- hi[i]
    return(out)
  }
  if (is.null(dim(m))) return(m[i])
  m[i, , drop = FALSE]
}

.resolve_index <- function(idx, n, names, what) {
  if (missing(idx) || is.null(idx)) return(seq_len(n))
  if (is.character(idx)) {
    pos <- match(idx, names)
    if (anyNA(pos)) stop("unknown ", what, " names: ", paste(utils::head(idx[is.na(pos)], 5L), collapse = ", "), call. = FALSE)
    return(pos)
  }
  if (is.logical(idx)) {
    if (length(idx) != n) stop("logical ", what, " index must have length ", n, call. = FALSE)
    return(which(idx))
  }
  idx <- as.integer(idx)
  if (any(idx < 0L, na.rm = TRUE)) idx <- seq_len(n)[idx]
  if (anyNA(idx) || any(idx < 1L | idx > n)) stop(what, " index out of bounds", call. = FALSE)
  idx
}

#' Subset an AnnR object
#'
#' `x[i, j]` subsets observations `i` and variables `j` consistently across
#' `X`, layers, `obs`, `var`, `obsm`, `varm`, `obsp`, `varp` and `raw`
#' (`raw` keeps all of its variables). Indices may be integer, logical or
#' names.
#'
#' @param x An [AnnR] object.
#' @param i,j Observation and variable indices.
#' @param ... Unused.
#' @param drop Unused.
#' @return An [AnnR] object.
#' @export
`[.AnnR` <- function(x, i, j, ..., drop = FALSE) {
  d <- dim(x)
  ii <- .resolve_index(if (missing(i)) NULL else i, d[1L], obs_names(x), "obs")
  jj <- .resolve_index(if (missing(j)) NULL else j, d[2L], var_names(x), "var")
  out <- unclass(x)
  if (!is.null(out$X)) out$X <- out$X[ii, jj, drop = FALSE]
  out$layers <- lapply(out$layers, function(l) l[ii, jj, drop = FALSE])
  out$obs <- .subset_rows(out$obs, ii)
  out$var <- .subset_rows(out$var, jj)
  out$obsm <- lapply(out$obsm, .subset_rows, i = ii)
  out$varm <- lapply(out$varm, .subset_rows, i = jj)
  out$obsp <- lapply(out$obsp, function(m) m[ii, ii, drop = FALSE])
  out$varp <- lapply(out$varp, function(m) m[jj, jj, drop = FALSE])
  if (!is.null(out$raw) && !is.null(out$raw$X)) out$raw$X <- out$raw$X[ii, , drop = FALSE]
  structure(out, class = "AnnR", h5ad_types = attr(x, "h5ad_types"))
}

#' @export
`$<-.AnnR` <- function(x, name, value) {
  y <- unclass(x)
  y[[name]] <- value
  structure(y, class = "AnnR", h5ad_types = attr(x, "h5ad_types"))
}

#' @export
`[[<-.AnnR` <- function(x, i, value) {
  y <- unclass(x)
  y[[i]] <- value
  structure(y, class = "AnnR", h5ad_types = attr(x, "h5ad_types"))
}
