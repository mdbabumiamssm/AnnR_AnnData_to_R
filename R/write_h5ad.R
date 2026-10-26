#' Write an .h5ad file
#'
#' Writes an [AnnR], `SingleCellExperiment` or `Seurat` object to an AnnData
#' `.h5ad` file that Python `anndata` reads natively.
#'
#' Data types: when the object was read with [read_h5ad()], the original
#' on-disk dtype of every element (e.g. `float32` counts, `int8` codes, `uint16`
#' indices) is reproduced, but only when the conversion is lossless; otherwise
#' the natural R type is written (`float64` for doubles, `int32` for integers,
#' `bool` for logicals). Vectors with `NA` are written with anndata's nullable
#' encodings.
#'
#' Scalars: a length-one vector in `uns` is written as a scalar unless it was
#' read as an array or is wrapped in [I()].
#'
#' The file is written to a temporary file in the same directory and renamed
#' into place only after it has been completed successfully. On Windows, an
#' existing destination is backed up during replacement; a failed replacement
#' attempts to restore that backup.
#'
#' @param x An [AnnR], `SingleCellExperiment` or `Seurat` object.
#' @param file Output path.
#' @param compression `"none"` (default, like anndata), `"gzip"`, `"lzf"` or
#'   `"zstd"` (zstd needs the hdf5plugin package to be read from Python).
#' @param compression_level Compression level (gzip 0-9, default 4; zstd 1-22,
#'   default 3).
#' @param shuffle Apply HDF5's byte-shuffle filter before compression.
#' @param threads Number of threads used to convert and compress.
#' @param overwrite Overwrite an existing file?
#' @param ... Passed to [as_AnnR()] when `x` is not an [AnnR] object.
#' @return The output path, invisibly.
#' @examples
#' \dontrun{
#' ad <- read_h5ad("in.h5ad")
#' write_h5ad(ad, "out.h5ad", compression = "gzip")
#' }
#' @export
write_h5ad <- function(x, file, compression = c("none", "gzip", "lzf", "zstd"),
                       compression_level = NULL, shuffle = FALSE, threads = annr_threads(),
                       overwrite = TRUE, ...) {
  compression <- match.arg(compression)
  ad <- if (inherits(x, "AnnR")) x else as_AnnR(x, ...)
  .validate_annr(ad)
  file <- path.expand(file)
  if (file.exists(file) && !isTRUE(overwrite)) stop("file exists: ", file, call. = FALSE)
  dir <- dirname(file)
  if (!dir.exists(dir)) stop("directory does not exist: ", dir, call. = FALSE)
  level <- compression_level %||% (if (compression == "zstd") 3L else 4L)
  comp <- list(compression = compression, level = as.integer(level), shuffle = isTRUE(shuffle))
  st <- .writer_state(attr(ad, "h5ad_types"), comp)

  tmp <- tempfile(pattern = ".annr_", tmpdir = dir, fileext = ".h5ad.tmp")
  w <- cpp_w_open(tmp, as.integer(threads))
  ok <- FALSE
  on.exit({
    if (!ok) {
      try(cpp_w_close(w), silent = TRUE)
      unlink(tmp)
    }
  }, add = TRUE)
  cpp_w_group(w, "/", "anndata", "0.1.0")
  .w_anndata(w, "", ad, st)
  cpp_w_close(w)
  # POSIX rename replaces the destination atomically. Never unlink it first.
  # On platforms that cannot replace an existing file, preserve a backup and
  # restore it if the final rename fails. Do not fall back to a partial copy.
  if (.Platform$OS.type == "windows" && file.exists(file)) {
    backup <- tempfile(pattern = ".annr_backup_", tmpdir = dir)
    if (!file.rename(file, backup)) stop("could not preserve the existing file", call. = FALSE)
    if (!file.rename(tmp, file)) {
      restored <- file.rename(backup, file)
      stop(if (restored) "could not replace the file; original restored" else
        paste0("could not replace the file; original preserved at ", backup), call. = FALSE)
    }
    unlink(backup)
  } else if (!file.rename(tmp, file)) {
    stop("could not move the written file into place; existing file preserved", call. = FALSE)
  }
  ok <- TRUE
  if (length(st$warnings)) {
    msgs <- unique(st$warnings)
    warning(paste(c("while writing:", utils::head(msgs, 10L)), collapse = "\n  "), call. = FALSE)
  }
  invisible(file)
}

# Writer state: dtype/encoding hints from the source file and warnings.
.writer_state <- function(types, comp) {
  hints <- new.env(hash = TRUE, parent = emptyenv())
  if (!is.null(types) && nrow(types)) {
    for (i in seq_len(nrow(types))) {
      assign(types$path[i], list(enc = types$encoding[i], ver = types$version[i],
                                 dtype = types$dtype[i], extra = types$extra[i]), envir = hints)
    }
  }
  st <- new.env(parent = emptyenv())
  st$hints <- hints
  st$comp <- comp
  st$warnings <- character()
  st
}

.hint <- function(st, path) get0(path, envir = st$hints, inherits = FALSE)

# dtype hint for a value; fixed-length string widths ("S7") are only reproduced
# for elements that were already in the modern encoding (legacy files are
# re-written with variable-length UTF-8, as Python anndata does).
.value_dtype <- function(h, x) {
  if (is.null(h)) return("")
  if (is.character(x) && startsWith(h$dtype, "S") && !nzchar(h$ver)) return("")
  h$dtype
}

.hint_dtype <- function(st, path) {
  h <- .hint(st, path)
  if (is.null(h)) "" else h$dtype
}

.w_warn <- function(st, msg) st$warnings <- c(st$warnings, msg)

.join <- function(prefix, name) if (nzchar(prefix)) paste0(prefix, "/", name) else name

.check_key <- function(k, path) {
  if (is.na(k) || !nzchar(k)) stop("empty element name in '", path, "'", call. = FALSE)
  if (grepl("/", k, fixed = TRUE)) {
    stop("element names may not contain '/': '", k, "' in '", path, "'", call. = FALSE)
  }
  k
}

.w_anndata <- function(w, prefix, ad, st) {
  ad <- unclass(ad)
  # Data frames in obsm/varm are indexed by the obs/var names in anndata.
  for (k in names(ad$obsm)) {
    if (is.data.frame(ad$obsm[[k]]) && .row_names_are_automatic(ad$obsm[[k]])) {
      rownames(ad$obsm[[k]]) <- rownames(ad$obs)
    }
  }
  for (k in names(ad$varm)) {
    if (is.data.frame(ad$varm[[k]]) && .row_names_are_automatic(ad$varm[[k]])) {
      rownames(ad$varm[[k]]) <- rownames(ad$var)
    }
  }
  if (!is.null(ad$X)) .w_elem(w, .join(prefix, "X"), ad$X, st)
  .w_dataframe(w, .join(prefix, "obs"), ad$obs, st)
  .w_dataframe(w, .join(prefix, "var"), ad$var, st)
  for (k in c("obsm", "varm", "obsp", "varp", "layers", "uns")) {
    .w_mapping(w, .join(prefix, k), ad[[k]] %||% list(), st)
  }
  if (!is.null(ad$raw)) {
    raw <- ad$raw
    p <- .join(prefix, "raw")
    cpp_w_group(w, p, "raw", "0.1.0")
    if (!is.null(raw$X)) .w_elem(w, paste0(p, "/X"), raw$X, st)
    rv <- raw$var %||% data.frame(row.names = .default_names(NCOL(raw$X)))
    .w_dataframe(w, paste0(p, "/var"), rv, st)
    .w_mapping(w, paste0(p, "/varm"), raw$varm %||% list(), st)
  }
}

.w_mapping <- function(w, path, m, st) {
  cpp_w_group(w, path, "dict", "0.1.0")
  if (!length(m)) return(invisible())
  nms <- names(m)
  if (is.null(nms)) nms <- as.character(seq_along(m) - 1L)
  for (i in seq_along(m)) {
    k <- .check_key(nms[[i]], path)
    .w_elem(w, paste0(path, "/", k), m[[i]], st)
  }
}

.shape_of <- function(x) {
  d <- dim(x)
  as.double(if (is.null(d)) length(x) else d)
}

.w_elem <- function(w, path, x, st, in_df = FALSE) {
  h <- .hint(st, path)
  if (is.null(x)) {
    cpp_w_null(w, path)
  } else if (inherits(x, "AnnR")) {
    cpp_w_group(w, path, "anndata", "0.1.0")
    .w_anndata(w, path, x, st)
  } else if (is.data.frame(x)) {
    if (!is.null(h) && identical(h$enc, "rec-array")) .w_recarray(w, path, x, st) else .w_dataframe(w, path, x, st)
  } else if (is.factor(x)) {
    .w_categorical(w, path, x, st)
  } else if (methods::is(x, "sparseMatrix")) {
    if (!is.null(h) && identical(h$enc, "array")) {
      # stored dense on disk (e.g. a dense obsp that went through colPairs)
      .w_elem(w, path, as.matrix(x), st, in_df)
    } else {
      .w_sparse(w, path, x, st)
    }
  } else if (methods::is(x, "Matrix")) {
    .w_elem(w, path, as.matrix(x), st, in_df)
  } else if (methods::is(x, "DelayedArray")) {
    .w_elem(w, path, .realize_delayed(x), st, in_df)
  } else if (inherits(x, "h5ad_awkward")) {
    .w_awkward(w, path, x, st)
  } else if (is.atomic(x)) {
    .w_atomic(w, path, x, st, h, in_df)
  } else if (is.list(x)) {
    .w_mapping(w, path, x, st)
  } else {
    .w_warn(st, sprintf("'%s': objects of class %s cannot be stored and were skipped", path,
                        paste(class(x), collapse = "/")))
  }
}

.realize_delayed <- function(x) {
  sparse <- tryCatch(isTRUE(DelayedArray::is_sparse(x)), error = function(e) FALSE)
  if (sparse) methods::as(x, "dgCMatrix") else as.matrix(x)
}

.as_plain_atomic <- function(x, path, st) {
  if (inherits(x, c("Date", "POSIXt", "difftime"))) {
    .w_warn(st, sprintf("'%s': %s values were written as strings", path, class(x)[1L]))
    return(structure(as.character(x), dim = dim(x)))
  }
  if (inherits(x, "integer64")) {
    .w_warn(st, sprintf("'%s': integer64 values were written as float64", path))
    return(as.double(x))
  }
  if (is.raw(x)) return(as.integer(x))
  d <- dim(x)
  bits <- attr(x, "h5ad_int64")
  attributes(x) <- NULL
  if (!is.null(d)) dim(x) <- d
  if (!is.null(bits)) attr(x, "h5ad_int64") <- bits
  x
}

.w_atomic <- function(w, path, x, st, h, in_df) {
  as_is <- inherits(x, "AsIs")
  x <- .as_plain_atomic(x, path, st)
  enc_hint <- if (is.null(h)) "" else h$enc
  scalar_hint <- enc_hint %in% c("numeric-scalar", "string")
  array_hint <- enc_hint %in% c("array", "string-array", "nullable-integer", "nullable-boolean",
                                "nullable-string-array") || (identical(enc_hint, "") && !is.null(h))
  is_scalar <- !in_df && is.null(dim(x)) && length(x) == 1L && !as_is &&
    (scalar_hint || (!array_hint))
  if (is_scalar && !is.na(x)) {
    enc <- if (is.character(x)) "string" else "numeric-scalar"
    cpp_w_array(w, path, x, numeric(0), TRUE, .value_dtype(h, x), enc, "0.2.0", NULL, FALSE)
    return(invisible())
  }
  has_na <- anyNA(x) && !is.double(x) && !is.complex(x)
  nullable_hint <- startsWith(enc_hint, "nullable-")
  if ((has_na || nullable_hint) && is.null(dim(x)) && (is.integer(x) || is.logical(x) || is.character(x))) {
    .w_nullable(w, path, x, st, h)
    return(invisible())
  }
  if (has_na) {
    # NA in a multi-dimensional integer/logical/character array: anndata has
    # no nullable N-d encoding.
    if (is.character(x)) {
      .w_warn(st, sprintf("'%s': NA strings written as empty strings", path))
    } else {
      .w_warn(st, sprintf("'%s': NA values written as NaN in a float64 array", path))
      x <- array(as.double(x), dim = dim(x))
    }
  }
  enc <- if (is.character(x)) "string-array" else "array"
  hint <- .value_dtype(h, x)
  cpp_w_array(w, path, x, .shape_of(x), TRUE, hint, enc, "0.2.0", st$comp, FALSE)
  invisible()
}

.w_nullable <- function(w, path, x, st, h) {
  mask <- is.na(x)
  if (is.character(x)) {
    enc <- "nullable-string-array"
    vals <- x
    vals[mask] <- ""
  } else if (is.logical(x)) {
    enc <- "nullable-boolean"
    vals <- x
    vals[mask] <- FALSE
  } else {
    enc <- "nullable-integer"
    vals <- x
    vals[mask] <- 0L
  }
  cpp_w_group(w, path, enc, "0.1.0")
  if (enc == "nullable-string-array") {
    na_value <- if (!is.null(h) && nzchar(h$extra)) h$extra else "NA"
    cpp_w_attr(w, path, "na-value", na_value, "string")
  }
  vhint <- .hint_dtype(st, paste0(path, "/values"))
  if (!nzchar(vhint) && enc == "nullable-integer") vhint <- "int64"
  cpp_w_array(w, paste0(path, "/values"), vals, as.double(length(vals)), TRUE, vhint,
              if (is.character(vals)) "string-array" else "array", "0.2.0", st$comp, FALSE)
  cpp_w_array(w, paste0(path, "/mask"), mask, as.double(length(mask)), TRUE, "bool", "array",
              "0.2.0", st$comp, FALSE)
}

# pandas picks the smallest signed code type that holds the category count.
.codes_dtype <- function(ncat) {
  if (ncat < 127L) "int8" else if (ncat < 32767L) "int16" else if (ncat < 2147483647) "int32" else "int64"
}

.w_categorical <- function(w, path, f, st) {
  cpp_w_group(w, path, "categorical", "0.2.0")
  cpp_w_attr(w, path, "ordered", is.ordered(f), "bool")
  lv <- levels(f)
  codes <- as.integer(unclass(f)) - 1L
  codes[is.na(codes)] <- -1L
  cpp_w_array(w, paste0(path, "/codes"), codes, as.double(length(codes)), TRUE,
              .codes_dtype(length(lv)), "array", "0.2.0", st$comp, FALSE)
  cats_path <- paste0(path, "/categories")
  cdtype <- .hint_dtype(st, cats_path)
  cats <- .restore_category_type(lv, cdtype)
  cpp_w_array(w, cats_path, cats, as.double(length(cats)), TRUE, cdtype,
              if (is.character(cats)) "string-array" else "array", "0.2.0", st$comp, FALSE)
}

# Categories that were numeric/bool on disk come back as R level strings.
.restore_category_type <- function(lv, dtype) {
  if (!nzchar(dtype) || dtype == "str" || startsWith(dtype, "S")) return(lv)
  if (dtype == "bool") {
    v <- as.logical(lv)
    if (!anyNA(v)) return(v)
    return(lv)
  }
  v <- suppressWarnings(as.numeric(lv))
  if (anyNA(v) || !identical(as.character(v), lv)) return(lv)
  if (grepl("^u?int", dtype) && all(v == trunc(v)) && all(abs(v) <= .Machine$integer.max)) {
    return(as.integer(v))
  }
  v
}

.w_dataframe <- function(w, path, df, st) {
  h <- .hint(st, path)
  idx_key <- if (!is.null(h) && nzchar(h$extra)) h$extra else "_index"
  cols <- names(df)
  if (anyDuplicated(cols)) stop("'", path, "': duplicated column names", call. = FALSE)
  if ("_index" %in% cols) stop("'", path, "': '_index' is a reserved column name", call. = FALSE)
  if (idx_key %in% cols) idx_key <- "_index"
  for (k in cols) .check_key(k, path)
  cpp_w_group(w, path, "dataframe", "0.2.0")
  cpp_w_attr(w, path, "_index", idx_key, "string")
  cpp_w_attr(w, path, "column-order", as.character(cols), "strings")
  index <- .df_index(df)
  ipath <- paste0(path, "/", idx_key)
  ih <- .hint(st, ipath)
  idtype <- if (is.null(ih)) "" else ih$dtype
  if (grepl("^u?int", idtype)) {
    iv <- suppressWarnings(as.numeric(index))
    if (!anyNA(iv) && identical(as.character(iv), index)) {
      cpp_w_array(w, ipath, iv, as.double(length(iv)), TRUE, idtype, "array", "0.2.0", st$comp, FALSE)
      index <- NULL
    }
  }
  if (!is.null(index)) {
    if (!is.null(ih) && identical(ih$enc, "nullable-string-array")) {
      .w_nullable(w, ipath, index, st, ih)
    } else {
      cpp_w_array(w, ipath, index, as.double(length(index)), TRUE, idtype, "string-array", "0.2.0",
                  st$comp, FALSE)
    }
  }
  for (k in cols) {
    col <- df[[k]]
    cp <- paste0(path, "/", k)
    if (!is.null(dim(col)) && !is.data.frame(col) && length(dim(col)) == 2L && ncol(col) != 1L) {
      .w_warn(st, sprintf("'%s': matrix columns are not supported in anndata data frames; skipped", cp))
      next
    }
    if (is.list(col) && !is.data.frame(col)) {
      flat <- tryCatch(unlist(col), error = function(e) NULL)
      if (is.atomic(flat) && length(flat) == length(col)) {
        col <- flat
      } else {
        .w_warn(st, sprintf("'%s': list columns cannot be stored; written as strings", cp))
        col <- vapply(col, function(v) paste(format(v), collapse = ","), character(1))
      }
    }
    if (!is.null(dim(col)) && !is.data.frame(col)) col <- as.vector(col)
    .w_elem(w, cp, col, st, in_df = TRUE)
  }
}

.df_index <- function(df) {
  idx <- attr(df, "h5ad_index")
  rn <- rownames(df)
  if (!is.null(idx) && length(idx) == nrow(df)) {
    idx2 <- idx
    idx2[is.na(idx2)] <- "NA"
    if (identical(make.unique(idx2), rn)) return(as.character(idx))
  }
  if (.row_names_are_automatic(df)) return(.default_names(nrow(df)))
  rn
}

.w_sparse <- function(w, path, m, st) {
  if (!methods::is(m, "generalMatrix")) m <- methods::as(m, "generalMatrix")
  if (methods::is(m, "TsparseMatrix")) m <- methods::as(m, "CsparseMatrix")
  if (methods::is(m, "nsparseMatrix")) m <- methods::as(m, "lMatrix")
  if (!(methods::is(m, "dsparseMatrix") || methods::is(m, "lsparseMatrix"))) m <- methods::as(m, "dMatrix")
  # Keep the CSR/CSC layout the element had on disk (e.g. after a detour
  # through SingleCellExperiment colPairs, which are column-compressed).
  h <- .hint(st, path)
  if (!is.null(h)) {
    if (identical(h$enc, "csr_matrix") && methods::is(m, "CsparseMatrix")) {
      m <- methods::as(m, "RsparseMatrix")
    } else if (identical(h$enc, "csc_matrix") && methods::is(m, "RsparseMatrix")) {
      m <- methods::as(m, "CsparseMatrix")
    }
  } else if (grepl("(^|/)(obsp|varp)/[^/]+$", path) && methods::is(m, "CsparseMatrix")) {
    # anndata/scanpy convention: pairwise graphs are CSR
    m <- methods::as(m, "RsparseMatrix")
  }
  if (methods::is(m, "RsparseMatrix")) {
    fmt <- "csr"
    idx <- m@j
  } else {
    fmt <- "csc"
    idx <- m@i
  }
  hints <- c(.hint_dtype(st, paste0(path, "/data")), .hint_dtype(st, paste0(path, "/indices")),
             .hint_dtype(st, paste0(path, "/indptr")))
  invisible(cpp_w_sparse(w, path, fmt, as.double(dim(m)), m@p, idx, m@x, hints, st$comp))
}

.w_recarray <- function(w, path, df, st) {
  cols <- lapply(df, function(c) if (is.factor(c)) as.character(c) else c)
  hints <- vapply(names(df), function(k) .hint_dtype(st, paste0(path, "/", k)), character(1))
  cpp_w_recarray(w, path, unname(cols), names(df), unname(hints))
}

.w_awkward <- function(w, path, x, st) {
  cpp_w_group(w, path, "awkward-array", "0.1.0")
  cpp_w_attr(w, path, "form", attr(x, "form"), "string")
  cpp_w_attr(w, path, "length", as.double(attr(x, "length")), "int64")
  buf <- unclass(x)
  for (k in names(buf)) {
    v <- buf[[k]]
    cpp_w_array(w, paste0(path, "/", k), v, .shape_of(v), TRUE, .hint_dtype(st, paste0(path, "/", k)),
                "array", "0.2.0", st$comp, FALSE)
  }
}
