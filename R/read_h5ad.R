#' Read an .h5ad file
#'
#' Reads an AnnData `.h5ad` file written by any version of Python `anndata`
#' (including pre-0.7 legacy layouts) with a parallel native reader.
#'
#' Orientation follows Python: in the returned [AnnR] object `X` is
#' `n_obs x n_vars`. CSR matrices become `dgRMatrix` and CSC matrices
#' `dgCMatrix` without any copy. When `as = "SingleCellExperiment"` or
#' `"Seurat"`, matrices are returned features x cells; sparse matrices are
#' reinterpreted (again without copying) and dense matrices are read directly
#' in that layout.
#'
#' Data types: floating point data become doubles, integers up to 32 bits
#' become integers, 64-bit integers become integers when every value fits and
#' doubles otherwise, booleans become logicals. The original on-disk dtypes are
#' recorded in `attr(x, "h5ad_types")` and reused by [write_h5ad()].
#'
#' @param file Path to an `.h5ad` file.
#' @param as Output class: `"AnnR"` (default), `"SingleCellExperiment"` or
#'   `"Seurat"`.
#' @param X Read `X`? (`TRUE`/`FALSE`)
#' @param layers,obsm,varm,obsp,varp,uns `TRUE` to read all entries, `FALSE`
#'   to skip the slot, or a character vector naming the entries to read.
#' @param raw Read `raw`? (`TRUE`/`FALSE`)
#' @param threads Number of threads for decoding.
#' @param ... Passed to the converter selected by `as`
#'   ([as_SingleCellExperiment()] or [as_Seurat()]).
#' @return An [AnnR] object, a `SingleCellExperiment` or a `Seurat` object.
#' @examples
#' \dontrun{
#' ad <- read_h5ad("pbmc.h5ad")
#' sce <- read_h5ad("pbmc.h5ad", as = "SingleCellExperiment")
#' small <- read_h5ad("pbmc.h5ad", layers = "counts", raw = FALSE, uns = FALSE)
#' }
#' @export
read_h5ad <- function(file, as = c("AnnR", "SingleCellExperiment", "Seurat"),
                      X = TRUE, layers = TRUE, obsm = TRUE, varm = TRUE, obsp = TRUE,
                      varp = TRUE, uns = TRUE, raw = TRUE, threads = annr_threads(), ...) {
  as <- match.arg(as)
  if (!is.character(file) || length(file) != 1L) stop("`file` must be a single path", call. = FALSE)
  path <- normalizePath(path.expand(file), mustWork = FALSE)
  if (!file.exists(path)) stop("file does not exist: ", file, call. = FALSE)
  slots <- list(
    X = .check_slot_sel(X, "X"), layers = .check_slot_sel(layers, "layers"),
    obsm = .check_slot_sel(obsm, "obsm"), varm = .check_slot_sel(varm, "varm"),
    obsp = .check_slot_sel(obsp, "obsp"), varp = .check_slot_sel(varp, "varp"),
    uns = .check_slot_sel(uns, "uns"), raw = .check_slot_sel(raw, "raw")
  )
  .reclaim_buffers()
  feature_major <- as != "AnnR"
  opts <- list(
    threads = as.integer(threads),
    strict = isTRUE(getOption("AnnR.strict", FALSE)),
    direct = !isFALSE(getOption("AnnR.direct", TRUE)),
    x_colmajor = !feature_major,
    slots = slots
  )
  res <- cpp_read_h5ad(path, opts)
  .emit_warnings(attr(res, "annr_warnings"), path)
  types <- as.data.frame(attr(res, "annr_types"), stringsAsFactors = FALSE)
  ad <- .finalize_anndata(res, types, feature_major = feature_major)
  if (isTRUE(getOption("AnnR.debug", FALSE))) {
    attr(ad, "annr_stats") <- attr(res, "annr_stats")
    attr(ad, "annr_fallbacks") <- attr(res, "annr_fallbacks")
  }
  switch(as,
    AnnR = ad,
    SingleCellExperiment = as_SingleCellExperiment(ad, ...),
    Seurat = as_Seurat(ad, ...)
  )
}

# Convert the raw list returned by the C++ reader into an AnnR object.
.finalize_anndata <- function(res, types, feature_major = FALSE) {
  r <- lapply(stats::setNames(names(res), names(res)), function(k) res[[k]])
  recarray_paths <- types$path[types$version == "legacy-recarray"]
  legacy_df <- FALSE
  for (k in c("obs", "var")) {
    if (k %in% recarray_paths && is.data.frame(r[[k]])) {
      r[[k]] <- .legacy_df_from_recarray(r[[k]])
      legacy_df <- TRUE
    }
  }
  if (!is.null(r$raw)) {
    attr(r$raw, "annr_anndata") <- NULL
    if ("raw/var" %in% recarray_paths && is.data.frame(r$raw$var)) {
      r$raw$var <- .legacy_df_from_recarray(r$raw$var)
    }
  }
  for (k in c("obsm", "varm")) {
    if (is.data.frame(r[[k]])) r[[k]] <- .split_compound_mapping(r[[k]])
  }
  r <- .finalize_tree(r)
  n_obs <- .infer_n(r, "obs", feature_major)
  n_var <- .infer_n(r, "var", feature_major)
  obs <- r$obs %||% data.frame(row.names = .default_names(n_obs))
  var <- r$var %||% data.frame(row.names = .default_names(n_var))
  uns <- r$uns %||% list()
  if (legacy_df) {
    cleaned <- .clean_uns_categories(obs, var, uns)
    obs <- cleaned$obs
    var <- cleaned$var
    uns <- cleaned$uns
  }
  obsp <- r$obsp %||% list()
  moved <- .move_adjacency(uns, obsp, types)
  uns <- moved$uns
  obsp <- moved$obsp
  types <- moved$types
  raw <- NULL
  if (!is.null(r$raw)) {
    raw <- list(X = r$raw$X, var = r$raw$var, varm = r$raw$varm %||% list())
    if (is.data.frame(raw$varm)) raw$varm <- .split_compound_mapping(raw$varm)
    if (is.null(raw$var) && !is.null(raw$X)) {
      raw$var <- data.frame(row.names = .default_names(if (feature_major) nrow(raw$X) else ncol(raw$X)))
    }
  }
  structure(
    list(X = r$X, obs = obs, var = var, obsm = r$obsm %||% list(), varm = r$varm %||% list(),
         obsp = obsp, varp = r$varp %||% list(), layers = r$layers %||% list(), uns = uns,
         raw = raw),
    class = "AnnR",
    h5ad_types = types,
    annr_feature_major = if (feature_major) TRUE else NULL
  )
}

.infer_n <- function(r, axis, feature_major) {
  df <- r[[axis]]
  if (is.data.frame(df)) return(nrow(df))
  x <- r$X
  if (is.null(x)) return(0L)
  d <- dim(x)
  if (feature_major && !.is_sparse(x)) d <- rev(d)
  if (axis == "obs") d[1L] else d[2L]
}

# Walk lists: repair data frames whose index was not unique and turn nested
# AnnData groups into AnnR objects.
.finalize_tree <- function(x) {
  if (is.data.frame(x)) return(.fix_df_index(x))
  if (is.list(x) && !is.object(x)) {
    nested <- isTRUE(attr(x, "annr_anndata"))
    at <- attributes(x)
    x <- lapply(x, .finalize_tree)
    attributes(x) <- at[intersect(names(at), c("names", "form", "length", "class"))]
    if (nested) {
      x$obs <- x$obs %||% data.frame(row.names = .default_names(if (is.null(x$X)) 0L else nrow(x$X)))
      x$var <- x$var %||% data.frame(row.names = .default_names(if (is.null(x$X)) 0L else ncol(x$X)))
      for (k in c("obsm", "varm", "obsp", "varp", "layers", "uns")) x[[k]] <- x[[k]] %||% list()
      x <- structure(x, class = "AnnR")
    }
    return(x)
  }
  if (inherits(x, "h5ad_awkward")) return(x)
  x
}

.fix_df_index <- function(df) {
  idx <- attr(df, "h5ad_index")
  for (k in seq_along(df)) {
    if (is.data.frame(df[[k]])) df[[k]] <- .fix_df_index(df[[k]])
  }
  if (is.null(idx)) return(df)
  idx2 <- idx
  idx2[is.na(idx2)] <- "NA"
  attr(df, "row.names") <- make.unique(idx2)
  df
}

# Pre-0.7 obs/var: first field of the record array is the index.
.legacy_df_from_recarray <- function(df) {
  idx <- df[[1L]]
  out <- df[-1L]
  if (anyDuplicated(idx) || anyNA(idx)) {
    attr(out, "h5ad_index") <- idx
    idx[is.na(idx)] <- "NA"
    attr(out, "row.names") <- make.unique(idx)
  } else {
    attr(out, "row.names") <- idx
  }
  out
}

# Pre-0.7 obsm/varm were record arrays; each field is one entry.
.split_compound_mapping <- function(df) {
  out <- lapply(df, function(col) col)
  names(out) <- names(df)
  out
}

# anndata.compat._clean_uns: categories stored in uns as "<col>_categories".
.clean_uns_categories <- function(obs, var, uns) {
  drop <- character()
  for (key in names(uns)) {
    if (!endsWith(key, "_categories")) next
    name <- sub("_categories$", "", key)
    cats <- as.character(uns[[key]])
    for (axis in c("obs", "var")) {
      df <- if (axis == "obs") obs else var
      if (!name %in% names(df)) next
      codes <- df[[name]]
      if (!is.numeric(codes) || any(codes >= length(cats), na.rm = TRUE)) next
      f <- codes + 1L
      f[codes < 0] <- NA
      df[[name]] <- structure(as.integer(f), levels = cats, class = "factor")
      if (axis == "obs") obs <- df else var <- df
      drop <- c(drop, key)
    }
  }
  uns[unique(drop)] <- NULL
  list(obs = obs, var = var, uns = uns)
}

# anndata moves uns["neighbors"]["distances"/"connectivities"] to obsp.
.move_adjacency <- function(uns, obsp, types) {
  nb <- uns$neighbors
  if (!is.list(nb) || is.data.frame(nb)) return(list(uns = uns, obsp = obsp, types = types))
  for (k in c("distances", "connectivities")) {
    m <- nb[[k]]
    if (!is.null(m) && length(dim(m)) == 2L) {
      warning("Moving element from uns$neighbors$", k, " to obsp$", k, call. = FALSE)
      obsp[[k]] <- m
      uns$neighbors[[k]] <- NULL
      # carry the dtype records along with the element
      from <- paste0("uns/neighbors/", k)
      hit <- types$path == from | startsWith(types$path, paste0(from, "/"))
      if (any(hit)) {
        add <- types[hit, , drop = FALSE]
        add$path <- sub(from, paste0("obsp/", k), add$path, fixed = TRUE)
        types <- rbind(types, add)
      }
    }
  }
  list(uns = uns, obsp = obsp, types = types)
}

#' Inspect the structure of an .h5ad file
#'
#' Lists every group and dataset with its anndata encoding, dtype, shape,
#' chunking, filters and stored size, without reading any data.
#'
#' @param file Path to an `.h5ad` (or any HDF5) file.
#' @return A data frame with one row per HDF5 object.
#' @export
h5ad_info <- function(file) {
  path <- normalizePath(path.expand(file), mustWork = TRUE)
  as.data.frame(cpp_h5_ls(path), stringsAsFactors = FALSE)
}

#' Read a single element of an .h5ad file
#'
#' Reads one element, such as `"obsm/X_umap"`, `"layers/counts"`, `"obs"` or
#' `"uns/neighbors"`, without loading the rest of the file. Matrices keep the
#' on-disk (Python) orientation.
#'
#' @param file Path to an `.h5ad` file.
#' @param element Path of the element inside the file.
#' @param threads Number of threads for decoding.
#' @return The element as an R object (see [read_h5ad()] for the mapping).
#' @examples
#' \dontrun{
#' umap <- read_h5ad_elem("pbmc.h5ad", "obsm/X_umap")
#' }
#' @export
read_h5ad_elem <- function(file, element, threads = annr_threads()) {
  .reclaim_buffers()
  path <- normalizePath(path.expand(file), mustWork = FALSE)
  if (!file.exists(path)) stop("file does not exist: ", file, call. = FALSE)
  opts <- list(threads = as.integer(threads), strict = isTRUE(getOption("AnnR.strict", FALSE)),
               direct = !isFALSE(getOption("AnnR.direct", TRUE)))
  res <- cpp_read_elem(path, sub("^/", "", element), opts)
  if (is.null(res)) return(NULL)
  .emit_warnings(attr(res, "annr_warnings"), path)
  for (a in c("annr_types", "annr_warnings", "annr_fallbacks", "annr_stats")) attr(res, a) <- NULL
  .finalize_tree(res)
}

# Large arrays live outside R's accounted heap, so R may not notice that a
# previous (now unreferenced) object holds gigabytes. Before a new read, run a
# full collection when more than `limit` bytes of such buffers are alive.
.reclaim_buffers <- function(limit = getOption("AnnR.gc_threshold", 2e9)) {
  if (cpp_external_bytes() > limit) gc(verbose = FALSE)
  invisible()
}
