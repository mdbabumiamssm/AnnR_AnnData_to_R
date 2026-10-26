# Dump an AnnR object read from an .h5ad file into NPY files + a JSON manifest
# so that an independent Python script can compare it with anndata.
#
# Usage: Rscript r_dump.R <file.h5ad> <out_dir> [threads]
suppressPackageStartupMessages({
  library(AnnR)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
infile <- args[[1]]
outdir <- args[[2]]
threads <- if (length(args) >= 3) as.integer(args[[3]]) else annr_threads()
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counter <- 0L
new_file <- function(ext) {
  counter <<- counter + 1L
  sprintf("f%06d.%s", counter, ext)
}

write_npy <- function(x, dtype) {
  fn <- new_file("npy")
  con <- file(file.path(outdir, fn), "wb")
  on.exit(close(con))
  header <- sprintf("{'descr': '%s', 'fortran_order': False, 'shape': (%.0f,), }", dtype, length(x))
  total <- 10 + nchar(header) + 1
  pad <- (64 - total %% 64) %% 64
  header <- paste0(header, strrep(" ", pad), "\n")
  writeBin(as.raw(0x93), con)
  writeBin(charToRaw("NUMPY"), con)
  writeBin(as.raw(c(1L, 0L)), con)
  writeBin(nchar(header), con, size = 2, endian = "little")
  writeBin(charToRaw(header), con)
  chunk <- 50000000
  n <- length(x)
  if (n > 0) {
    for (s in seq(1, n, by = chunk)) {
      e <- min(n, s + chunk - 1)
      piece <- x[s:e]
      switch(dtype,
        "<f8" = writeBin(as.double(piece), con, size = 8, endian = "little"),
        "<i4" = writeBin(as.integer(piece), con, size = 4, endian = "little"),
        "<c16" = writeBin(as.complex(piece), con, endian = "little")
      )
    }
  }
  fn
}

write_strings <- function(x) {
  fn <- new_file("json")
  x <- as.vector(x)  # column-major flattening, like the numeric dumps
  writeLines(enc2utf8(as.character(jsonlite::toJSON(unname(enc2utf8(x)), na = "null"))),
             file.path(outdir, fn), useBytes = TRUE)
  fn
}

dump_atomic <- function(x) {
  node <- list(kind = "array", rtype = typeof(x), dim = if (is.null(dim(x))) NULL else as.list(dim(x)),
               length = length(x))
  node$data <- switch(typeof(x),
    double = write_npy(x, "<f8"),
    integer = write_npy(x, "<i4"),
    logical = write_npy(as.integer(x), "<i4"),
    character = write_strings(x),
    complex = write_npy(x, "<c16"),
    raw = write_npy(as.integer(x), "<i4"),
    stop("unsupported atomic type ", typeof(x))
  )
  node
}

dump_value <- function(x) {
  if (is.null(x)) return(list(kind = "null"))
  if (inherits(x, "AnnR")) {
    slots <- lapply(c(X = "X", obs = "obs", var = "var", obsm = "obsm", varm = "varm",
                      obsp = "obsp", varp = "varp", layers = "layers", uns = "uns", raw = "raw"),
                    function(k) dump_value(unclass(x)[[k]]))
    return(list(kind = "anndata", slots = slots))
  }
  if (is.data.frame(x)) {
    idx <- attr(x, "h5ad_index")
    rn <- if (is.null(idx)) attr(x, "row.names") else idx
    cols <- lapply(seq_along(x), function(i) dump_value(x[[i]]))
    names(cols) <- NULL
    return(list(kind = "dataframe", rownames = write_strings(as.character(rn)),
                rownames_type = typeof(attr(x, "row.names")),
                colnames = write_strings(names(x)), columns = cols, nrow = nrow(x)))
  }
  if (is.factor(x)) {
    return(list(kind = "factor", codes = write_npy(as.integer(unclass(x)), "<i4"),
                levels = write_strings(levels(x)), ordered = is.ordered(x)))
  }
  if (methods::is(x, "sparseMatrix")) {
    cls <- class(x)[1]
    node <- list(kind = "sparse", class = cls, dim = as.list(dim(x)))
    node$p <- write_npy(x@p, "<i4")
    node$idx <- write_npy(if (cls %in% c("dgRMatrix", "lgRMatrix")) x@j else x@i, "<i4")
    node$x <- if (is.logical(x@x)) write_npy(as.integer(x@x), "<i4") else write_npy(x@x, "<f8")
    node$xtype <- typeof(x@x)
    return(node)
  }
  if (inherits(x, "h5ad_awkward")) {
    return(list(kind = "awkward", form = attr(x, "form"), length = attr(x, "length"),
                buffers = lapply(unclass(x), dump_value)))
  }
  if (is.atomic(x)) return(dump_atomic(x))
  if (is.list(x)) {
    items <- lapply(x, dump_value)
    return(list(kind = "dict", keys = write_strings(names(x) %||% character()), items = unname(items)))
  }
  stop("cannot dump object of class ", paste(class(x), collapse = "/"))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

options(AnnR.debug = TRUE)
t0 <- proc.time()[["elapsed"]]
warns <- character()
ad <- withCallingHandlers(read_h5ad(infile, threads = threads),
                          warning = function(w) {
                            warns <<- c(warns, conditionMessage(w))
                            invokeRestart("muffleWarning")
                          })
elapsed <- proc.time()[["elapsed"]] - t0
stats <- attr(ad, "annr_stats")
manifest <- list(file = infile, read_seconds = elapsed, warnings = warns,
                 direct_reads = stats[[1]], fallback_reads = stats[[2]],
                 fallback_notes = unique(attr(ad, "annr_fallbacks")), root = dump_value(ad))
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null", digits = NA, pretty = FALSE),
           file.path(outdir, "manifest.json"), useBytes = TRUE)
cat(sprintf("dumped %s in %.2fs (%d files)\n", basename(infile), elapsed, counter))
