#' Number of threads used by AnnR
#'
#' Returns `getOption("AnnR.threads")` when set, otherwise the number of
#' available cores capped at 16.
#'
#' @return An integer.
#' @export
annr_threads <- function() {
  opt <- getOption("AnnR.threads")
  if (!is.null(opt)) return(max(1L, as.integer(opt)))
  n <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  if (is.na(n)) n <- 1L
  max(1L, min(16L, as.integer(n)))
}

#' Build information
#'
#' Reports the HDF5 library version AnnR is linked against and which codecs
#' are compiled in.
#'
#' @return A named character vector.
#' @export
annr_build_info <- function() {
  cpp_build_info()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.is_sparse <- function(x) methods::is(x, "sparseMatrix")

.check_slot_sel <- function(x, nm) {
  ok <- (is.logical(x) && length(x) == 1L && !is.na(x)) || is.character(x)
  if (!ok) stop("`", nm, "` must be TRUE, FALSE or a character vector of keys", call. = FALSE)
  x
}

.default_names <- function(n) as.character(seq_len(n) - 1L)

.emit_warnings <- function(msgs, file) {
  msgs <- unique(msgs)
  if (!length(msgs)) return(invisible())
  shown <- utils::head(msgs, 10L)
  more <- if (length(msgs) > 10L) sprintf("\n  ... and %d more", length(msgs) - 10L) else ""
  warning(sprintf("while reading '%s':\n  %s%s", basename(file),
                  paste(shown, collapse = "\n  "), more), call. = FALSE)
}
