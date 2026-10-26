# One benchmark measurement in a fresh R process.
# Usage: Rscript bench_r.R <tool> <op> <file> <out> <threads> <compression>
#   tool: annr | anndataR ; op: read | read_sce | write
a <- commandArgs(trailingOnly = TRUE)
tool <- a[[1]]; op <- a[[2]]; file <- a[[3]]; out <- a[[4]]
threads <- as.integer(a[[5]]); comp <- a[[6]]
status <- function(key) {
  s <- readLines("/proc/self/status")
  as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", grep(paste0("^", key), s, value = TRUE))) / 1024^2
}
suppressPackageStartupMessages({
  library(Matrix)
  if (tool == "annr") library(AnnR) else library(anndataR)
  library(SingleCellExperiment)
})
base <- status("VmRSS")
rd <- function(as) {
  if (tool == "annr") {
    AnnR::read_h5ad(file, as = if (as == "native") "AnnR" else "SingleCellExperiment", threads = threads)
  } else {
    anndataR::read_h5ad(file, as = if (as == "native") "InMemoryAnnData" else "SingleCellExperiment")
  }
}
res <- list(tool = tool, op = op, file = file, threads = threads, compression = comp)
if (op %in% c("read", "read_sce")) {
  t <- system.time(obj <- suppressWarnings(rd(if (op == "read") "native" else "sce")))[["elapsed"]]
} else {
  obj <- suppressWarnings(rd("native"))
  invisible(gc())
  if (file.exists(out)) unlink(out)
  t <- system.time({
    if (tool == "annr") AnnR::write_h5ad(obj, out, compression = comp, threads = threads)
    else anndataR::write_h5ad(obj, out, compression = comp, mode = "w")
  })[["elapsed"]]
  res$out_bytes <- file.size(out)
  unlink(out)
}
res$seconds <- t
res$peak_gb <- status("VmHWM")
res$base_gb <- base
cat("RESULT", jsonlite::toJSON(res, auto_unbox = TRUE, digits = 6), "\n")
