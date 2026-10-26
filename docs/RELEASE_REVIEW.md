# AnnR 0.1.0 release review

Reviewed on **26 October 2026**. Author and maintainer: **MD BABU MIA**.

## Changes made during release preparation

* Corrected the unpublished Git history to use MD BABU MIA as author and
  committer, using the GitHub account-linked noreply address. Removed the AI
  co-author trailer before first publication. Package metadata, both licenses,
  and citation files identify MD BABU MIA.
* Fixed validation of pairwise matrix dimensions and raw feature annotations.
* Removed deletion of the destination before file replacement. POSIX uses
  rename; Windows uses a backup and restoration attempt. Regression checks
  verify that serialization errors and a simulated failed final rename leave
  the existing destination intact.
* Added a safeguard against applying original names to changed cells or genes
  when exporting converted Seurat objects. Subset AnnR before conversion.
* Made validation commands return failure for mismatches or empty corpora;
  unexpected output elements also fail round-trip validation. Temporary Python
  reference directories are cleaned up, and stale result files are removed.
* Strengthened Awkward-array read comparisons from form/length checks to
  comparisons of the form and every data buffer. Added comparator negative
  controls for small value changes, index changes, and transposition.
* Added a complete illustrated tutorial, executable R examples, synthetic
  `.h5ad` and CSV data, reproducible generator, citation metadata, release
  notes, reporting guidance, and Linux GitHub Actions checks.

## Fresh verification

Environment: Linux x86_64; R 4.6.1; HDF5 1.14.6; libdeflate and Zstandard
available. Python reference: anndata 0.13.4, pandas 3.1.0, h5py 3.16.0,
NumPy 2.5.4, SciPy 1.18.1, Awkward 2.14.0, hdf5plugin 7.1.0.

| Check | Result |
|---|---|
| Source package build and `R CMD check --no-manual` | 0 errors, 0 warnings, 0 notes |
| R tests | 207 expectations passed across 23 test blocks; no failures, skips or warnings |
| Real-data read comparisons | 41/41 files; 1,534 elements; no differences |
| Strengthened synthetic read comparisons | 9/9 files; 523 elements; no differences |
| Real-data round trips, direct / SCE / Seurat | 1,222 equal elements per route, same compared dtype/encoding |
| Synthetic round trips, direct / SCE / Seurat | 448 equal elements per route, same compared dtype/encoding |
| Direct gzip / LZF / Zstandard round trips | 1,222 equal real-data elements per codec |
| Thread determinism | 50 files; read thread counts 1, 2, 3, 7, 16, 32; write 1 versus 32; no differences |
| Comparator negative controls | 4/4 passed; empty-corpus driver returns failure |
| External-buffer regression | Copy-on-modify, saveRDS/readRDS, GC lifetime, and threaded re-read passed |
| Tutorial | All non-installation R code blocks executed; AnnR, SCE and Seurat exports verified in Python |

Read comparisons total **2,057 elements**. The historical total of 2,051
predates six additional Awkward buffer checks. The source's native C++ engine
was not changed in this review. The Seurat identity safeguard was also checked
against all 41 real and 9 synthetic files using the final installed package.
Aggregate results are in [validation-results.json](validation-results.json).
Private local input lists, dataset paths, and raw logs are not distributed.
Checks of the HDF5 engine and unmodified conversion routes used the release
candidate; the final Seurat safeguard was additionally validated against both
full small-file corpora. The native engine and writer were unchanged between
those candidate and final checks.

## Limits of this review

* Large files (0.2–7.3 GB) were not repeated, respecting the small-file scope.
  Their README results belong to a previous development build.
* Historical benchmark timings and sanitizer runs were inspected but not
  repeated. The review does not establish a universal speed ranking.
* The PDF reference manual was not built (`--no-manual`). Documentation,
  examples, installed HTML tutorial, native build, and tests were checked.
* Fresh local checks were on Linux. The Windows backup path and other
  operating systems were not independently tested here.
* Numerical equality is NaN-aware. It does not imply equal HDF5 bytes, NaN
  payloads, or signed-zero bit patterns. Legacy encodings may be normalized
  to the representation written by Python anndata.
* Immediate conversion round trips are validated. Arbitrary downstream edits
  still need validation of raw snapshots, auxiliary matrices, and metadata.
* Zarr, backed mode, and streaming row/column reads are unsupported; external
  buffers can make R's own memory accounting underestimate process memory.

GitHub distribution is separate from CRAN or Bioconductor submission. This
release does not claim acceptance by either package registry.
