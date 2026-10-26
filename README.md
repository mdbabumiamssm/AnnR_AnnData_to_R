# AnnR

**Fast, faithful AnnData (`.h5ad`) input/output for R.**

AnnR reads and writes the on-disk format of Python's
[`anndata`](https://anndata.readthedocs.io) through a native C++ engine linked
directly against the HDF5 C library. It was built to be the most efficient and
the most faithful way to move single-cell data between Python and R:

* **Parallel direct-chunk reader.** HDF5 is used only for metadata; worker
  threads `pread()` the raw chunks, decompress them (gzip via
  [libdeflate](https://github.com/ebiggers/libdeflate), LZF, Zstandard) and
  convert/transpose them straight into R memory. Anything outside this
  envelope (big-endian data, fletcher32, scale-offset, unallocated chunks,
  user blocks, plugin filters) transparently falls back to `H5Dread`.
* **Zero-copy sparse matrices.** CSR becomes `dgRMatrix` and CSC becomes
  `dgCMatrix` with the on-disk arrays used as-is; conversion to the
  features x cells orientation of `SingleCellExperiment`/`Seurat` is a free
  reinterpretation (CSR of cells x genes == CSC of genes x cells).
* **Lossless dtypes.** R has no `float32`, `int8`, `uint64`, ... AnnR records
  the on-disk dtype and encoding of every element and `write_h5ad()` restores
  them whenever that is lossless (values are checked), so a round trip
  reproduces the original file's dtypes and encodings. 64-bit integers beyond
  2^53 keep their exact bit patterns.
* **Every anndata encoding**, current and legacy: `csr_matrix`, `csc_matrix`,
  dense arrays of any rank, `dataframe` (0.1.0 and 0.2.0), `categorical`,
  `nullable-integer`/`-boolean`/`-string-array`, `string`, `numeric-scalar`,
  `rec-array`, `awkward-array`, `null`, nested `AnnData`, plus files written
  by anndata < 0.7 (compound `obs`/`var`, `h5sparse` matrices, `raw.X`,
  `uns` categories, neighbour graphs in `uns`).
* **Converters** to and from `SingleCellExperiment` and `Seurat` (v5) that
  are themselves lossless (validated by round trips, see below).

## Installation

System requirements: the HDF5 C library (>= 1.10.5) and zlib development
files; libdeflate and libzstd are used when available.

```sh
# Debian/Ubuntu
sudo apt install libhdf5-dev zlib1g-dev libdeflate-dev libzstd-dev
# Fedora:  sudo dnf install hdf5-devel zlib-devel libdeflate-devel libzstd-devel
# macOS:   brew install hdf5 libdeflate zstd
```

```r
# install.packages("remotes")
remotes::install_github("mdbabumiamssm/AnnR")
```

If no system HDF5 is found, the configure script falls back to the
Bioconductor package `Rhdf5lib`.

## Usage

```r
library(AnnR)

ad <- read_h5ad("pbmc68k_reduced.h5ad") # AnnR object, cells x genes like Python
ad
#> AnnR object with n_obs x n_vars = 700 x 765
#>     obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain'
#>     var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
#>     uns: 'bulk_labels_colors', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'PCs'
#>     obsp: 'connectivities', 'distances'
#>     raw: 765 vars
#>     X: matrix 700 x 765 (stored as float32)

ad$X; ad$obs; ad$layers$counts; ad$obsm$X_umap; ad$uns
ad[ad$obs$louvain == "1", ]             # subsets every slot consistently

sce <- read_h5ad("pbmc.h5ad", as = "SingleCellExperiment")
seu <- read_h5ad("pbmc.h5ad", as = "Seurat")

# read only what you need
small <- read_h5ad("big.h5ad", layers = "counts", obsm = "X_umap", raw = FALSE, uns = FALSE)
umap  <- read_h5ad_elem("big.h5ad", "obsm/X_umap")
h5ad_info("big.h5ad")                   # structure, dtypes, chunking, compression

write_h5ad(ad,  "out.h5ad")             # AnnR, SingleCellExperiment or Seurat
write_h5ad(sce, "out.h5ad", compression = "gzip")
```

The number of threads defaults to the available cores (max 16) and can be
set per call (`threads =`) or with `options(AnnR.threads = n)`.

### Type mapping

| anndata / numpy                       | R                                            |
|---------------------------------------|----------------------------------------------|
| `float16/32/64`                       | `double`                                     |
| `int8/16/32`, `uint8/16`              | `integer`                                    |
| `int64`, `uint32/64`                  | `integer` if all values fit, else `double` (exact bits kept) |
| `bool`                                | `logical`                                    |
| strings                               | `character` (UTF-8)                          |
| `categorical`                         | `factor` (`ordered` when ordered)            |
| nullable integer/boolean/string       | vector with `NA`                             |
| `csr_matrix` / `csc_matrix`           | `dgRMatrix` / `dgCMatrix` (`lgRMatrix`/`lgCMatrix` for bool) |
| `dataframe`                           | `data.frame` (row names = index)             |
| `dict`                                | named `list`                                 |
| `rec-array`                           | `data.frame`                                 |
| complex                               | `complex`                                    |

When writing, R vectors with `NA` use the nullable encodings; a length-one
vector in `uns` is written as a scalar unless it was read as an array or is
wrapped in `I()`.

### Conversions

`as_SingleCellExperiment()`: `X` -> assay `"X"`, layers -> assays, `obs`/`var`
-> `colData`/`rowData`, `obsm` -> `reducedDims`, `varm` -> `rowData(sce)$varm`,
`obsp`/`varp` -> `colPairs`/`rowPairs`, `uns` -> `metadata()`, `raw` -> altExp
`"raw"`.

`as_Seurat()`: `X` -> layer `data`, `layers$counts` -> layer `counts`, other
layers -> layers of the same name, `obs` -> meta data, 2-D `obsm` -> reductions
(`X_pca` -> `pca`), `obsp` -> graphs, `raw` -> assay `"raw"`, `uns` -> `misc`.
Names Seurat cannot hold (underscores, duplicates) are adjusted and the
originals restored by `write_h5ad()`.

`as_AnnR()` converts `SingleCellExperiment`, `SummarizedExperiment` and
`Seurat` objects back.

## Validation

Every claim below was checked against the reference implementation, Python
`anndata` 0.13.4 (pandas 3.1, h5py 3.16), with scripts in `dev/validation/`.

**Read fidelity.** AnnR's in-memory result is dumped and compared, element by
element and bit for bit (NaN-aware), with what Python `anndata.read_h5ad`
returns, including R type choice, orientation, factor levels/codes, row names
and NA positions. The comparator shares no code with AnnR and was itself
mutation-tested (a single `1e-7` change, a transposed embedding or one moved
sparse index each make it fail).

| Corpus | Files | Pass | Elements compared |
|---|---:|---:|---:|
| Real-world files < 200 MB (scanpy, scvi-tools, squidpy and anndata test data, archives written by anndata 0.5-0.12, 10x/CITE-seq pipeline outputs) | 41 | 41 | 1,534 |
| Synthetic edge cases (all dtypes and encodings, NaN/Inf/-0.0, unicode, nullable, int64/uint64 beyond 2^53, rec-arrays, awkward, nested AnnData, big-endian, fletcher32, scale-offset, unallocated chunks, user block, zstd, LZF, 2-D chunking) | 9 | 9 | 517 |
| Real-world files 0.2-7.3 GB (up to 203,496 cells and 2 x 484M non-zeros) * | 11 | 11 | 681 |

**Round trips.** Each file is read by AnnR, written back, read by Python and
classified per element (same categories as the
[anndataR comparison post](https://lazappi.id.au/posts/2026-10-25-anndata-r-claude/)).
On-disk encoding and dtype are compared with the original file.

| Corpus | Route | Compression | Equal, same encoding & dtype | Equal, other encoding/dtype | Values differ | Python cannot read | Not written |
|---|---|---|---:|---:|---:|---:|---:|
| real < 200 MB | AnnR | none | 1,222 | 0 | 0 | 0 | 0 |
| real < 200 MB | AnnR | gzip | 1,222 | 0 | 0 | 0 | 0 |
| real < 200 MB | AnnR | lzf | 1,222 | 0 | 0 | 0 | 0 |
| real < 200 MB | AnnR | zstd | 1,222 | 0 | 0 | 0 | 0 |
| real < 200 MB | via SingleCellExperiment | none | 1,222 | 0 | 0 | 0 | 0 |
| real < 200 MB | via Seurat | none | 1,222 | 0 | 0 | 0 | 0 |
| synthetic | AnnR / SCE / Seurat | none | 448 each | 0 | 0 | 0 | 0 |
| real 0.2-7.3 GB * | AnnR | none | 561 | 0 | 0 | 0 | 0 |

\* The large files were validated with the previous build; the final build
changed only how large vectors are allocated (see below) and was re-validated
on the small and synthetic corpora.

**Robustness.**

* AddressSanitizer + UndefinedBehaviorSanitizer: 50 files x (read, write with
  4 codecs, read back, SCE and Seurat conversion, subsetting): no findings.
* Determinism: results are bit-identical with 1, 2, 3, 7, 16 and 32 threads
  (50 files), and files written with 1 or 32 threads read back identically.
* Truncated, corrupted (damaged gzip chunk) and non-HDF5 files produce R
  errors naming the element, never crashes; failed writes leave no file.
* `R CMD check`: 0 errors, 0 warnings, 0 notes; 181 testthat expectations.


## Benchmarks

Fresh process per run, median of 3 (warm page cache) or 2 (cold: file evicted
from the page cache). AnnR 0.1.0 (16 threads, and 1 thread), anndataR 1.2.3
(Bioconductor), Python anndata 0.13.4. 32-core Linux workstation, NVMe.
Memory is the peak RSS of the whole process, including the R/Python runtime
and loaded packages (~0.5 GB for R with Bioconductor).

Datasets: **A** scanpy `pbmc68k_reduced` (700 x 765, dense, gzip);
**B** 1,500 x 33,538 cells x genes, CSR 3.5M non-zeros, gzip (10 MB);
**C** scvi-hub Tabula Sapiens bone marrow, 27,112 x 3,000, metadata-heavy (32 MB);
**D** 18,786 x 20,686, two CSR matrices of 31M non-zeros with unsorted
indices, gzip (171 MB).

Read into the native object (seconds, cold cache):

| Dataset | AnnR | AnnR, 1 thread | anndataR | Python anndata | vs anndataR | vs Python |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.02 | 0.02 | 0.39 | 0.04 | 19x | 2.1x |
| B | 0.04 | 0.07 | 0.71 | 0.11 | 16x | 2.4x |
| C | 0.07 | 0.07 | 1.08 | 0.11 | 15x | 1.5x |
| D | 0.20 | 0.91 | 5.13 | 1.25 | 25x | 6.2x |

Read into a `SingleCellExperiment` (seconds, warm):

| Dataset | AnnR | anndataR | speed-up |
|---|---:|---:|---:|
| A | 0.14 | 0.84 | 6x |
| B | 0.08 | 1.32 | 16x |
| C | 0.40 | 12.90 | 33x |
| D | 0.17 | 11.88 | 69x |

Write (seconds):

| Dataset | AnnR | anndataR | Python | AnnR gzip | anndataR gzip | Python gzip |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.02 | 0.72 | 0.02 | 0.03 | 0.81 | 0.06 |
| B | 0.02 | 0.54 | 0.04 | 0.03 | 1.01 | 0.32 |
| C | 0.06 | 1.35 | 0.10 | 0.16 | 2.02 | 0.47 |
| D | 0.15 | 1.01 | 0.15 | 0.28 | 10.79 | 5.63 |

Peak memory for D: AnnR 1.2 GB, anndataR 1.8 GB, Python 0.6 GB (R keeps
`float32` values as `double`).

With the previous build, a 7.3 GB uncompressed file (203,496 cells, two
matrices of 484M non-zeros) was read in about 2 s and written in about 2.3 s
(warm cache).

### Why it is fast

* Chunks are fetched with `pread()` and decoded by all cores; HDF5 is only
  asked where the chunks are (`H5Dchunk_iter`).
* Values are converted (float32 -> double, int8 codes -> factors, ...) and,
  for dense matrices, transposed in cache-sized tiles directly into the
  destination; there is no intermediate copy.
* Sparse matrices are not transposed: CSR is a `dgRMatrix` as-is, and the
  `SingleCellExperiment`/`Seurat` orientation is a reinterpretation.
* Unsorted sparse indices (common after `anndata.concat`) are sorted with a
  per-row radix sort, which `Matrix` requires.
* Large vectors are allocated outside R's managed heap (ALTREP), so a read
  does not trigger repeated full garbage collections. In a session with
  Bioconductor loaded those collections were 92% of the read time. Memory
  is released by the garbage collector as usual; `read_h5ad()` runs a
  collection itself when more than 2 GB of such buffers are alive.
* Writes compress chunks in parallel and hand them to HDF5 with
  `H5Dwrite_chunk`.


## Limitations

* `.h5ad` (HDF5) only; Zarr stores are not supported yet.
* Objects are read into memory; there is no backed/lazy mode.
* R's `Matrix` package limits sparse matrices to 2^31 - 1 non-zeros.
* R stores `float32` data as `double`, so in-memory matrices take more space
  than in Python (the file is written back as `float32`).
* Blosc-compressed datasets need the HDF5 blosc plugin (`HDF5_PLUGIN_PATH`).
  gzip, LZF and Zstandard are built in.
* On Windows the parallel direct path is disabled (reads use `H5Dread`).
* Large vectors live outside R's managed heap, so `object.size()` and
  `gc()` under-report them.

## Reproducing the validation

`dev/validation/` contains everything used to validate AnnR against Python
`anndata` 0.13.4: `make_synthetic.py` (edge-case files), `r_dump.R` +
`compare_read.py` (exact element-wise read comparison),
`compare_roundtrip.py` (per-element round-trip classification against the
original file), `run_*_corpus.py` drivers, sanitizer and determinism suites.
`dev/bench/` contains the benchmark harness.

## License

MIT.
