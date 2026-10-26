# Reproducing AnnR validation

The R tests ship with small synthetic fixtures. Optional, independent Python
comparisons use the scripts under `dev/validation/`. They compare values,
identities, categorical levels, missing values, sparse structure, and
encoding/dtype metadata. They do not prove compatibility with every possible
future encoding or arbitrary edits in downstream packages.

## R package checks

From a checkout with AnnR's dependencies installed:

```sh
Rscript -e 'testthat::test_local(".", stop_on_failure = TRUE)'
R CMD build .
R CMD check --no-manual AnnR_0.1.0.tar.gz
Rscript -e 'source(system.file("examples", "tutorial.R", package="AnnR"))'
```

`--no-manual` skips the PDF reference manual, not tests or examples. Install
LaTeX and omit this flag to check PDF manual creation. GitHub Actions runs
R package checks on Linux for pushes and pull requests.

## Python reference environment

Use an isolated environment and record exact installed versions. Development
results used Python anndata 0.13.4, pandas 3.1 and h5py 3.16. The release review
records the actual environment used for its new comparisons.

```sh
python -m venv .oracle/venv
.oracle/venv/bin/python -m pip install anndata==0.13.4 numpy pandas scipy h5py awkward hdf5plugin
.oracle/venv/bin/python -m pip freeze
mkdir -p dev/work/syn
.oracle/venv/bin/python -I dev/validation/make_synthetic.py dev/work/syn > dev/work/syn/list.txt
.oracle/venv/bin/python -I dev/validation/run_read_corpus.py dev/work/syn/list.txt dev/work/read_syn
.oracle/venv/bin/python -I dev/validation/run_roundtrip_corpus.py dev/work/syn/list.txt dev/work/rt_syn_annr
.oracle/venv/bin/python -I dev/validation/run_roundtrip_corpus.py dev/work/syn/list.txt dev/work/rt_syn_sce --via sce
.oracle/venv/bin/python -I dev/validation/run_roundtrip_corpus.py dev/work/syn/list.txt dev/work/rt_syn_seurat --via seurat
```

These drivers use the AnnR installed in the R library visible to `Rscript`.
Install the checkout first. To target a separate test library, set `R_LIBS`
for all commands. Read and round-trip drivers now exit with failure when
comparisons fail or the input list is empty. `awkward` and `hdf5plugin` enable
optional edge cases; confirm that the expected nine files were generated.

For real files, create `dev/validation/corpus_small.txt` with one local path
per line, then run the same drivers using that list. The original development
corpus mixes public test files and local pipeline outputs, so its exact
41-file manifest is intentionally not distributed. Users can reproduce the
method on their own files; the generated fixtures are fully reproducible.
Never commit machine-specific corpus lists or confidential inputs.

Round-trip comparisons use original metadata for modern files and Python's
normalized rewrite for legacy files or elements without comparable original
metadata. They compare decoded element values and dtype/encoding, not whole
HDF5 bytes, compression ratio, chunk placement, or ordering of equivalent
sparse indices. Unmodified conversion routes are tested; changes made inside
Seurat or SingleCellExperiment require additional checks appropriate to those
changes.

`run_all_small.sh` runs both prepared corpora through all conversion routes,
plus direct gzip/LZF/Zstandard output and thread determinism. The script
requires Bash and propagates command failures. `determinism.R` uses both local
list files and exits with failure on a mismatch.

## Comparator negative controls

```sh
.oracle/venv/bin/python -I dev/validation/test_comparators.py
```

These four controls reject a 1e-7 change, an index change, and a transposed
matrix, while accepting an unchanged matrix. The read comparator also checks
Awkward forms and buffer contents; historical totals predate the extra buffer
checks. Equality is numeric and NaN-aware, not a signed-zero or NaN-payload
bit-pattern guarantee.

## Sanitizers and benchmarks

`asan_suite.R` exercises an installed sanitizer build; setting up a compiler
with AddressSanitizer and UndefinedBehaviorSanitizer is platform dependent.
Historical sanitizer results are distinguished from fresh release tests.
`dev/bench/` contains the timing harness. Compare equivalent matrices,
conversion targets, thread counts and compression; record cold/warm cache
conditions and process peak RSS. Timings are workload and machine dependent.

See [RELEASE_REVIEW.md](RELEASE_REVIEW.md) for precisely what was repeated for
this release. Large-file development results are not a fresh validation of
the final release build.
