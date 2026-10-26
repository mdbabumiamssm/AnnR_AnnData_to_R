#!/usr/bin/env bash
# Run from a checkout after preparing the local corpus lists; see docs/VALIDATION.md.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${ANNR_PYTHON:-.oracle/venv/bin/python}
W=dev/work
"$PY" -I dev/validation/make_fixtures.py tests/testthat/fixtures >/dev/null
"$PY" -I dev/validation/run_read_corpus.py dev/validation/corpus_small.txt "$W/read_small"
"$PY" -I dev/validation/run_read_corpus.py "$W/syn/list.txt" "$W/read_syn"
for via in annr sce seurat; do
  "$PY" -I dev/validation/run_roundtrip_corpus.py dev/validation/corpus_small.txt "$W/rt_small_$via" --via "$via"
  "$PY" -I dev/validation/run_roundtrip_corpus.py "$W/syn/list.txt" "$W/rt_syn_$via" --via "$via"
done
for c in gzip lzf zstd; do
  "$PY" -I dev/validation/run_roundtrip_corpus.py dev/validation/corpus_small.txt "$W/rt_small_$c" --compression "$c"
done
Rscript dev/validation/determinism.R
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)'
