#!/bin/sh
# Full validation on the small real files + synthetic edge cases.
set -u
cd "$(dirname "$0")/../.."
PY=.oracle/venv/bin/python
W=dev/work
$PY -I dev/validation/make_fixtures.py tests/testthat/fixtures > /dev/null
$PY -I dev/validation/run_read_corpus.py dev/validation/corpus_small.txt $W/read_small | tail -n 2
$PY -I dev/validation/run_read_corpus.py $W/syn/list.txt $W/read_syn | tail -n 2
for via in annr sce seurat; do
  $PY -I dev/validation/run_roundtrip_corpus.py dev/validation/corpus_small.txt $W/rt_small_$via --via $via | tail -n 1
  $PY -I dev/validation/run_roundtrip_corpus.py $W/syn/list.txt $W/rt_syn_$via --via $via | tail -n 1
done
for c in gzip lzf zstd; do
  $PY -I dev/validation/run_roundtrip_corpus.py dev/validation/corpus_small.txt $W/rt_small_$c --compression $c | tail -n 1
done
Rscript dev/validation/determinism.R 2>&1 | tail -n 1
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = FALSE)' 2>&1 | tail -n 15
