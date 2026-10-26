#!/bin/sh
# Clean, optimised (-O2) build + install of AnnR. Regenerates cpp11 bindings
# and roxygen docs first; roxygen compiles with -O0 in src/, so objects are
# removed before the real install.
set -e
cd "$(dirname "$0")/.."
Rscript -e 'cpp11::cpp_register()' >/dev/null
if [ "${1:-}" != "--no-docs" ]; then
  Rscript -e 'suppressMessages(roxygen2::roxygenise())' >/dev/null 2>&1 || Rscript -e 'roxygen2::roxygenise()'
fi
rm -f src/*.o src/*.so
MAKEFLAGS="-j16" R CMD INSTALL --preclean --no-test-load . > dev/logs/install.log 2>&1 || { tail -50 dev/logs/install.log; exit 1; }
grep -E "warning:" dev/logs/install.log || true
rm -f src/*.o src/*.so
tail -1 dev/logs/install.log
