#!/bin/sh
# Requires pandoc. Run from the repository root.
set -eu
pandoc docs/TUTORIAL.md --standalone --embed-resources --toc --toc-depth=2 \
  --metadata title='AnnR Tutorial — MD BABU MIA' \
  --css data-raw/tutorial.css --resource-path=docs -o inst/doc/AnnR-tutorial.html
# Normalize whitespace emitted by pandoc's default highlighting stylesheet.
rendered=$(mktemp)
sed 's/[[:blank:]]*$//' inst/doc/AnnR-tutorial.html > "$rendered"
cat "$rendered" > inst/doc/AnnR-tutorial.html
rm -f "$rendered"
