"""Turn bench_results.json into Markdown tables.

Usage: python summarize_bench.py <bench_results.json> > benchmarks.md
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import h5py

res = json.loads(Path(sys.argv[1]).read_text())
IMPLS = ["AnnR (16 threads)", "AnnR (1 thread)", "anndataR", "Python anndata"]


def describe(f):
    size = Path(f).stat().st_size
    with h5py.File(f, "r") as h:
        x = h["X"]
        if isinstance(x, h5py.Group):
            shape = tuple(x.attrs["shape"])
            nnz = x["data"].shape[0]
            comp = x["data"].compression or "none"
            kind = f"{x.attrs['encoding-type'].split('_')[0].upper()}, {nnz/1e6:,.0f}M nnz"
        else:
            shape = x.shape
            comp = x.compression or "none"
            kind = "dense"
    name = Path(f).parent.name if Path(f).name in ("cells.h5ad", "data.h5ad", "adata.h5ad") else Path(f).stem
    if Path(f).name == "data.h5ad":
        name = Path(f).parents[2].name
    return f"{name[:34]} ({size/1e9:.2f} GB; {shape[0]:,} x {shape[1]:,}; {kind}; {comp})"


by = defaultdict(dict)
for r in res:
    by[(r["file"], r["op"], r["cache"])][r["impl"]] = r

files = list(dict.fromkeys(r["file"] for r in res))
for op, cache, title in (("read", "warm", "Read into the native object (warm page cache)"),
                         ("read", "cold", "Read into the native object (cold: file evicted from page cache)"),
                         ("read_sce", "warm", "Read into a SingleCellExperiment (warm)"),
                         ("write_none", "warm", "Write .h5ad, no compression"),
                         ("write_gzip", "warm", "Write .h5ad, gzip")):
    impls = [i for i in IMPLS if any(i in by[(f, op, cache)] for f in files)]
    if not impls:
        continue
    print(f"\n### {title}\n")
    print("| File | " + " | ".join(impls) + " | AnnR speed-up vs anndataR | vs Python |")
    print("|---|" + "---:|" * (len(impls) + 2))
    for f in files:
        row = by[(f, op, cache)]
        if not row:
            continue
        cells = []
        for i in impls:
            r = row.get(i)
            if r is None:
                cells.append("–")
            elif "error" in r:
                cells.append("error")
            else:
                cells.append(f"{r['seconds']:.2f} s<br>{r['peak_gb']:.1f} GB")
        a = row.get("AnnR (16 threads)", {}).get("seconds")
        b = row.get("anndataR", {}).get("seconds")
        c = row.get("Python anndata", {}).get("seconds")
        su1 = f"{b / a:.1f}x" if a and b else "–"
        su2 = f"{c / a:.1f}x" if a and c else "–"
        print(f"| {describe(f)} | " + " | ".join(cells) + f" | {su1} | {su2} |")
print("\nTimes are medians of fresh-process runs; memory is the process peak RSS (VmHWM), "
      "including the R or Python runtime and loaded packages.")
