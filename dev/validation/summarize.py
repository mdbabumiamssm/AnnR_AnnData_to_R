"""Summarise validation results (read fidelity + round trips) as Markdown.

Usage: python summarize.py <work_dir> > validation_summary.md
"""
import json
import sys
from pathlib import Path

work = Path(sys.argv[1])


def load(name):
    p = work / name / "summary.json"
    return json.loads(p.read_text()) if p.exists() else None


print("## Read fidelity (AnnR vs Python anndata 0.13.4, exact element comparison)\n")
print("| Corpus | Files | Pass | Fail | Errors | Elements compared |")
print("|---|---:|---:|---:|---:|---:|")
for label, name in (("Real files < 200 MB", "read_small"), ("Real files 0.2-7.3 GB", "read_big"),
                    ("Synthetic edge cases", "read_syn")):
    d = load(name)
    if not d:
        continue
    s = d["summary"]
    n = sum(s.values())
    el = sum(r.get("checked", 0) for r in d["results"])
    print(f"| {label} | {n} | {s['PASS']} | {s['FAIL']} | {s['R-ERROR'] + s['PY-ERROR']} | {el:,} |")

print("\n## Round trips (Python file -> AnnR -> .h5ad -> Python anndata), per element\n")
print("| Corpus | Path | Compression | Equal, same encoding & dtype | Equal, different encoding/dtype "
      "| Values differ | Python cannot read | Not written |")
print("|---|---|---|---:|---:|---:|---:|---:|")
for name in sorted(p.name for p in work.iterdir() if p.is_dir() and p.name.startswith("rt_")):
    d = load(name)
    if not d:
        continue
    t = d["total"]
    corpus = {"small": "real < 200 MB", "big": "real 0.2-7.3 GB", "syn": "synthetic"}[name.split("_")[1]]
    errors = sum(1 for r in d["results"] if r.get("status") in ("R-ERROR", "PY-ERROR"))
    print(f"| {corpus} | {d['via']} | {d['compression']} | {t['equal_same']:,} | {t['equal_diff']} | "
          f"{t['values_differ']} | {t['unreadable'] + errors} | {t['missing']} |")
