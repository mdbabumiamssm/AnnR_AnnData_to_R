"""Round-trip every file through AnnR (read -> write) and classify elements.

Usage: python run_roundtrip_corpus.py <list.txt> <work_dir> [--via annr|sce|seurat] [--compression none|gzip|lzf|zstd]
"""
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


def main():
    files = [l.strip() for l in Path(sys.argv[1]).read_text().splitlines() if l.strip()]
    work = Path(sys.argv[2])
    via = arg("--via", "annr")
    comp = arg("--compression", "none")
    work.mkdir(parents=True, exist_ok=True)
    total = {k: 0 for k in ("equal_same", "equal_diff", "values_differ", "unreadable", "missing")}
    rows = []
    for i, f in enumerate(files):
        out = work / f"rt_{i:03d}.h5ad"
        t0 = time.time()
        r = subprocess.run(["Rscript", str(HERE / "r_roundtrip.R"), f, str(out), via, comp],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"R-ERROR {Path(f).name}: {(r.stderr or r.stdout)[-800:]}", flush=True)
            rows.append({"file": f, "status": "R-ERROR", "detail": (r.stderr or r.stdout)[-2000:]})
            continue
        js = work / f"rt_{i:03d}.json"
        js.unlink(missing_ok=True)
        c = subprocess.run([sys.executable, "-I", str(HERE / "compare_roundtrip.py"), f, str(out),
                            "--json", str(js)], capture_output=True, text=True)
        if not js.exists():
            print(f"PY-ERROR {Path(f).name}: {c.stderr[-1500:]}", flush=True)
            rows.append({"file": f, "status": "PY-ERROR", "detail": c.stderr[-2000:]})
            continue
        res = json.loads(js.read_text())
        for k, v in res["counts"].items():
            total[k] += v
        bad = {k: v for k, v in res["counts"].items() if k != "equal_same" and v}
        if res.get("extra_elements"):
            bad["extra_elements"] = len(res["extra_elements"])
        status = "OK" if not bad else "DIFF"
        print(f"{status:4s} {Path(f).parent.name}/{Path(f).name}  {res['counts']}  ({r.stdout.strip()})", flush=True)
        for d in res["details"][:8]:
            print("      ", d, flush=True)
        rows.append({"file": f, "status": status, "counts": res["counts"], "details": res["details"],
                     "timing": r.stdout.strip(), "wall": time.time() - t0})
        out.unlink(missing_ok=True)
    (work / "summary.json").write_text(json.dumps({"via": via, "compression": comp, "total": total,
                                                   "results": rows}, indent=1))
    print("\nTOTAL", via, comp, total)
    sys.exit(1 if any(r["status"] != "OK" for r in rows) or not rows else 0)


if __name__ == "__main__":
    main()
