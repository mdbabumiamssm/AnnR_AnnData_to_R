"""Run the read-fidelity check (r_dump.R + compare_read.py) over many files.

Usage: python run_read_corpus.py <list.txt> <work_dir> [--keep]
Each line of list.txt is a path to an .h5ad file.
"""
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    files = [l.strip() for l in Path(sys.argv[1]).read_text().splitlines() if l.strip()]
    work = Path(sys.argv[2])
    keep = "--keep" in sys.argv
    work.mkdir(parents=True, exist_ok=True)
    results = []
    for i, f in enumerate(files):
        dump = work / f"dump_{i:03d}"
        if dump.exists():
            shutil.rmtree(dump)
        t0 = time.time()
        r = subprocess.run(["Rscript", str(HERE / "r_dump.R"), f, str(dump)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            results.append({"file": f, "status": "R-ERROR", "detail": (r.stderr or r.stdout)[-2000:]})
            print(f"R-ERROR {Path(f).name}\n{(r.stderr or r.stdout)[-1500:]}", flush=True)
            continue
        out_json = work / f"result_{i:03d}.json"
        c = subprocess.run([sys.executable, "-I", str(HERE / "compare_read.py"), f, str(dump),
                            "--json", str(out_json)], capture_output=True, text=True)
        print(c.stdout.strip(), flush=True)
        if c.returncode != 0 or not out_json.exists():
            print(c.stderr[-3000:], flush=True)
            results.append({"file": f, "status": "PY-ERROR", "detail": c.stderr[-2000:]})
        else:
            res = json.loads(out_json.read_text())
            res["status"] = "PASS" if not res["issues"] else "FAIL"
            res["wall_seconds"] = time.time() - t0
            results.append(res)
        if not keep:
            shutil.rmtree(dump, ignore_errors=True)
    summary = {s: sum(r["status"] == s for r in results) for s in ("PASS", "FAIL", "R-ERROR", "PY-ERROR")}
    (work / "summary.json").write_text(json.dumps({"summary": summary, "results": results}, indent=1))
    print("\nSUMMARY", summary)
    total = sum(r.get("checked", 0) for r in results)
    print(f"elements compared: {total}")


if __name__ == "__main__":
    main()
