"""Benchmark AnnR vs anndataR vs Python anndata (fresh process per run).
Usage: python run_bench.py <files.txt> <out.json> [--reps 3] [--cold-reps 2]
"""
import json, statistics, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
WORK = HERE.parent / "work" / "bench_out"
WORK.mkdir(parents=True, exist_ok=True)

def arg(name, default):
    return int(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default

def run(cmd, timeout=3600):
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"error": "timeout"}
    for line in r.stdout.splitlines():
        if line.startswith("RESULT "):
            d = json.loads(line[7:]); d["wall"] = time.time() - t0; return d
    return {"error": (r.stderr or r.stdout)[-600:]}

def evict(f):
    subprocess.run([PY, "-I", str(HERE / "evict.py"), f], check=False)

def main():
    files = [l.strip() for l in Path(sys.argv[1]).read_text().splitlines() if l.strip() and not l.startswith("#")]
    out_json = Path(sys.argv[2])
    reps, cold_reps = arg("--reps", 3), arg("--cold-reps", 2)
    results = json.loads(out_json.read_text()) if out_json.exists() else []
    done = {(r["file"], r["impl"], r["op"], r["cache"]) for r in results}
    out = str(WORK / "bench_tmp.h5ad")
    impls = [
        ("AnnR (16 threads)", lambda op, f, c: ["Rscript", str(HERE / "bench_r.R"), "annr", op, f, out, "16", c]),
        ("AnnR (1 thread)", lambda op, f, c: ["Rscript", str(HERE / "bench_r.R"), "annr", op, f, out, "1", c]),
        ("anndataR", lambda op, f, c: ["Rscript", str(HERE / "bench_r.R"), "anndataR", op, f, out, "1", c]),
        ("Python anndata", lambda op, f, c: [PY, "-I", str(HERE / "bench_py.py"), op.replace("_sce", ""), f, out, c]),
    ]
    ops = [("read", "none"), ("read_sce", "none"), ("write", "none"), ("write", "gzip")]
    for f in files:
        for impl, mk in impls:
            for op, comp in ops:
                if impl == "Python anndata" and op == "read_sce":
                    continue
                if impl == "AnnR (1 thread)" and op == "read_sce":
                    continue
                opname = op if op != "write" else f"write_{comp}"
                for cache, n in (("warm", reps), ("cold", cold_reps if op.startswith("read") else 0)):
                    if n == 0 or (f, impl, opname, cache) in done:
                        continue
                    times, peaks, last = [], [], None
                    if cache == "warm":
                        run(mk(op, f, comp))  # warm-up
                    for _ in range(n):
                        if cache == "cold":
                            evict(f)
                        d = run(mk(op, f, comp))
                        if "error" in d:
                            last = d; break
                        times.append(d["seconds"]); peaks.append(d["peak_gb"]); last = d
                    rec = {"file": f, "impl": impl, "op": opname, "cache": cache}
                    if times:
                        rec.update(seconds=statistics.median(times), runs=times, peak_gb=max(peaks),
                                   base_gb=last.get("base_gb"), out_bytes=last.get("out_bytes"))
                    else:
                        rec["error"] = (last or {}).get("error", "no result")
                    results.append(rec)
                    out_json.write_text(json.dumps(results, indent=1))
                    print(f"{Path(f).parent.name[:28]:28s} {impl:18s} {opname:11s} {cache:4s} "
                          f"{rec.get('seconds', float('nan')):8.2f}s  peak {rec.get('peak_gb', float('nan')):6.2f} GB "
                          f"{rec.get('error', '')[:120]}", flush=True)

if __name__ == "__main__":
    main()
