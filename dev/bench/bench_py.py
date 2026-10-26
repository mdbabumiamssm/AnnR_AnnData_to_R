"""One benchmark measurement with Python anndata (reference implementation).
Usage: python bench_py.py <op> <file> <out> <compression>
"""
import json, sys, time, warnings
op, file, out, comp = sys.argv[1:5]

def status(key):
    for line in open("/proc/self/status"):
        if line.startswith(key):
            return int(line.split()[1]) / 1024**2
import anndata as ad
base = status("VmRSS")
warnings.simplefilter("ignore")
res = {"tool": "python-anndata", "op": op, "file": file, "compression": comp}
if op == "read":
    t0 = time.perf_counter(); a = ad.read_h5ad(file); res["seconds"] = time.perf_counter() - t0
else:
    a = ad.read_h5ad(file)
    t0 = time.perf_counter()
    a.write_h5ad(out, compression=None if comp == "none" else comp)
    res["seconds"] = time.perf_counter() - t0
    import os
    res["out_bytes"] = os.path.getsize(out); os.unlink(out)
res["peak_gb"] = status("VmHWM"); res["base_gb"] = base
print("RESULT", json.dumps(res))
