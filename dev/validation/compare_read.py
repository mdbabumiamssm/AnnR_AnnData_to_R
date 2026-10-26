"""Compare an AnnR dump (r_dump.R) against Python anndata, element by element.

The reference is `anndata.read_h5ad` from the official Python package; no code
is shared with AnnR. Every value is compared exactly (NaN == NaN), shapes and
orientation are checked, and the R type chosen for each dtype is verified.

Usage: python compare_read.py <file.h5ad> <dump_dir> [--json out.json]
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad
try:
    import hdf5plugin  # noqa: F401  (zstd/blosc filters)
except ImportError:
    pass

INT32_MIN = np.iinfo(np.int32).min


class Checker:
    def __init__(self, dump_dir: Path):
        self.dir = dump_dir
        self.issues: list[str] = []
        self.n_checked = 0

    # -- loading helpers -------------------------------------------------
    def npy(self, name):
        return np.load(self.dir / name, mmap_mode="r")

    def strings(self, name):
        with open(self.dir / name, encoding="utf-8") as fh:
            return json.load(fh)

    def fail(self, path, msg):
        self.issues.append(f"{path}: {msg}")

    # -- comparisons -----------------------------------------------------
    def same_numeric(self, path, got, exp):
        """Exact comparison of float64/int arrays, NaN-aware, chunked."""
        got = np.asarray(got)
        exp = np.asarray(exp)
        if got.shape != exp.shape:
            self.fail(path, f"shape {got.shape} != {exp.shape}")
            return False
        n = got.size
        step = 50_000_000
        g_flat = got.reshape(-1)
        e_flat = exp.reshape(-1)
        for s in range(0, max(n, 1), step):
            g = np.asarray(g_flat[s : s + step])
            e = np.asarray(e_flat[s : s + step])
            if g.dtype.kind == "f" or e.dtype.kind == "f":
                g64 = g.astype(np.float64)
                e64 = e.astype(np.float64)
                ok = (g64 == e64) | (np.isnan(g64) & np.isnan(e64))
            else:
                ok = g.astype(np.int64) == e.astype(np.int64)
            if not ok.all():
                bad = np.flatnonzero(~ok)[:5] + s
                self.fail(path, f"{(~ok).sum()} values differ, e.g. at {bad.tolist()}: "
                                f"got {g_flat[bad].tolist()} expected {e_flat[bad].tolist()}")
                return False
        return True

    def expected_rtype(self, arr: np.ndarray) -> str:
        k = arr.dtype.kind
        if k == "b":
            return "logical"
        if k == "f":
            return "double"
        if k == "c":
            return "complex"
        if k in "iu":
            small = (k == "i" and arr.dtype.itemsize <= 4) or (k == "u" and arr.dtype.itemsize <= 2)
            if small:
                return "integer"
            if arr.size == 0:
                return "integer"
            a = np.asarray(arr)
            lo, hi = a.min(), a.max()
            if k == "u":
                return "integer" if int(hi) <= 2**31 - 1 else "double"
            return "integer" if (int(lo) >= -(2**31 - 1) and int(hi) <= 2**31 - 1) else "double"
        return "?"

    def check_array(self, path, node, arr):
        self.n_checked += 1
        if node["kind"] != "array":
            self.fail(path, f"expected array node, got {node['kind']}")
            return
        arr = np.asarray(arr)
        if arr.dtype.kind in "OUST":
            if node["rtype"] != "character":
                self.fail(path, f"string array read as {node['rtype']}")
                return
            got = self.strings(node["data"])
            exp = [x.decode("utf-8") if isinstance(x, bytes) else (None if x is None else str(x))
                   for x in arr.reshape(-1, order="F")]
            if node["dim"] is not None and list(node["dim"]) != list(arr.shape):
                self.fail(path, f"dim {node['dim']} != shape {arr.shape}")
            if got != exp:
                diffs = [i for i, (g, e) in enumerate(zip(got, exp)) if g != e][:3]
                self.fail(path, f"strings differ (n={len(got)} vs {len(exp)}), e.g. {[(got[i], exp[i]) for i in diffs]}")
            return
        exp_type = self.expected_rtype(arr)
        if node["rtype"] != exp_type:
            self.fail(path, f"dtype {arr.dtype} read as {node['rtype']}, expected {exp_type}")
            return
        if arr.ndim >= 2:
            if node["dim"] is None or list(node["dim"]) != list(arr.shape):
                self.fail(path, f"dim {node['dim']} != shape {arr.shape}")
                return
        elif node["length"] != arr.size:
            self.fail(path, f"length {node['length']} != {arr.size}")
            return
        got = self.npy(node["data"])
        if arr.ndim >= 2:
            got = np.asarray(got).reshape(arr.shape, order="F")
        else:
            got = np.asarray(got).reshape(arr.shape)
        if arr.dtype.kind == "c":
            self.same_numeric(path + "[re]", got.real, arr.real)
            self.same_numeric(path + "[im]", got.imag, arr.imag)
            return
        exp = arr.astype(np.int64) if arr.dtype.kind == "b" else arr
        self.same_numeric(path, got, exp)

    def check_masked(self, path, node, values, mask):
        """Nullable arrays: values where mask is False, NA elsewhere."""
        self.n_checked += 1
        if node["kind"] != "array":
            self.fail(path, f"expected array node for nullable, got {node['kind']}")
            return
        if node["rtype"] == "character":
            got = self.strings(node["data"])
            exp = [None if m else str(v) for v, m in zip(values, mask)]
            if got != exp:
                self.fail(path, "nullable strings differ")
            return
        got = np.asarray(self.npy(node["data"]))
        if node["rtype"] in ("integer", "logical"):
            na = got == INT32_MIN
        else:
            na = np.isnan(got)
        if not np.array_equal(na, np.asarray(mask, dtype=bool)):
            self.fail(path, "NA positions differ from mask")
            return
        keep = ~np.asarray(mask, dtype=bool)
        v = np.asarray(values)
        self.same_numeric(path, got[keep], v[keep].astype(np.int64) if v.dtype.kind == "b" else v[keep])

    def check_sparse(self, path, node, m):
        self.n_checked += 1
        if node["kind"] != "sparse":
            self.fail(path, f"expected sparse node, got {node['kind']}")
            return
        fmt = m.format
        is_bool = m.dtype == bool
        want = {"csr": "lgRMatrix" if is_bool else "dgRMatrix",
                "csc": "lgCMatrix" if is_bool else "dgCMatrix"}[fmt]
        if node["class"] != want:
            self.fail(path, f"class {node['class']} != {want}")
            return
        if list(node["dim"]) != list(m.shape):
            self.fail(path, f"dim {node['dim']} != {m.shape}")
            return
        c = m.copy()
        c.sum_duplicates()  # also sorts indices, exactly what AnnR guarantees
        self.same_numeric(path + "@p", self.npy(node["p"]), c.indptr)
        self.same_numeric(path + "@idx", self.npy(node["idx"]), c.indices)
        data = c.data.astype(np.int64) if is_bool else c.data
        self.same_numeric(path + "@x", self.npy(node["x"]), data)

    def check_factor(self, path, node, cat):
        self.n_checked += 1
        if node["kind"] != "factor":
            self.fail(path, f"categorical read as {node['kind']}")
            return
        levels = self.strings(node["levels"])
        exp_levels = [str(c) if not isinstance(c, (float, np.floating)) else None for c in cat.categories]
        if any(l is None for l in exp_levels):
            exp_levels = [repr_r_double(c) for c in cat.categories]
        if cat.categories.dtype.kind == "b":
            exp_levels = ["TRUE" if c else "FALSE" for c in cat.categories]
        if levels != exp_levels:
            self.fail(path, f"levels {levels[:5]} != {exp_levels[:5]}")
        if bool(node["ordered"]) != bool(cat.ordered):
            self.fail(path, "ordered flag differs")
        codes = np.asarray(self.npy(node["codes"])).astype(np.int64)
        exp = np.asarray(cat.codes).astype(np.int64) + 1
        exp[np.asarray(cat.codes) < 0] = INT32_MIN
        self.same_numeric(path + "@codes", codes, exp)

    def check_dataframe(self, path, node, df: pd.DataFrame):
        self.n_checked += 1
        if node["kind"] != "dataframe":
            self.fail(path, f"dataframe read as {node['kind']}")
            return
        rn = self.strings(node["rownames"])
        exp_rn = [str(x) for x in df.index]
        if node["nrow"] != len(df):
            self.fail(path, f"nrow {node['nrow']} != {len(df)}")
        if rn != exp_rn and not (len(df) == 0):
            self.fail(path, f"row names differ, e.g. {rn[:3]} vs {exp_rn[:3]}")
        cols = self.strings(node["colnames"])
        if cols != [str(c) for c in df.columns]:
            self.fail(path, f"columns {cols} != {list(df.columns)}")
            return
        for c, cnode in zip(df.columns, node["columns"]):
            self.check_value(f"{path}/{c}", cnode, df[c])

    def check_mapping(self, path, node, mapping):
        self.n_checked += 1
        if node["kind"] != "dict":
            self.fail(path, f"mapping read as {node['kind']}")
            return
        keys = self.strings(node["keys"])
        exp_keys = list(mapping.keys())
        if keys != exp_keys:
            if sorted(keys) == sorted(exp_keys):
                self.fail(path, f"key order differs: {keys[:6]} vs {exp_keys[:6]}")
            else:
                self.fail(path, f"keys {keys} != {exp_keys}")
                return
        lookup = dict(zip(keys, node["items"]))
        for k in exp_keys:
            self.check_value(f"{path}/{k}", lookup[k], mapping[k])

    def check_value(self, path, node, val):
        if val is None:
            self.n_checked += 1
            if node["kind"] != "null":
                self.fail(path, f"None read as {node['kind']}")
            return
        if isinstance(val, ad.AnnData):
            self.check_anndata(path, node, val)
        elif isinstance(val, pd.DataFrame):
            self.check_dataframe(path, node, val)
        elif isinstance(val, pd.Series):
            self.check_value(path, node, val.array)
        elif isinstance(val, pd.Categorical) or isinstance(getattr(val, "dtype", None), pd.CategoricalDtype):
            self.check_factor(path, node, pd.Categorical(val))
        elif isinstance(val, (pd.arrays.IntegerArray, pd.arrays.BooleanArray)):
            self.check_masked(path, node, val._data, val._mask)
        elif isinstance(val, pd.api.extensions.ExtensionArray) and isinstance(val.dtype, pd.StringDtype):
            mask = np.asarray(pd.isna(val))
            vals = np.asarray(val.astype(object))
            self.check_masked(path, node, vals, mask)
        elif isinstance(val, pd.api.extensions.ExtensionArray) and val.dtype.kind in "biuf":
            self.check_array(path, node, np.asarray(val))
        elif isinstance(val, pd.api.extensions.ExtensionArray):
            self.check_array(path, node, np.asarray(val, dtype=object))
        elif sp.issparse(val):
            self.check_sparse(path, node, val)
        elif isinstance(val, np.ndarray) and val.dtype.names is not None:
            self.check_recarray(path, node, val)
        elif isinstance(val, np.ndarray):
            self.check_array(path, node, val)
        elif isinstance(val, dict) or hasattr(val, "keys"):
            self.check_mapping(path, node, val)
        elif isinstance(val, (str, bytes)):
            self.check_array(path, node, np.array([val.decode() if isinstance(val, bytes) else val], dtype=object))
        elif isinstance(val, (bool, np.bool_)):
            self.check_array(path, node, np.array([val], dtype=bool))
        elif isinstance(val, np.integer):
            self.check_array(path, node, np.array([val], dtype=val.dtype))
        elif isinstance(val, int):
            self.check_array(path, node, np.array([val], dtype=np.uint64 if val > 2**63 - 1 else np.int64))
        elif isinstance(val, (float, np.floating)):
            self.check_array(path, node, np.array([val], dtype=np.asarray(val).dtype if isinstance(val, np.floating) else np.float64))
        elif isinstance(val, (complex, np.complexfloating)):
            self.check_array(path, node, np.array([val]))
        elif type(val).__module__.startswith("awkward"):
            self.n_checked += 1
            if node["kind"] != "awkward":
                self.fail(path, f"awkward array read as {node['kind']}")
            else:
                import awkward as ak
                form, length, buffers = ak.to_buffers(val)
                if int(node["length"]) != length:
                    self.fail(path, "awkward length differs")
                if json.loads(node["form"]) != json.loads(form.to_json()):
                    self.fail(path, "awkward form differs")
                got_buffers = node["buffers"]
                if set(got_buffers) != set(buffers):
                    self.fail(path, "awkward buffer names differ")
                else:
                    for name, values in buffers.items():
                        self.check_array(path + "/" + name, got_buffers[name], values)
        else:
            self.fail(path, f"no comparison rule for {type(val)}")

    def check_recarray(self, path, node, arr):
        self.n_checked += 1
        if node["kind"] != "dataframe":
            self.fail(path, f"rec-array read as {node['kind']}")
            return
        cols = self.strings(node["colnames"])
        if cols != list(arr.dtype.names):
            self.fail(path, f"fields {cols} != {arr.dtype.names}")
            return
        for name, cnode in zip(arr.dtype.names, node["columns"]):
            col = arr[name]
            if col.dtype.kind == "O" or col.dtype.kind == "S":
                col = np.array([x.decode() if isinstance(x, bytes) else x for x in col], dtype=object)
            self.check_value(f"{path}/{name}", cnode, np.asarray(col))

    def check_anndata(self, path, node, a: ad.AnnData):
        self.n_checked += 1
        if node["kind"] != "anndata":
            self.fail(path, f"AnnData read as {node['kind']}")
            return
        s = node["slots"]
        p = (path + "/") if path else ""
        self.check_value(p + "X", s["X"], a.X)
        self.check_dataframe(p + "obs", s["obs"], a.obs)
        self.check_dataframe(p + "var", s["var"], a.var)
        for k in ("obsm", "varm", "obsp", "varp", "layers"):
            m = {kk: vv for kk, vv in dict(getattr(a, k)).items() if kk is not None}
            self.check_mapping(p + k, s[k], m)
        self.check_mapping(p + "uns", s["uns"], dict(a.uns))
        if a.raw is None:
            if s["raw"]["kind"] != "null":
                self.fail(p + "raw", "R has raw but Python does not")
        else:
            rn = s["raw"]
            if rn["kind"] != "dict":
                self.fail(p + "raw", f"raw read as {rn['kind']}")
                return
            keys = self.strings(rn["keys"])
            items = dict(zip(keys, rn["items"]))
            self.check_value(p + "raw/X", items["X"], a.raw.X)
            self.check_dataframe(p + "raw/var", items["var"], a.raw.var)
            self.check_mapping(p + "raw/varm", items["varm"], dict(a.raw.varm))


def repr_r_double(x):
    """as.character() of a double in R (15 significant digits)."""
    v = float(x)
    if v.is_integer() and abs(v) < 1e15:
        return str(int(v))
    return format(float(f"{v:.15g}"), ".15g")


def main():
    h5ad_path, dump_dir = sys.argv[1], Path(sys.argv[2])
    out_json = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None
    manifest = json.loads((dump_dir / "manifest.json").read_text(encoding="utf-8"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = ad.read_h5ad(h5ad_path)
    chk = Checker(dump_dir)
    chk.check_anndata("", manifest["root"], a)
    result = {"file": h5ad_path, "checked": chk.n_checked, "issues": chk.issues,
              "r_read_seconds": manifest["read_seconds"], "r_warnings": manifest["warnings"]}
    status = "PASS" if not chk.issues else "FAIL"
    result["direct_reads"] = manifest.get("direct_reads")
    result["fallback_reads"] = manifest.get("fallback_reads")
    result["fallback_notes"] = manifest.get("fallback_notes")
    print(f"{status} {Path(h5ad_path).name}: {chk.n_checked} elements checked, {len(chk.issues)} issues "
          f"[direct={manifest.get('direct_reads')} fallback={manifest.get('fallback_reads')} "
          f"{manifest.get('fallback_notes') or ''}]")
    for i in chk.issues[:25]:
        print("   -", i)
    if out_json:
        Path(out_json).write_text(json.dumps(result, indent=1))
    sys.exit(1 if chk.issues else 0)


if __name__ == "__main__":
    main()
