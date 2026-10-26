"""Classify every element of an R-written .h5ad against Python anndata.

Reference on-disk metadata: the original file itself (modern files), or what
Python anndata writes for it (legacy files, `py_ref`). For every leaf element (each X/layer/obsm/... entry, every
obs/var column and index, every uns leaf, raw parts) the R-written file is
classified like the anndataR blog comparison:

  equal_same        values equal, same encoding-type and dtype as Python's own write
  equal_diff        values equal, but encoding or dtype differ
  values_differ     Python reads it but the values differ
  unreadable        Python cannot read the file
  missing           R did not write the element

Usage: python compare_roundtrip.py <original.h5ad> <r_written.h5ad> [--json out.json] [--keep-ref]
"""
from __future__ import annotations

import json
import sys
import tempfile
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad
try:
    import hdf5plugin  # noqa: F401  (zstd/blosc filters for verification)
except ImportError:
    pass


def leaves(obj, prefix, out):
    """Flatten an AnnData into {path: value} at the granularity of the post's table."""
    if isinstance(obj, ad.AnnData):
        p = (prefix + "/") if prefix else ""
        if obj.X is not None:
            out[p + "X"] = obj.X
        frame_leaves(obj.obs, p + "obs", out)
        frame_leaves(obj.var, p + "var", out)
        for k in ("obsm", "varm", "obsp", "varp", "layers"):
            for kk, v in dict(getattr(obj, k)).items():
                if kk is None:
                    continue
                leaves(v, f"{p}{k}/{kk}", out)
        for kk, v in obj.uns.items():
            leaves(v, f"{p}uns/{kk}", out)
        if obj.raw is not None:
            out[p + "raw/X"] = obj.raw.X
            frame_leaves(obj.raw.var, p + "raw/var", out)
            for kk, v in dict(obj.raw.varm).items():
                leaves(v, f"{p}raw/varm/{kk}", out)
        return
    if isinstance(obj, pd.DataFrame):
        frame_leaves(obj, prefix, out)
        return
    if isinstance(obj, dict) or (hasattr(obj, "keys") and not isinstance(obj, (pd.Series, np.ndarray))):
        if len(obj) == 0:
            out[prefix] = {}
        for k, v in obj.items():
            leaves(v, f"{prefix}/{k}", out)
        return
    out[prefix] = obj


def frame_leaves(df, prefix, out):
    out[prefix + "/<index>"] = df.index
    out[prefix + "/<columns>"] = list(df.columns)
    for c in df.columns:
        out[f"{prefix}/{c}"] = df[c]


def as_np(v):
    if isinstance(v, (pd.Series, pd.Index)):
        return v.array
    return v


def nan_equal(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind in "OUST" or b.dtype.kind in "OUST":
        fa = [x.decode() if isinstance(x, bytes) else x for x in a.reshape(-1)]
        fb = [x.decode() if isinstance(x, bytes) else x for x in b.reshape(-1)]
        return fa == fb
    if a.dtype.kind == "c" or b.dtype.kind == "c":
        return nan_equal(a.real, b.real) and nan_equal(a.imag, b.imag)
    if a.dtype.kind == "b" or b.dtype.kind == "b":
        return np.array_equal(a.astype(bool), b.astype(bool)) and a.dtype.kind == b.dtype.kind
    af = a.astype(np.float64) if a.dtype.kind == "f" else a
    bf = b.astype(np.float64) if b.dtype.kind == "f" else b
    if af.dtype.kind == "f" or bf.dtype.kind == "f":
        af = af.astype(np.float64)
        bf = bf.astype(np.float64)
        return bool(np.all((af == bf) | (np.isnan(af) & np.isnan(bf))))
    return np.array_equal(af.astype(object) if af.dtype.kind == "u" else af,
                          bf.astype(object) if bf.dtype.kind == "u" else bf)


def values_equal(a, b):
    a = as_np(a)
    b = as_np(b)
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, dict) or isinstance(b, dict):
        return isinstance(a, dict) and isinstance(b, dict) and len(a) == len(b) == 0
    if isinstance(a, list) or isinstance(b, list):
        return list(a) == list(b)
    if sp.issparse(a) or sp.issparse(b):
        if not (sp.issparse(a) and sp.issparse(b)) or a.shape != b.shape:
            return False
        ca = sp.csr_matrix(a)
        cb = sp.csr_matrix(b)
        ca.sum_duplicates()
        cb.sum_duplicates()
        return (np.array_equal(ca.indptr, cb.indptr) and np.array_equal(ca.indices, cb.indices)
                and nan_equal(ca.data, cb.data))
    if isinstance(a, pd.Categorical) or isinstance(getattr(a, "dtype", None), pd.CategoricalDtype):
        if not (isinstance(b, pd.Categorical) or isinstance(getattr(b, "dtype", None), pd.CategoricalDtype)):
            return False
        a = pd.Categorical(a)
        b = pd.Categorical(b)
        return (a.ordered == b.ordered and nan_equal(np.asarray(a.categories), np.asarray(b.categories))
                and np.array_equal(a.codes, b.codes))
    if isinstance(a, pd.api.extensions.ExtensionArray) or isinstance(b, pd.api.extensions.ExtensionArray):
        na_a = np.asarray(pd.isna(a))
        na_b = np.asarray(pd.isna(b))
        if na_a.shape != na_b.shape or not np.array_equal(na_a, na_b):
            return False
        va = np.asarray(pd.array(a)[~na_a].to_numpy() if isinstance(a, pd.api.extensions.ExtensionArray) else np.asarray(a)[~na_a])
        vb = np.asarray(pd.array(b)[~na_b].to_numpy() if isinstance(b, pd.api.extensions.ExtensionArray) else np.asarray(b)[~na_b])
        return nan_equal(va, vb)
    if isinstance(a, ad.AnnData) or isinstance(b, ad.AnnData):
        return False  # nested AnnData are flattened by leaves()
    if isinstance(a, np.ndarray) and a.dtype.names is not None:
        if not (isinstance(b, np.ndarray) and b.dtype.names == a.dtype.names) or a.shape != b.shape:
            return False
        return all(nan_equal(a[n], b[n]) for n in a.dtype.names)
    if np.isscalar(a) or np.isscalar(b) or isinstance(a, (str, bytes)) or isinstance(b, (str, bytes)):
        if type(a) is not type(b) and not (isinstance(a, (np.generic, int, float, bool, str)) and isinstance(b, (np.generic, int, float, bool, str))):
            return False
        if isinstance(a, str) or isinstance(b, str):
            return isinstance(a, str) and isinstance(b, str) and a == b
        return nan_equal(np.asarray(a), np.asarray(b)) and np.asarray(a).shape == np.asarray(b).shape
    if type(a).__module__.startswith("awkward"):
        import awkward as ak
        return ak.to_list(a) == ak.to_list(b)
    try:
        return nan_equal(np.asarray(a), np.asarray(b))
    except Exception:
        return False


def h5_meta(f: h5py.File, path: str):
    """Encoding + dtypes of an element on disk (including its components)."""
    if path.endswith("/<index>"):
        g = f[path[: -len("/<index>")]]
        path = g.name.lstrip("/") + "/" + g.attrs.get("_index", "_index")
    elif path.endswith("/<columns>"):
        return ("column-order",)
    if path not in f:
        return None
    o = f[path]

    def dt(d):
        s = h5py.check_string_dtype(d.dtype)
        if s is not None:
            return f"str:{s.encoding}:{'vlen' if s.length is None else s.length}"
        if d.dtype.names:
            return "compound(" + ",".join(f"{n}:{dt_field(d.dtype[n])}" for n in d.dtype.names) + ")"
        return str(d.dtype)

    def dt_field(t):
        s = h5py.check_string_dtype(t)
        if s is not None:
            return f"str:{s.encoding}"
        return str(t)

    enc = o.attrs.get("encoding-type", "")
    ver = o.attrs.get("encoding-version", "")
    if isinstance(o, h5py.Dataset):
        return (enc, ver, dt(o), o.shape)
    parts = []
    for k in sorted(o.keys()):
        c = o[k]
        if isinstance(c, h5py.Dataset) and enc in ("csr_matrix", "csc_matrix", "categorical",
                                                    "nullable-integer", "nullable-boolean",
                                                    "nullable-string-array", "awkward-array"):
            parts.append((k, c.attrs.get("encoding-type", ""), dt(c)))
    extra = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in o.attrs.items()
             if k in ("ordered", "na-value", "shape")}
    return (enc, ver, tuple(parts), json.dumps(extra, sort_keys=True, default=str))


def main():
    orig, rfile = sys.argv[1], sys.argv[2]
    out_json = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        A = ad.read_h5ad(orig)
        tmpdir = tempfile.mkdtemp(prefix="annr_ref_")
        ref = Path(tmpdir) / "py_ref.h5ad"
        # Classic (anndata < 0.11 compatible) encoding for strings without
        # missing values, independent of pandas' string inference.
        try:
            try:
                ad.settings.allow_write_nullable_strings = False
                A.write_h5ad(ref)
            except ValueError:  # strings with missing values need the nullable encoding
                ad.settings.allow_write_nullable_strings = True
                A.write_h5ad(ref)
        except Exception as e:  # noqa: BLE001  (e.g. anndata cannot re-write complex scalars)
            ref_error = repr(e)
            with h5py.File(ref, "w"):
                pass
        else:
            ref_error = None
        res = {"file": orig, "rfile": rfile, "counts": {}, "details": []}
        try:
            B = ad.read_h5ad(rfile)
        except Exception as e:  # noqa: BLE001
            la = {}
            leaves(A, "", la)
            res["counts"] = {"unreadable": len(la)}
            res["details"] = [f"python cannot read R file: {e!r}"]
            print(json.dumps(res["counts"]))
            if out_json:
                Path(out_json).write_text(json.dumps(res, indent=1))
            return
    la, lb = {}, {}
    leaves(A, "", la)
    leaves(B, "", lb)
    counts = {k: 0 for k in ("equal_same", "equal_diff", "values_differ", "unreadable", "missing")}
    with h5py.File(ref, "r") as fr, h5py.File(rfile, "r") as fw, h5py.File(orig, "r") as fo:
        modern = "encoding-type" in fo.attrs
        for p, va in la.items():
            if p not in lb:
                counts["missing"] += 1
                res["details"].append(f"missing: {p}")
                continue
            vb = lb[p]
            try:
                eq = values_equal(va, vb)
            except Exception as e:  # noqa: BLE001
                eq = False
                res["details"].append(f"compare error {p}: {e!r}")
            if not eq:
                counts["values_differ"] += 1
                res["details"].append(f"values_differ: {p}: {type(as_np(va)).__name__} vs {type(as_np(vb)).__name__}")
                continue
            # Reference metadata: the original file itself when it is in the
            # modern format, otherwise what Python anndata writes.
            mo = h5_meta(fo, p) if modern else None
            ma = mo if mo is not None else h5_meta(fr, p)
            mb = h5_meta(fw, p)
            if ma == mb:
                counts["equal_same"] += 1
            else:
                counts["equal_diff"] += 1
                res["details"].append(f"equal_diff: {p}: python={ma} R={mb}")
    if ref_error:
        res["details"].append(f"note: python could not write a reference copy ({ref_error[:120]})")
    extra = sorted(set(lb) - set(la))
    if extra:
        res["details"].append(f"extra elements written by R: {extra[:10]}")
    res["counts"] = counts
    res["n_leaves"] = len(la)
    print(json.dumps(counts))
    for d in res["details"][:40]:
        print("   ", d)
    if out_json:
        Path(out_json).write_text(json.dumps(res, indent=1, default=str))
    if "--keep-ref" not in sys.argv:
        ref.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
