"""Generate synthetic .h5ad files covering every anndata encoding and dtype
edge case, plus hand-built HDF5 layouts anndata itself never writes.

Usage: python make_synthetic.py <out_dir>
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad

try:
    import awkward as ak
except ImportError:  # pragma: no cover
    ak = None
try:
    import hdf5plugin
except ImportError:  # pragma: no cover
    hdf5plugin = None

rng = np.random.default_rng(20261026)
UNICODE = ["α-β", "细胞_1", "naïve T", "🧬cell", "", " spaced ", "tab\tsep", "new\nline", "ünïcödé", "Ж"]


def obs_frame(n, prefix):
    names = [f"{prefix}{i}_{UNICODE[i % len(UNICODE)].strip() or 'x'}" for i in range(n)]
    df = pd.DataFrame(index=pd.Index(names, dtype=object))
    df["f64"] = rng.normal(size=n)
    df.loc[df.index[::7], "f64"] = np.nan
    df["f32"] = rng.normal(size=n).astype(np.float32)
    df["f32_inf"] = np.where(np.arange(n) % 5 == 0, np.inf, -np.inf).astype(np.float32)
    df["i8"] = rng.integers(-128, 127, size=n).astype(np.int8)
    df["i16"] = rng.integers(-32768, 32767, size=n).astype(np.int16)
    df["i32"] = rng.integers(-(2**31) + 1, 2**31 - 1, size=n).astype(np.int32)
    df["i64_small"] = rng.integers(-1000, 1000, size=n).astype(np.int64)
    df["i64_big"] = (rng.integers(0, 1000, size=n) + 2**40).astype(np.int64)
    df["u8"] = rng.integers(0, 255, size=n).astype(np.uint8)
    df["u16"] = rng.integers(0, 65535, size=n).astype(np.uint16)
    df["u32_big"] = (rng.integers(0, 1000, size=n) + 2**31 + 5).astype(np.uint32)
    df["u64"] = rng.integers(0, 2**20, size=n).astype(np.uint64)
    df["boolean"] = rng.random(n) > 0.5
    df["cat_str"] = pd.Categorical(rng.choice(["a", "b", "c", "ü"], size=n))
    df["cat_ordered"] = pd.Categorical(rng.choice(["low", "mid", "high"], size=n),
                                       categories=["low", "mid", "high"], ordered=True)
    cm = pd.Categorical(rng.choice(["x", "y"], size=n))
    cm[::3] = np.nan
    df["cat_missing"] = cm
    df["cat_int"] = pd.Categorical(rng.choice([1, 5, 10], size=n))
    df["cat_float"] = pd.Categorical(rng.choice([0.5, 1.25, 3.0], size=n))
    df["cat_bool"] = pd.Categorical(rng.choice([True, False], size=n))
    df["cat_many"] = pd.Categorical([f"c{i % 300}" for i in range(n)])
    df["strings"] = pd.array([UNICODE[i % len(UNICODE)] for i in range(n)], dtype=object)
    ni = pd.array(rng.integers(0, 100, size=n), dtype="Int64")
    ni[::4] = pd.NA
    df["nullable_int"] = ni
    df["nullable_int_full"] = pd.array(rng.integers(0, 100, size=n), dtype="Int32")
    nb = pd.array(rng.random(n) > 0.5, dtype="boolean")
    nb[1::5] = pd.NA
    df["nullable_bool"] = nb
    ns = pd.array([UNICODE[i % len(UNICODE)] for i in range(n)], dtype=pd.StringDtype("python", na_value=pd.NA))
    ns[2::6] = pd.NA
    df["nullable_str"] = ns
    return df


def make_main(path: Path, compression=None, n_obs=257, n_var=131):
    ad.settings.allow_write_nullable_strings = True
    obs = obs_frame(n_obs, "c")
    var = obs_frame(n_var, "g")
    var.index.name = "gene_ids"
    X = sp.random(n_obs, n_var, density=0.1, format="csr", dtype=np.float32, random_state=1)
    X.data[::17] = np.nan
    X.data[::23] = np.inf
    X.data[1::29] = 0.0  # explicit zeros are kept
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.layers["dense_f64"] = rng.normal(size=(n_obs, n_var))
    a.layers["dense_f64"][0, :5] = [np.nan, np.inf, -np.inf, -0.0, 1e-308]
    a.layers["dense_f32"] = rng.normal(size=(n_obs, n_var)).astype(np.float32)
    a.layers["dense_f16"] = rng.normal(size=(n_obs, n_var)).astype(np.float16)
    a.layers["dense_i8"] = rng.integers(-128, 127, size=(n_obs, n_var)).astype(np.int8)
    a.layers["dense_u16"] = rng.integers(0, 65535, size=(n_obs, n_var)).astype(np.uint16)
    a.layers["dense_i64_big"] = rng.integers(2**40, 2**41, size=(n_obs, n_var)).astype(np.int64)
    a.layers["dense_u64"] = rng.integers(0, 2**62, size=(n_obs, n_var), dtype=np.uint64) * np.uint64(3)
    a.layers["dense_bool"] = rng.random((n_obs, n_var)) > 0.7
    a.layers["csr_i64"] = sp.random(n_obs, n_var, density=0.05, format="csr", random_state=2,
                                    data_rvs=lambda k: rng.integers(1, 2**40, size=k)).astype(np.int64)
    a.layers["csc_f32"] = sp.random(n_obs, n_var, density=0.05, format="csc", dtype=np.float32, random_state=3)
    a.layers["csr_bool"] = sp.random(n_obs, n_var, density=0.05, format="csr", random_state=4).astype(bool)
    a.layers["csr_empty"] = sp.csr_matrix((n_obs, n_var), dtype=np.float64)
    big = sp.random(n_obs, n_var, density=0.05, format="csr", random_state=5)
    big.indptr = big.indptr.astype(np.int64)
    big.indices = big.indices.astype(np.int64)
    a.layers["csr_int64_index"] = big
    a.obsm["dense_f32"] = rng.normal(size=(n_obs, 7)).astype(np.float32)
    a.obsm["dense_i32"] = rng.integers(-5, 5, size=(n_obs, 3)).astype(np.int32)
    a.obsm["one_col"] = rng.normal(size=(n_obs, 1))
    a.obsm["zero_col"] = np.zeros((n_obs, 0))
    a.obsm["sparse"] = sp.random(n_obs, 40, density=0.1, format="csr", random_state=6)
    odf = pd.DataFrame({"num": rng.normal(size=n_obs),
                        "cat": pd.Categorical(rng.choice(["p", "q"], size=n_obs))}, index=obs.index)
    a.obsm["frame"] = odf
    a.varm["loadings"] = rng.normal(size=(n_var, 11))
    a.obsp["graph"] = sp.random(n_obs, n_obs, density=0.02, format="csr", dtype=np.float32, random_state=7)
    a.obsp["dense"] = rng.normal(size=(n_obs, n_obs))
    a.varp["csc"] = sp.random(n_var, n_var, density=0.05, format="csc", random_state=8)
    u = a.uns
    u["s_int"] = 42
    u["s_float"] = 3.25
    u["s_bool"] = True
    u["s_str"] = "héllo 🌍"
    u["s_empty_str"] = ""
    u["s_np_f32"] = np.float32(1.5)
    u["s_np_i8"] = np.int8(-7)
    u["s_np_u64"] = np.uint64(2**63 + 11)
    u["s_complex"] = np.complex128(1.5 - 2.0j)
    u["a_f64"] = rng.normal(size=9)
    u["a_one"] = np.array([5])
    u["a_one_str"] = np.array(["#1f77b4"], dtype=object)
    u["a_empty"] = np.array([], dtype=np.float64)
    u["a_empty_str"] = np.array([], dtype=object)
    u["a_2d_int"] = rng.integers(0, 9, size=(4, 6))
    u["a_3d"] = rng.normal(size=(2, 3, 4))
    u["a_3d_i16"] = rng.integers(-300, 300, size=(3, 2, 5)).astype(np.int16)
    u["a_str"] = np.array(UNICODE, dtype=object)
    u["a_str_2d"] = np.array([["a", "bb"], ["ccc", "🧪"], ["", "x"]], dtype=object)
    u["a_bool"] = np.array([True, False, True])
    u["a_complex"] = rng.normal(size=4) + 1j * rng.normal(size=4)
    u["a_u32_big"] = np.array([2**32 - 1, 0, 2**31], dtype=np.uint32)
    u["nested"] = {"level1": {"level2": {"leaf": np.arange(3), "word": "deep"}}, "x": 1.0}
    u["empty_dict"] = {}
    u["frame_rangeindex"] = pd.DataFrame({"a": np.arange(5), "b": list("vwxyz")})
    u["frame_strindex"] = pd.DataFrame({"v": rng.normal(size=3)}, index=["r1", "r2", "r3"])
    u["recarray"] = np.rec.fromarrays(
        [np.array(["g1", "g2", "g3"], dtype=object), np.array([1.5, 2.5, 3.5], dtype=np.float32),
         np.array([1, 2, 3], dtype=np.int64)], names=["name", "score", "rank"])
    u["categorical"] = pd.Categorical(["u", "v", "u", np.nan])
    u["sparse"] = sp.random(5, 6, density=0.3, format="csr", random_state=9)
    u["none"] = None
    u["nullable_int"] = pd.array([1, None, 3], dtype="Int64")
    u["nullable_str"] = pd.array(["a", None, "c"], dtype=pd.StringDtype("python", na_value=pd.NA))
    u["colors"] = ["#aabbcc", "#ddeeff"]
    if ak is not None:
        u["awkward"] = ak.Array([[1, 2, 3], [], [4, 5]])
    sub = ad.AnnData(X=np.arange(6, dtype=np.float32).reshape(2, 3))
    sub.obs["k"] = ["p", "q"]
    u["nested_anndata"] = sub
    raw_var = obs_frame(n_var + 20, "rg")[["f64", "cat_str"]]
    a.raw = ad.AnnData(X=sp.random(n_obs, n_var + 20, density=0.1, format="csr", dtype=np.float32,
                                   random_state=10), var=raw_var)
    a.raw.varm["raw_loadings"] = rng.normal(size=(n_var + 20, 3))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a.write_h5ad(path, compression=compression,
                     compression_opts=4 if compression == "gzip" else None)
    return a


def make_min_dtype_indices(path: Path):
    with ad.settings.override(write_csr_csc_indices_with_min_possible_dtype=True):
        a = ad.AnnData(X=sp.random(300, 200, density=0.1, format="csr", dtype=np.float32, random_state=11))
        a.obsm["wide"] = sp.random(300, 70000, density=0.001, format="csr", random_state=12)
        a.write_h5ad(path)


def make_unsorted(path: Path):
    """CSR with unsorted and duplicate indices written verbatim."""
    with h5py.File(path, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        g = f.create_group("X")
        g.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": (3, 4)})
        g["data"] = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        g["indices"] = np.array([3, 1, 1, 0, 2, 2], dtype=np.int32)  # row0: 3,1 ; row1: 1,0,2 ; row2: 2 (dup across? no)
        g["indptr"] = np.array([0, 2, 5, 6], dtype=np.int32)
        dup = f.create_group("layers").create_group("dups")
        f["layers"].attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
        dup.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": (3, 4)})
        dup["data"] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        dup["indices"] = np.array([2, 2, 0, 3, 1], dtype=np.int32)  # (0,2) duplicated, row 1 unsorted
        dup["indptr"] = np.array([0, 2, 4, 5], dtype=np.int32)
        for k, n in (("obs", 3), ("var", 4)):
            d = f.create_group(k)
            d.attrs.update({"encoding-type": "dataframe", "encoding-version": "0.2.0", "_index": "_index",
                            "column-order": np.array([], dtype=float)})
            ds = d.create_dataset("_index", data=np.array([f"{k}{i}" for i in range(n)], dtype=object),
                                  dtype=h5py.string_dtype())
            ds.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})


def make_layouts(path: Path):
    """Valid h5ad files using HDF5 features anndata does not produce."""
    n_obs, n_var = 503, 307
    dense = rng.normal(size=(n_obs, n_var))
    counts = sp.random(n_obs, n_var, density=0.08, format="csr", dtype=np.float32, random_state=13)
    a = ad.AnnData(X=counts.copy())
    a.layers["dense_chunked"] = dense
    a.layers["dense_be"] = dense.copy()
    a.layers["dense_fletcher"] = dense.astype(np.float32)
    a.layers["dense_scaleoffset"] = rng.integers(0, 1000, size=(n_obs, n_var)).astype(np.int32)
    a.layers["dense_unalloc"] = np.zeros((n_obs, n_var), dtype=np.float32)
    a.layers["dense_compact"] = np.zeros((n_obs, n_var))  # replaced below (LZF, column chunks)
    a.obs["be_int"] = rng.integers(0, 10**12, size=n_obs).astype(np.int64)
    a.obs["fixed_str"] = [f"s{i}" for i in range(n_obs)]
    a.obsm["gzip_2d_chunks"] = rng.normal(size=(n_obs, 37)).astype(np.float32)
    a.obsm["shuffle_lzf"] = rng.integers(0, 2**20, size=(n_obs, 13)).astype(np.int64)
    a.uns["small_compact"] = np.arange(10, dtype=np.int16)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a.write_h5ad(path)
    with h5py.File(path, "r+") as f:
        def replace(name, data, **kw):
            attrs = dict(f[name].attrs)
            del f[name]
            d = f.create_dataset(name, data=data, **kw)
            d.attrs.update(attrs)
        replace("layers/dense_chunked", dense, chunks=(64, 50), compression="gzip", shuffle=True)
        replace("layers/dense_be", dense.astype(">f8"), chunks=(100, 100), compression="gzip")
        replace("layers/dense_fletcher", dense.astype(np.float32), chunks=(128, 64), fletcher32=True)
        replace("layers/dense_scaleoffset", f["layers/dense_scaleoffset"][()], chunks=(128, 128), scaleoffset=0)
        # chunked but only partially written: unwritten chunks hold the fill value
        attrs = dict(f["layers/dense_unalloc"].attrs)
        del f["layers/dense_unalloc"]
        d = f.create_dataset("layers/dense_unalloc", shape=(n_obs, n_var), dtype=np.float32, chunks=(100, 100),
                             fillvalue=np.float32(2.5))
        d[:100, :100] = 1.0
        d.attrs.update(attrs)
        # column-oriented chunks with LZF (chunk edges on both axes)
        replace("layers/dense_compact", rng.normal(size=(n_obs, n_var)), chunks=(n_obs, 16), compression="lzf")
        replace("obs/be_int", f["obs/be_int"][()].astype(">i8"))
        del f["obs/fixed_str"]
        d = f.create_dataset("obs/fixed_str", data=np.array([f"s{i}".encode() for i in range(n_obs)], dtype="S8"))
        d.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
        replace("obsm/gzip_2d_chunks", f["obsm/gzip_2d_chunks"][()], chunks=(97, 11), compression="gzip",
                compression_opts=9)
        replace("obsm/shuffle_lzf", f["obsm/shuffle_lzf"][()], chunks=(64, 5), compression="lzf", shuffle=True)
        attrs = dict(f["uns/small_compact"].attrs)
        data = f["uns/small_compact"][()]
        del f["uns/small_compact"]
        dcpl = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
        dcpl.set_layout(h5py.h5d.COMPACT)
        sid = h5py.h5s.create_simple(data.shape)
        dsid = h5py.h5d.create(f.id, b"uns/small_compact", h5py.h5t.STD_I16LE, sid, dcpl=dcpl)
        dsid.write(h5py.h5s.ALL, h5py.h5s.ALL, data)
        f["uns/small_compact"].attrs.update(attrs)


def make_userblock(path: Path):
    a = ad.AnnData(X=sp.random(50, 40, density=0.2, format="csr", dtype=np.float32, random_state=14))
    a.layers["d"] = rng.normal(size=(50, 40))
    tmp = path.with_suffix(".tmp.h5ad")
    a.write_h5ad(tmp)
    with h5py.File(tmp, "r") as src, h5py.File(path, "w", userblock_size=512) as dst:
        for k, v in src.attrs.items():
            dst.attrs[k] = v
        for k in src:
            src.copy(src[k], dst, name=k)
    tmp.unlink()


def make_zstd(path: Path):
    if hdf5plugin is None:
        return
    a = ad.AnnData(X=sp.random(400, 300, density=0.1, format="csr", dtype=np.float32, random_state=15))
    a.layers["dense"] = rng.normal(size=(400, 300))
    a.write_h5ad(path, compression=hdf5plugin.Zstd(clevel=5))


def make_dense_x(path: Path):
    """Dense X with gzip 2-D chunking (row/column edges not multiples of chunks)."""
    a = ad.AnnData(X=rng.normal(size=(1001, 777)).astype(np.float32))
    a.obs["g"] = pd.Categorical(rng.choice(list("abc"), size=1001))
    a.write_h5ad(path, compression="gzip")


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    make_main(out / "syn_all_types.h5ad")
    make_main(out / "syn_all_types_gzip.h5ad", compression="gzip")
    make_main(out / "syn_all_types_lzf.h5ad", compression="lzf")
    make_min_dtype_indices(out / "syn_min_dtype_indices.h5ad")
    make_unsorted(out / "syn_unsorted_dups.h5ad")
    make_layouts(out / "syn_hdf5_layouts.h5ad")
    make_userblock(out / "syn_userblock.h5ad")
    make_zstd(out / "syn_zstd.h5ad")
    make_dense_x(out / "syn_dense_x_gzip.h5ad")
    for p in sorted(out.glob("syn_*.h5ad")):
        print(p)


if __name__ == "__main__":
    main()
