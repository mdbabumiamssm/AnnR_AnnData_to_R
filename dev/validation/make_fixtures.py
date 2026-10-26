"""Create the small test fixtures shipped in tests/testthat/fixtures.

Each fixture comes with expected values computed by Python anndata, stored in
expected.json, so the R tests check AnnR against the reference implementation.

Usage: python make_fixtures.py <fixtures_dir>
"""
import json
import sys
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad

rng = np.random.default_rng(7)


def modern(path):
    n_obs, n_var = 23, 17
    obs = pd.DataFrame(index=[f"cell-{i}-ü" for i in range(n_obs)])
    obs["group"] = pd.Categorical(rng.choice(["b", "a", "c"], size=n_obs), categories=["b", "a", "c"])
    obs["level"] = pd.Categorical(rng.choice(["lo", "hi"], size=n_obs), categories=["lo", "hi"], ordered=True)
    obs["score"] = rng.normal(size=n_obs).astype(np.float32)
    obs["count"] = rng.integers(0, 1000, size=n_obs).astype(np.int64)
    obs["flag"] = rng.random(n_obs) > 0.5
    ni = pd.array(rng.integers(0, 5, size=n_obs), dtype="Int64")
    ni[::4] = pd.NA
    obs["maybe"] = ni
    obs["label"] = np.array([f"L{i % 3}" for i in range(n_obs)], dtype=object)
    var = pd.DataFrame(index=[f"G{i}" for i in range(n_var)])
    var["hv"] = rng.random(n_var) > 0.5
    X = sp.random(n_obs, n_var, density=0.3, format="csr", dtype=np.float32, random_state=3)
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.layers["counts"] = sp.random(n_obs, n_var, density=0.3, format="csc", random_state=4,
                                   data_rvs=lambda k: rng.integers(1, 50, size=k)).astype(np.int32)
    a.layers["dense"] = rng.normal(size=(n_obs, n_var))
    a.obsm["X_pca"] = rng.normal(size=(n_obs, 5)).astype(np.float32)
    a.varm["PCs"] = rng.normal(size=(n_var, 5))
    a.obsp["distances"] = sp.random(n_obs, n_obs, density=0.1, format="csr", random_state=5)
    a.uns["title"] = "fixture ✓"
    a.uns["n"] = 3
    a.uns["colors"] = np.array(["#000000"], dtype=object)
    a.uns["params"] = {"k": 15, "metric": "cosine", "weights": np.array([0.5, 1.5])}
    a.uns["big"] = np.array([2**62 + 1, -(2**62) - 3], dtype=np.int64)
    a.raw = ad.AnnData(X=sp.random(n_obs, n_var + 3, density=0.3, format="csr", dtype=np.float32,
                                   random_state=6), var=pd.DataFrame(index=[f"R{i}" for i in range(n_var + 3)]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a.write_h5ad(path, compression="gzip")
    b = ad.read_h5ad(path)
    Xd = b.X.toarray()
    return {
        "shape": list(b.shape),
        "obs_names": list(b.obs_names),
        "var_names": list(b.var_names),
        "X_dense": Xd.astype(float).tolist(),
        "counts_dense": b.layers["counts"].toarray().astype(int).tolist(),
        "dense": b.layers["dense"].tolist(),
        "group": [str(x) for x in b.obs["group"]],
        "group_levels": list(b.obs["group"].cat.categories),
        "level_levels": list(b.obs["level"].cat.categories),
        "score": b.obs["score"].astype(float).tolist(),
        "count": b.obs["count"].astype(int).tolist(),
        "flag": b.obs["flag"].tolist(),
        "maybe": [None if pd.isna(x) else int(x) for x in b.obs["maybe"]],
        "X_pca": b.obsm["X_pca"].astype(float).tolist(),
        "PCs": b.varm["PCs"].tolist(),
        "distances_dense": b.obsp["distances"].toarray().tolist(),
        "raw_shape": list(b.raw.X.shape),
        "raw_X_sum": float(b.raw.X.sum()),
        "big": [str(v) for v in b.uns["big"]],
    }


def legacy(path):
    """A pre-0.7 style file: compound obs/var, h5sparse X, raw.X, uns categories."""
    n_obs, n_var = 6, 4
    with h5py.File(path, "w") as f:
        obs = np.zeros(n_obs, dtype=[("index", "S6"), ("louvain", "i1"), ("n_genes", "<i8")])
        obs["index"] = [f"c{i}".encode() for i in range(n_obs)]
        obs["louvain"] = [0, 1, 1, 0, 2, 2]
        obs["n_genes"] = np.arange(n_obs) * 10
        f["obs"] = obs
        var = np.zeros(n_var, dtype=[("index", "S6"), ("means", "<f4")])
        var["index"] = [f"g{i}".encode() for i in range(n_var)]
        var["means"] = [0.5, 1.5, 2.5, 3.5]
        f["var"] = var
        X = sp.random(n_obs, n_var, density=0.5, format="csr", random_state=1)
        g = f.create_group("X")
        g.attrs["h5sparse_format"] = "csr"
        g.attrs["h5sparse_shape"] = X.shape
        g["data"], g["indices"], g["indptr"] = X.data, X.indices, X.indptr
        rX = sp.random(n_obs, n_var + 1, density=0.5, format="csr", random_state=2)
        r = f.create_group("raw.X")
        r.attrs["h5sparse_format"] = "csr"
        r.attrs["h5sparse_shape"] = rX.shape
        r["data"], r["indices"], r["indptr"] = rX.data, rX.indices, rX.indptr
        rv = np.zeros(n_var + 1, dtype=[("index", "S6")])
        rv["index"] = [f"r{i}".encode() for i in range(n_var + 1)]
        f["raw.var"] = rv
        u = f.create_group("uns")
        u["louvain_categories"] = np.array([b"0", b"1", b"2"])
        u["scalar"] = 2.5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = ad.read_h5ad(path)
    return {
        "obs_names": list(b.obs_names),
        "louvain": [str(x) for x in b.obs["louvain"]],
        "louvain_levels": [str(x) for x in b.obs["louvain"].cat.categories],
        "n_genes": b.obs["n_genes"].astype(int).tolist(),
        "means": b.var["means"].astype(float).tolist(),
        "X_dense": b.X.toarray().tolist(),
        "raw_shape": list(b.raw.X.shape),
        "uns_keys": sorted(b.uns.keys()),
    }


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    expected = {"modern": modern(out / "modern.h5ad"), "legacy": legacy(out / "legacy.h5ad")}
    (out / "expected.json").write_text(json.dumps(expected, indent=0, ensure_ascii=False))
    for p in sorted(out.iterdir()):
        print(p, p.stat().st_size)


if __name__ == "__main__":
    main()
