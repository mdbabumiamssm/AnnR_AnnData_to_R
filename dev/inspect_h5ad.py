"""Inventory encodings, dtypes, filters and chunking across h5ad files."""
import sys, collections, h5py, numpy as np
FILT = {1: "gzip", 2: "shuffle", 3: "fletcher32", 4: "szip", 32000: "lzf", 32001: "blosc", 32015: "zstd"}
enc_counter = collections.Counter(); dtype_counter = collections.Counter(); filt_counter = collections.Counter()
def filters(ds):
    p = ds.id.get_create_plist(); out = []
    for i in range(p.get_nfilters()):
        out.append(FILT.get(p.get_filter(i)[0], str(p.get_filter(i)[0])))
    return "+".join(out) or "none"
for path in sys.argv[1:]:
    try:
        f = h5py.File(path, "r")
    except Exception as e:
        print("### OPEN FAIL", path, e); continue
    print(f"\n### {path}  root-attrs={dict(f.attrs)}")
    def visit(name, obj):
        enc = obj.attrs.get("encoding-type", None)
        enc = enc.decode() if isinstance(enc, bytes) else enc
        enc_counter[enc] += 1
        if isinstance(obj, h5py.Dataset):
            dt = str(obj.dtype) if obj.dtype.kind != "O" else ("vlen-str" if h5py.check_string_dtype(obj.dtype) else "O")
            if obj.dtype.names: dt = "compound" + str(obj.dtype.names)[:60]
            dtype_counter[dt] += 1
            fl = filters(obj); filt_counter[fl] += 1
            depth = name.count("/")
            if depth <= 2 or (obj.size or 0) > 1e6 or obj.shape is None:
                print(f"  D {name:55s} enc={enc} shape={obj.shape} dtype={dt} chunks={obj.chunks} filt={fl}")
        else:
            if name.count("/") <= 1:
                print(f"  G {name:55s} enc={enc} attrs={ {k: (v.tolist() if hasattr(v,'tolist') else v) for k,v in obj.attrs.items() if k!='column-order'} }")
    f.visititems(visit)
print("\nENCODINGS", enc_counter.most_common()); print("DTYPES", dtype_counter.most_common()); print("FILTERS", filt_counter.most_common())
