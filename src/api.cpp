// R-callable entry points.
#include "common.h"

namespace annr {
SEXP read_h5ad_impl(const std::string& path, SEXP opts);
SEXP read_elem_impl(const std::string& path, const std::string& elem, SEXP opts);
SEXP w_open(const std::string& path, int threads);
void w_close(SEXP xp);
void w_group(SEXP xp, const std::string& path, const std::string& enc, const std::string& ver);
void w_attr(SEXP xp, const std::string& path, const std::string& name, SEXP value, const std::string& kind);
std::string w_array(SEXP xp, const std::string& path, SEXP x, SEXP shape, bool colmajor,
                    const std::string& hint, const std::string& enc, const std::string& ver,
                    SEXP comp_opts, bool unlimited);
void w_null(SEXP xp, const std::string& path);
SEXP w_sparse(SEXP xp, const std::string& path, const std::string& fmt, SEXP shape, SEXP p, SEXP idx,
              SEXP x, SEXP hints, SEXP comp_opts);
void w_recarray(SEXP xp, const std::string& path, SEXP cols, SEXP names, SEXP hints);
}

using namespace annr;

[[cpp11::init]]
void annr_init(DllInfo* dll) {
  register_altrep_classes(dll);
  H5open();
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  register_filters();
}

[[cpp11::register]]
SEXP cpp_read_h5ad(std::string path, SEXP opts) {
  return read_h5ad_impl(path, opts);
}

[[cpp11::register]]
SEXP cpp_read_elem(std::string path, std::string elem, SEXP opts) {
  return read_elem_impl(path, elem, opts);
}

[[cpp11::register]]
double cpp_external_bytes() { return external_live_bytes(); }

[[cpp11::register]]
cpp11::writable::strings cpp_build_info() {
  unsigned maj = 0, min = 0, rel = 0;
  H5get_libversion(&maj, &min, &rel);
  std::string h5 = std::to_string(maj) + "." + std::to_string(min) + "." + std::to_string(rel);
  cpp11::writable::strings out({h5,
#ifdef ANNR_HAVE_LIBDEFLATE
                                "libdeflate",
#else
                                "zlib",
#endif
#ifdef ANNR_HAVE_ZSTD
                                "yes"
#else
                                "no"
#endif
  });
  out.names() = {"hdf5", "gzip_codec", "zstd"};
  return out;
}

namespace {

struct LsRow {
  std::string path, kind, encoding, dtype, shape, chunks, filters;
  double stored = 0;
};

std::string filter_name(H5Z_filter_t id) {
  switch (id) {
    case H5Z_FILTER_DEFLATE: return "gzip";
    case H5Z_FILTER_SHUFFLE: return "shuffle";
    case H5Z_FILTER_FLETCHER32: return "fletcher32";
    case H5Z_FILTER_SZIP: return "szip";
    case H5Z_FILTER_NBIT: return "nbit";
    case H5Z_FILTER_SCALEOFFSET: return "scaleoffset";
    case FILTER_LZF: return "lzf";
    case FILTER_ZSTD: return "zstd";
    case 32001: return "blosc";
    default: return std::to_string(id);
  }
}

std::string dims_str(const std::vector<hsize_t>& d) {
  std::string s;
  for (size_t i = 0; i < d.size(); ++i) {
    if (i) s += " x ";
    s += std::to_string(d[i]);
  }
  return s.empty() ? "scalar" : s;
}

void walk(hid_t g, const std::string& prefix, std::vector<LsRow>& rows) {
  for (const auto& name : group_children(g)) {
    const std::string path = join_path(prefix, name);
    if (!link_exists(g, name)) {
      rows.push_back({path, "broken-link", "", "", "", "", "", 0});
      continue;
    }
    H5Handle obj = open_object(g, name);
    LsRow r;
    r.path = path;
    r.encoding = attr_string(obj, "encoding-type", "");
    if (object_kind(obj) == H5I_GROUP) {
      r.kind = "group";
      rows.push_back(r);
      walk(obj, path, rows);
      continue;
    }
    if (object_kind(obj) != H5I_DATASET) {
      r.kind = "other";
      rows.push_back(r);
      continue;
    }
    r.kind = "dataset";
    H5Handle t = dataset_type(obj);
    r.dtype = describe_type(t).dtype;
    H5Handle sp = dataset_space(obj);
    r.shape = space_is_null(sp) ? "null" : dims_str(space_dims(sp));
    H5Handle dcpl(H5Dget_create_plist(obj), H5Pclose, "dcpl");
    if (H5Pget_layout(dcpl) == H5D_CHUNKED) {
      int rank = H5Pget_chunk(dcpl, 0, nullptr);
      std::vector<hsize_t> c(static_cast<size_t>(std::max(rank, 0)));
      if (rank > 0) H5Pget_chunk(dcpl, rank, c.data());
      r.chunks = dims_str(c);
    } else if (H5Pget_layout(dcpl) == H5D_CONTIGUOUS) {
      r.chunks = "contiguous";
    } else {
      r.chunks = "compact";
    }
    int nf = H5Pget_nfilters(dcpl);
    for (int i = 0; i < nf; ++i) {
      unsigned flags = 0, cfg = 0;
      size_t ncd = 0;
      H5Z_filter_t id = H5Pget_filter2(dcpl, static_cast<unsigned>(i), &flags, &ncd, nullptr, 0, nullptr, &cfg);
      if (!r.filters.empty()) r.filters += "+";
      r.filters += filter_name(id);
    }
    r.stored = static_cast<double>(H5Dget_storage_size(obj));
    H5Eclear2(H5E_DEFAULT);
    rows.push_back(r);
  }
}

}  // namespace

[[cpp11::register]]
cpp11::writable::list cpp_h5_ls(std::string path) {
  H5Handle fapl(H5Pcreate(H5P_FILE_ACCESS), H5Pclose, "H5Pcreate");
#if H5_VERSION_GE(1, 12, 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 10 && H5_VERS_RELEASE >= 7)
  H5Pset_file_locking(fapl, false, true);
#endif
  H5Handle f(H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl), H5Fclose,
             ("cannot open '" + path + "' as an HDF5 file").c_str());
  std::vector<LsRow> rows;
  walk(f, "", rows);
  cpp11::writable::strings p, k, e, d, s, c, fl;
  cpp11::writable::doubles st;
  for (const auto& r : rows) {
    p.push_back(r.path);
    k.push_back(r.kind);
    e.push_back(r.encoding);
    d.push_back(r.dtype);
    s.push_back(r.shape);
    c.push_back(r.chunks);
    fl.push_back(r.filters);
    st.push_back(r.stored);
  }
  cpp11::writable::list out({p, k, e, d, s, c, fl, st});
  out.names() = {"path", "kind", "encoding", "dtype", "shape", "chunks", "filters", "stored_bytes"};
  return out;
}

[[cpp11::register]]
SEXP cpp_w_open(std::string path, int threads) { return w_open(path, threads); }

[[cpp11::register]]
void cpp_w_close(SEXP xp) { w_close(xp); }

[[cpp11::register]]
void cpp_w_group(SEXP xp, std::string path, std::string enc, std::string ver) {
  w_group(xp, path, enc, ver);
}

[[cpp11::register]]
void cpp_w_attr(SEXP xp, std::string path, std::string name, SEXP value, std::string kind) {
  w_attr(xp, path, name, value, kind);
}

[[cpp11::register]]
std::string cpp_w_array(SEXP xp, std::string path, SEXP x, SEXP shape, bool colmajor, std::string hint,
                        std::string enc, std::string ver, SEXP comp, bool unlimited) {
  return w_array(xp, path, x, shape, colmajor, hint, enc, ver, comp, unlimited);
}

[[cpp11::register]]
void cpp_w_null(SEXP xp, std::string path) { w_null(xp, path); }

[[cpp11::register]]
SEXP cpp_w_sparse(SEXP xp, std::string path, std::string fmt, SEXP shape, SEXP p, SEXP idx, SEXP x,
                  SEXP hints, SEXP comp) {
  return w_sparse(xp, path, fmt, shape, p, idx, x, hints, comp);
}

[[cpp11::register]]
void cpp_w_recarray(SEXP xp, std::string path, SEXP cols, SEXP names, SEXP hints) {
  w_recarray(xp, path, cols, names, hints);
}
