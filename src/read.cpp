// Element readers: turn anndata-encoded HDF5 objects into R objects.
#include "common.h"

#include <algorithm>
#include <type_traits>
#include <unordered_set>
#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

namespace annr {

namespace {

struct ElemOpts {
  bool colmajor = true;  // R dims equal the on-disk (Python) dims
};

struct SlotSel {
  bool read = true;
  bool all = true;
  std::vector<std::string> keys;
  bool wants(const std::string& k) const {
    if (!read) return false;
    if (all) return true;
    return std::find(keys.begin(), keys.end(), k) != keys.end();
  }
};

cpp11::sexp read_elem(hid_t loc, const std::string& name, const std::string& path, Ctx& ctx,
                      const ElemOpts& o);

size_t product(const std::vector<hsize_t>& v) {
  size_t p = 1;
  for (hsize_t x : v) p *= static_cast<size_t>(x);
  return p;
}

std::vector<hsize_t> reversed(std::vector<hsize_t> v) {
  std::reverse(v.begin(), v.end());
  return v;
}

// Position in a column-major array (dims D) of the element with C-order index k.
size_t corder_to_forder(size_t k, const std::vector<hsize_t>& D) {
  const size_t rank = D.size();
  size_t f = 0, stride = 1;
  std::vector<size_t> idx(rank);
  for (size_t j = rank; j-- > 0;) {
    idx[j] = k % D[j];
    k /= D[j];
  }
  for (size_t j = 0; j < rank; ++j) {
    f += idx[j] * stride;
    stride *= D[j];
  }
  return f;
}

// ---------------------------------------------------------------------------
// Strings

bool valid_utf8(const unsigned char* s, size_t n) {
  size_t i = 0;
  while (i < n) {
    unsigned char c = s[i];
    if (c < 0x80) {
      ++i;
      continue;
    }
    size_t extra;
    if ((c & 0xE0) == 0xC0) {
      if (c < 0xC2) return false;
      extra = 1;
    } else if ((c & 0xF0) == 0xE0) {
      extra = 2;
    } else if ((c & 0xF8) == 0xF0) {
      if (c > 0xF4) return false;
      extra = 3;
    } else {
      return false;
    }
    if (i + extra >= n) return false;
    for (size_t k = 1; k <= extra; ++k) {
      if ((s[i + k] & 0xC0) != 0x80) return false;
    }
    i += extra + 1;
  }
  return true;
}

// Must run inside unwind_protect.
SEXP mkchar(const char* s, size_t n) {
  if (n > static_cast<size_t>(INT_MAX)) Rf_error("string too long");
  const unsigned char* u = reinterpret_cast<const unsigned char*>(s);
  bool ascii = true;
  for (size_t i = 0; i < n; ++i) {
    if (u[i] & 0x80) {
      ascii = false;
      break;
    }
  }
  cetype_t enc = (ascii || valid_utf8(u, n)) ? CE_UTF8 : CE_LATIN1;
  return Rf_mkCharLenCE(s, static_cast<int>(n), enc);
}

size_t fixed_len(const char* s, size_t size, H5T_str_t pad) {
  size_t len = strnlen(s, size);
  if (pad == H5T_STR_SPACEPAD) {
    while (len > 0 && s[len - 1] == ' ') --len;
  }
  return len;
}

struct VlenBuffer {
  hid_t memtype;
  hid_t space;
  std::vector<char*> ptrs;
  ~VlenBuffer() {
    if (!ptrs.empty()) {
#if H5_VERSION_GE(1, 12, 0)
      H5Treclaim(memtype, space, H5P_DEFAULT, ptrs.data());
#else
      H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, ptrs.data());
#endif
    }
  }
};

cpp11::sexp read_strings(hid_t dset, const TypeDesc& t, bool colmajor) {
  H5Handle space = dataset_space(dset);
  if (space_is_null(space)) return cpp11::sexp(alloc_vector(STRSXP, 0));
  std::vector<hsize_t> dims = space_dims(space);
  const size_t rank = dims.size();
  const size_t n = rank == 0 ? 1 : product(dims);
  cpp11::sexp out(alloc_vector(STRSXP, static_cast<R_xlen_t>(n)));
  const bool permute = colmajor && rank >= 2;
  auto pos = [&](size_t k) { return permute ? corder_to_forder(k, dims) : k; };

  if (t.is_vlen_string) {
    H5Handle mt(H5Tcopy(H5T_C_S1), H5Tclose, "H5Tcopy");
    H5Tset_size(mt, H5T_VARIABLE);
    H5Tset_cset(mt, t.cset);
    if (rank == 1 && n > 0) {
      const size_t B = size_t(1) << 18;
      for (size_t s = 0; s < n; s += B) {
        hsize_t start = s, count = std::min(B, n - s);
        H5Handle fs = dataset_space(dset);
        if (H5Sselect_hyperslab(fs, H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0) {
          h5_fail("H5Sselect_hyperslab");
        }
        H5Handle ms(H5Screate_simple(1, &count, nullptr), H5Sclose, "H5Screate_simple");
        VlenBuffer buf{mt, ms, std::vector<char*>(count, nullptr)};
        if (H5Dread(dset, mt, ms, fs, H5P_DEFAULT, buf.ptrs.data()) < 0) h5_fail("cannot read strings");
        cpp11::unwind_protect([&] {
          for (size_t i = 0; i < count; ++i) {
            const char* p = buf.ptrs[i];
            SET_STRING_ELT(out, static_cast<R_xlen_t>(s + i), p ? mkchar(p, strlen(p)) : R_BlankString);
          }
        });
      }
    } else if (n > 0) {
      VlenBuffer buf{mt, space, std::vector<char*>(n, nullptr)};
      if (H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.ptrs.data()) < 0) {
        h5_fail("cannot read strings");
      }
      cpp11::unwind_protect([&] {
        for (size_t i = 0; i < n; ++i) {
          const char* p = buf.ptrs[i];
          SET_STRING_ELT(out, static_cast<R_xlen_t>(pos(i)), p ? mkchar(p, strlen(p)) : R_BlankString);
        }
      });
    }
  } else if (n > 0) {
    H5Handle ftype = dataset_type(dset);
    H5Handle mt(H5Tcopy(ftype), H5Tclose, "H5Tcopy");
    const size_t sz = t.size;
    std::vector<char> buf(sz * n + 1, '\0');
    if (H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) < 0) h5_fail("cannot read strings");
    cpp11::unwind_protect([&] {
      for (size_t i = 0; i < n; ++i) {
        const char* p = buf.data() + i * sz;
        SET_STRING_ELT(out, static_cast<R_xlen_t>(pos(i)), mkchar(p, fixed_len(p, sz, t.strpad)));
      }
    });
  }
  if (rank >= 2) set_dim(out, colmajor ? dims : reversed(dims));
  return out;
}

// ---------------------------------------------------------------------------
// Numbers

// Turns a buffer of raw int64 (or uint64) values stored in a REALSXP into an
// integer vector when every value fits, otherwise into doubles in place.
cpp11::sexp narrow_int64(cpp11::sexp x, bool is_unsigned64, Ctx& ctx, const std::string& path) {
  const size_t n = static_cast<size_t>(Rf_xlength(x));
  int64_t* v = reinterpret_cast<int64_t*>(REAL(x));
  std::atomic<bool> fits{true}, big{false};
  const size_t nb = block_count(n, 1 << 16, ctx.nthreads);
  parallel_for(nb, ctx.nthreads, [&](size_t b, int) {
    size_t lo = n * b / nb, hi = n * (b + 1) / nb;
    bool f = true, g = false;
    for (size_t i = lo; i < hi; ++i) {
      if (is_unsigned64) {
        uint64_t u = static_cast<uint64_t>(v[i]);
        if (u > static_cast<uint64_t>(INT32_MAX)) f = false;
        if (u > (uint64_t(1) << 53)) g = true;
      } else {
        int64_t s = v[i];
        if (s > INT32_MAX || s < -INT32_MAX) f = false;
        if (s > (int64_t(1) << 53) || s < -(int64_t(1) << 53)) g = true;
      }
    }
    if (!f) fits = false;
    if (g) big = true;
  });
  if (fits) {
    cpp11::sexp y(alloc_vector(INTSXP, static_cast<R_xlen_t>(n)));
    int* yp = INTEGER(y);
    parallel_for(nb, ctx.nthreads, [&](size_t b, int) {
      size_t lo = n * b / nb, hi = n * (b + 1) / nb;
      for (size_t i = lo; i < hi; ++i) yp[i] = static_cast<int>(v[i]);
    });
    return y;
  }
  // Values beyond 2^53 are not exactly representable as doubles: keep the
  // original 64-bit patterns so that write_h5ad() can restore them exactly.
  cpp11::sexp bits = R_NilValue;
  if (big) {
    bits = cpp11::sexp(alloc_vector(REALSXP, static_cast<R_xlen_t>(n)));
    std::memcpy(REAL(bits), REAL(x), n * sizeof(double));
  }
  double* d = REAL(x);
  parallel_for(nb, ctx.nthreads, [&](size_t b, int) {
    size_t lo = n * b / nb, hi = n * (b + 1) / nb;
    for (size_t i = lo; i < hi; ++i) {
      int64_t s = v[i];
      d[i] = is_unsigned64 ? static_cast<double>(static_cast<uint64_t>(s)) : static_cast<double>(s);
    }
  });
  if (big) {
    cpp11::unwind_protect([&] { Rf_setAttrib(x, Rf_install("h5ad_int64"), bits); });
    ctx.warn("'" + path + "': 64-bit integers beyond 2^53 are shown as approximate doubles "
             "(exact values are kept for write_h5ad())");
  }
  return x;
}

cpp11::sexp read_numeric_array(hid_t dset, const TypeDesc& t, bool colmajor, Ctx& ctx,
                               const std::string& path) {
  H5Handle space = dataset_space(dset);
  std::vector<hsize_t> dims = space_dims(space);
  const size_t rank = dims.size();
  const size_t n = space_is_null(space) ? 0 : (rank == 0 ? 1 : product(dims));
  Dst dst;
  SEXPTYPE st;
  if (t.is_bool) {
    dst = Dst::LGL;
    st = LGLSXP;
  } else if (t.cls == H5T_FLOAT) {
    dst = Dst::F64;
    st = REALSXP;
  } else if (t.cls == H5T_INTEGER) {
    bool small = t.is_signed ? t.size <= 4 : t.size <= 2;
    dst = small ? Dst::I32 : Dst::I64;
    st = small ? INTSXP : REALSXP;
  } else {
    throw std::runtime_error("'" + path + "': unsupported numeric type " + t.dtype);
  }
  cpp11::sexp x(alloc_vector(st, static_cast<R_xlen_t>(n)));
  if (n > 0) {
    void* p = st == REALSXP ? static_cast<void*>(REAL(x))
                            : (st == INTSXP ? static_cast<void*>(INTEGER(x)) : static_cast<void*>(LOGICAL(x)));
    read_numeric(dset, t, dst, p, colmajor, ctx);
  }
  if (dst == Dst::I64) x = narrow_int64(x, !t.is_signed && t.size == 8, ctx, path);
  if (rank >= 2) set_dim(x, colmajor ? dims : reversed(dims));
  return x;
}

cpp11::sexp read_complex_array(hid_t dset, const TypeDesc&, bool colmajor) {
  H5Handle space = dataset_space(dset);
  std::vector<hsize_t> dims = space_dims(space);
  const size_t rank = dims.size();
  const size_t n = rank == 0 ? 1 : product(dims);
  H5Handle mt(H5Tcreate(H5T_COMPOUND, 16), H5Tclose, "H5Tcreate");
  H5Tinsert(mt, "r", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert(mt, "i", 8, H5T_NATIVE_DOUBLE);
  cpp11::sexp x(alloc_vector(CPLXSXP, static_cast<R_xlen_t>(n)));
  if (n > 0) {
    std::vector<double> buf(2 * n);
    if (H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) < 0) h5_fail("cannot read complex data");
    Rcomplex* c = COMPLEX(x);
    const bool permute = colmajor && rank >= 2;
    for (size_t k = 0; k < n; ++k) {
      size_t p = permute ? corder_to_forder(k, dims) : k;
      c[p].r = buf[2 * k];
      c[p].i = buf[2 * k + 1];
    }
  }
  if (rank >= 2) set_dim(x, colmajor ? dims : reversed(dims));
  return x;
}

// ---------------------------------------------------------------------------
// Compound datasets (rec-arrays and pre-0.7 dataframes)

cpp11::sexp make_data_frame(cpp11::sexp cols, const std::vector<std::string>& names, size_t nrow) {
  return cpp11::unwind_protect([&] {
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(names.size())));
    for (size_t i = 0; i < names.size(); ++i) SET_STRING_ELT(nm, static_cast<R_xlen_t>(i), mkchar(names[i].data(), names[i].size()));
    Rf_setAttrib(cols, R_NamesSymbol, nm);
    SEXP rn = PROTECT(Rf_allocVector(INTSXP, nrow == 0 ? 0 : 2));
    if (nrow > 0) {
      INTEGER(rn)[0] = NA_INTEGER;
      INTEGER(rn)[1] = -static_cast<int>(nrow);
    }
    Rf_setAttrib(cols, R_RowNamesSymbol, rn);
    Rf_setAttrib(cols, R_ClassSymbol, Rf_mkString("data.frame"));
    UNPROTECT(2);
    return static_cast<SEXP>(cols);
  });
}

cpp11::sexp read_compound(hid_t dset, const std::string& path, Ctx& ctx) {
  H5Handle ftype = dataset_type(dset);
  H5Handle space = dataset_space(dset);
  std::vector<hsize_t> dims = space_dims(space);
  const size_t n = dims.empty() ? 1 : product(dims);
  const int nm = H5Tget_nmembers(ftype);
  std::vector<std::string> names;
  std::vector<cpp11::sexp> cols;
  for (int i = 0; i < nm; ++i) {
    char* cname = H5Tget_member_name(ftype, static_cast<unsigned>(i));
    std::string name(cname ? cname : "");
    if (cname) H5free_memory(cname);
    H5Handle mtype(H5Tget_member_type(ftype, static_cast<unsigned>(i)), H5Tclose, "member type");
    TypeDesc mt = describe_type(mtype);
    const std::string fpath = path + "/" + name;
    H5T_class_t cls = H5Tget_class(mtype);
    std::vector<hsize_t> adims;
    H5Handle base;
    hid_t elem_type = mtype;
    if (cls == H5T_ARRAY) {
      int ar = H5Tget_array_ndims(mtype);
      adims.resize(static_cast<size_t>(ar));
      H5Tget_array_dims2(mtype, adims.data());
      base = H5Handle(H5Tget_super(mtype), H5Tclose, "array base");
      elem_type = base;
      mt = describe_type(base);
    }
    const size_t k = adims.empty() ? 1 : product(adims);
    // Memory type for this field.
    H5Handle melem;
    size_t esz = 0;
    enum { NUM_F64, NUM_I32, NUM_I64, NUM_U64, BOOL8, STR_V, STR_F, CPLX, SKIP } kind = SKIP;
    if (mt.is_bool) {
      melem = H5Handle(H5Tget_native_type(elem_type, H5T_DIR_ASCEND), H5Tclose, "native");
      kind = BOOL8;
    } else if (mt.cls == H5T_FLOAT) {
      melem = H5Handle(H5Tcopy(H5T_NATIVE_DOUBLE), H5Tclose, "copy");
      kind = NUM_F64;
    } else if (mt.cls == H5T_INTEGER) {
      bool small = mt.is_signed ? mt.size <= 4 : mt.size <= 2;
      if (small) {
        melem = H5Handle(H5Tcopy(H5T_NATIVE_INT32), H5Tclose, "copy");
        kind = NUM_I32;
      } else if (!mt.is_signed && mt.size == 8) {
        melem = H5Handle(H5Tcopy(H5T_NATIVE_UINT64), H5Tclose, "copy");
        kind = NUM_U64;
      } else {
        melem = H5Handle(H5Tcopy(H5T_NATIVE_INT64), H5Tclose, "copy");
        kind = NUM_I64;
      }
    } else if (mt.is_string) {
      if (mt.is_vlen_string) {
        melem = H5Handle(H5Tcopy(H5T_C_S1), H5Tclose, "copy");
        H5Tset_size(melem, H5T_VARIABLE);
        H5Tset_cset(melem, mt.cset);
        kind = STR_V;
      } else {
        melem = H5Handle(H5Tcopy(elem_type), H5Tclose, "copy");
        kind = STR_F;
      }
    } else if (mt.is_complex) {
      melem = H5Handle(H5Tcreate(H5T_COMPOUND, 16), H5Tclose, "create");
      H5Tinsert(melem, "r", 0, H5T_NATIVE_DOUBLE);
      H5Tinsert(melem, "i", 8, H5T_NATIVE_DOUBLE);
      kind = CPLX;
    }
    if (kind == SKIP) {
      ctx.warn("'" + fpath + "': field type " + mt.dtype + " is not supported and was skipped");
      continue;
    }
    H5Handle mfield;
    if (!adims.empty()) {
      mfield = H5Handle(H5Tarray_create2(melem, static_cast<unsigned>(adims.size()), adims.data()),
                        H5Tclose, "array create");
    } else {
      mfield = H5Handle(H5Tcopy(melem), H5Tclose, "copy");
    }
    esz = H5Tget_size(mfield);
    H5Handle rec(H5Tcreate(H5T_COMPOUND, esz), H5Tclose, "H5Tcreate");
    if (H5Tinsert(rec, name.c_str(), 0, mfield) < 0) h5_fail("H5Tinsert");
    std::vector<uint8_t> buf(std::max<size_t>(1, esz * n));
    if (H5Dread(dset, rec, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) < 0) {
      h5_fail("cannot read field '" + fpath + "'");
    }
    const size_t total = n * k;
    cpp11::sexp col;
    switch (kind) {
      case NUM_F64: {
        col = cpp11::sexp(alloc_vector(REALSXP, static_cast<R_xlen_t>(total)));
        double* d = REAL(col);
        for (size_t r = 0; r < n; ++r)
          for (size_t c = 0; c < k; ++c) std::memcpy(&d[c * n + r], buf.data() + (r * k + c) * 8, 8);
        break;
      }
      case NUM_I32: {
        col = cpp11::sexp(alloc_vector(INTSXP, static_cast<R_xlen_t>(total)));
        int* d = INTEGER(col);
        for (size_t r = 0; r < n; ++r)
          for (size_t c = 0; c < k; ++c) std::memcpy(&d[c * n + r], buf.data() + (r * k + c) * 4, 4);
        break;
      }
      case NUM_I64:
      case NUM_U64: {
        col = cpp11::sexp(alloc_vector(REALSXP, static_cast<R_xlen_t>(total)));
        double* d = REAL(col);
        for (size_t r = 0; r < n; ++r)
          for (size_t c = 0; c < k; ++c) std::memcpy(&d[c * n + r], buf.data() + (r * k + c) * 8, 8);
        col = narrow_int64(col, kind == NUM_U64, ctx, fpath);
        break;
      }
      case BOOL8: {
        const size_t bs = H5Tget_size(melem);
        col = cpp11::sexp(alloc_vector(LGLSXP, static_cast<R_xlen_t>(total)));
        int* d = LOGICAL(col);
        for (size_t r = 0; r < n; ++r)
          for (size_t c = 0; c < k; ++c) d[c * n + r] = buf[(r * k + c) * bs] != 0;
        break;
      }
      case CPLX: {
        col = cpp11::sexp(alloc_vector(CPLXSXP, static_cast<R_xlen_t>(total)));
        Rcomplex* d = COMPLEX(col);
        for (size_t r = 0; r < n; ++r)
          for (size_t c = 0; c < k; ++c) {
            std::memcpy(&d[c * n + r].r, buf.data() + (r * k + c) * 16, 8);
            std::memcpy(&d[c * n + r].i, buf.data() + (r * k + c) * 16 + 8, 8);
          }
        break;
      }
      case STR_V: {
        col = cpp11::sexp(alloc_vector(STRSXP, static_cast<R_xlen_t>(total)));
        struct Reclaim {
          hid_t t, s;
          void* b;
          ~Reclaim() {
#if H5_VERSION_GE(1, 12, 0)
            H5Treclaim(t, s, H5P_DEFAULT, b);
#else
            H5Dvlen_reclaim(t, s, H5P_DEFAULT, b);
#endif
          }
        } reclaim{rec, space, buf.data()};
        cpp11::unwind_protect([&] {
          for (size_t r = 0; r < n; ++r)
            for (size_t c = 0; c < k; ++c) {
              char* p;
              std::memcpy(&p, buf.data() + r * esz + c * sizeof(char*), sizeof(char*));
              SET_STRING_ELT(col, static_cast<R_xlen_t>(c * n + r), p ? mkchar(p, strlen(p)) : R_BlankString);
            }
        });
        break;
      }
      case STR_F: {
        const size_t ss = H5Tget_size(melem);
        H5T_str_t pad = H5Tget_strpad(melem);
        col = cpp11::sexp(alloc_vector(STRSXP, static_cast<R_xlen_t>(total)));
        cpp11::unwind_protect([&] {
          for (size_t r = 0; r < n; ++r)
            for (size_t c = 0; c < k; ++c) {
              const char* p = reinterpret_cast<const char*>(buf.data() + r * esz + c * ss);
              SET_STRING_ELT(col, static_cast<R_xlen_t>(c * n + r), mkchar(p, fixed_len(p, ss, pad)));
            }
        });
        break;
      }
      default:
        break;
    }
    if (!adims.empty()) {
      std::vector<hsize_t> md{static_cast<hsize_t>(n), static_cast<hsize_t>(k)};
      set_dim(col, md);
    }
    ctx.record(fpath, "rec-array-field", "", (kind == STR_V ? "str" : mt.dtype),
               adims.empty() ? "" : std::to_string(k));
    names.push_back(name);
    cols.push_back(col);
  }
  cpp11::sexp lst(alloc_vector(VECSXP, static_cast<R_xlen_t>(cols.size())));
  for (size_t i = 0; i < cols.size(); ++i) SET_VECTOR_ELT(lst, static_cast<R_xlen_t>(i), cols[i]);
  return make_data_frame(lst, names, n);
}

// ---------------------------------------------------------------------------
// Sparse matrices

template <typename V>
size_t canonicalize_sparse(int* p, int* idx, V* x, size_t major, int64_t minor, Ctx& ctx,
                           const std::string& path) {
  std::atomic<bool> bad{false};
  std::vector<uint8_t> unsorted_row(major, 0);
  std::atomic<size_t> n_unsorted{0};
  const size_t nb = block_count(major, 64, ctx.nthreads);
  parallel_for(nb, ctx.nthreads, [&](size_t b, int) {
    size_t lo = major * b / nb, hi = major * (b + 1) / nb, cnt = 0;
    for (size_t m = lo; m < hi; ++m) {
      bool u = false;
      for (int k = p[m]; k < p[m + 1]; ++k) {
        int v = idx[k];
        if (v < 0 || v >= minor) {
          bad = true;
          return;
        }
        if (k > p[m] && v <= idx[k - 1]) u = true;
      }
      if (u) {
        unsorted_row[m] = 1;
        ++cnt;
      }
    }
    n_unsorted += cnt;
  });
  if (bad) throw std::runtime_error("'" + path + "': sparse matrix has out-of-range indices (corrupt file?)");
  if (n_unsorted == 0) return 0;
  // Stable sort of each unsorted row: LSD radix sort on the column index
  // (8-bit digits, as many passes as the minor dimension needs), insertion
  // sort for short rows. Duplicate indices are counted on the way.
  int passes = 1;
  while (passes < 4 && (static_cast<int64_t>(1) << (8 * passes)) < minor) ++passes;
  std::atomic<size_t> n_dups{0};
  parallel_for(nb, ctx.nthreads, [&](size_t b, int) {
    size_t lo = major * b / nb, hi = major * (b + 1) / nb, dups = 0;
    std::vector<uint32_t> ka, kb;
    std::vector<V> va, vb;
    for (size_t m = lo; m < hi; ++m) {
      if (!unsorted_row[m]) continue;
      const int start = p[m], len = p[m + 1] - p[m];
      int* ri = idx + start;
      V* rx = x + start;
      if (len <= 48) {
        for (int i = 1; i < len; ++i) {
          int key = ri[i];
          V val = rx[i];
          int j = i - 1;
          while (j >= 0 && ri[j] > key) {
            ri[j + 1] = ri[j];
            rx[j + 1] = rx[j];
            --j;
          }
          ri[j + 1] = key;
          rx[j + 1] = val;
        }
      } else {
        const size_t n = static_cast<size_t>(len);
        ka.assign(ri, ri + n);
        va.assign(rx, rx + n);
        kb.resize(n);
        vb.resize(n);
        for (int pass = 0; pass < passes; ++pass) {
          const int shift = 8 * pass;
          size_t count[257] = {0};
          for (size_t i = 0; i < n; ++i) ++count[((ka[i] >> shift) & 0xffu) + 1];
          for (int d = 0; d < 256; ++d) count[d + 1] += count[d];
          for (size_t i = 0; i < n; ++i) {
            size_t pos = count[(ka[i] >> shift) & 0xffu]++;
            kb[pos] = ka[i];
            vb[pos] = va[i];
          }
          ka.swap(kb);
          va.swap(vb);
        }
        for (size_t i = 0; i < n; ++i) {
          ri[i] = static_cast<int>(ka[i]);
          rx[i] = va[i];
        }
      }
      for (int k = 1; k < len; ++k) {
        if (ri[k] == ri[k - 1]) ++dups;
      }
    }
    n_dups += dups;
  });
  const size_t dups = n_dups;
  if (dups == 0) return 0;
  // Sum duplicates (scipy semantics), compacting the arrays.
  int w = 0;
  for (size_t m = 0; m < major; ++m) {
    int start = p[m], end = p[m + 1];
    p[m] = w;
    for (int k = start; k < end; ++k) {
      if (w > p[m] && idx[k] == idx[w - 1]) {
        if constexpr (std::is_same_v<V, int>) {
          x[w - 1] = (x[w - 1] != 0 || x[k] != 0) ? 1 : 0;
        } else {
          x[w - 1] += x[k];
        }
      } else {
        idx[w] = idx[k];
        x[w] = x[k];
        ++w;
      }
    }
  }
  p[major] = w;
  ctx.warn("'" + path + "': " + std::to_string(dups) +
           " duplicate sparse entries were summed (scipy semantics)");
  return dups;
}

cpp11::sexp truncate_vector(cpp11::sexp x, size_t n) {
  if (static_cast<size_t>(Rf_xlength(x)) == n) return x;
  return cpp11::sexp(cpp11::safe[Rf_xlengthgets](x, static_cast<R_xlen_t>(n)));
}

cpp11::sexp read_sparse(hid_t g, const std::string& path, bool csr, std::vector<int64_t> shape,
                        Ctx& ctx) {
  if (shape.size() != 2) throw std::runtime_error("'" + path + "': sparse matrix needs a 2-D shape");
  const int64_t nrow = shape[0], ncol = shape[1];
  if (nrow < 0 || ncol < 0 || nrow > INT_MAX || ncol > INT_MAX) {
    throw std::runtime_error("'" + path + "': sparse matrix dimensions exceed R's limits");
  }
  const int64_t major = csr ? nrow : ncol, minor = csr ? ncol : nrow;
  H5Handle dd = open_dataset(g, "data");
  H5Handle di = open_dataset(g, "indices");
  H5Handle dp = open_dataset(g, "indptr");
  TypeDesc td = describe_type(dataset_type(dd));
  TypeDesc ti = describe_type(dataset_type(di));
  TypeDesc tp = describe_type(dataset_type(dp));
  if (!td.is_numeric || ti.cls != H5T_INTEGER || tp.cls != H5T_INTEGER) {
    throw std::runtime_error("'" + path + "': unsupported sparse component types");
  }
  auto npoints = [](hid_t d) {
    H5Handle s = dataset_space(d);
    return static_cast<size_t>(H5Sget_simple_extent_npoints(s));
  };
  const size_t ndata = npoints(dd), nidx = npoints(di), nptr = npoints(dp);
  if (nptr != static_cast<size_t>(major) + 1) {
    throw std::runtime_error("'" + path + "': indptr has length " + std::to_string(nptr) +
                             ", expected " + std::to_string(major + 1));
  }
  std::vector<int64_t> indptr(nptr);
  read_numeric(dp, tp, Dst::I64, indptr.data(), false, ctx);
  if (indptr[0] != 0) throw std::runtime_error("'" + path + "': indptr must start at 0");
  for (size_t m = 0; m < static_cast<size_t>(major); ++m) {
    if (indptr[m + 1] < indptr[m]) throw std::runtime_error("'" + path + "': indptr is not non-decreasing");
  }
  const int64_t nnz = indptr[static_cast<size_t>(major)];
  if (nnz > INT_MAX) {
    throw std::runtime_error("'" + path + "': " + std::to_string(nnz) +
                             " non-zero values exceed the 2^31-1 limit of R's Matrix package");
  }
  if (static_cast<size_t>(nnz) > ndata || static_cast<size_t>(nnz) > nidx) {
    throw std::runtime_error("'" + path + "': data/indices are shorter than indptr implies");
  }
  cpp11::sexp p(alloc_vector(INTSXP, major + 1));
  int* pp = INTEGER(p);
  for (size_t m = 0; m < nptr; ++m) pp[m] = static_cast<int>(indptr[m]);
  indptr.clear();
  indptr.shrink_to_fit();

  ProfileScope prof_all{"sparse", static_cast<size_t>(nnz)};
  const bool lgl = td.is_bool;
  cpp11::sexp idx, x;
  {
    ProfileScope prof_alloc{"alloc", static_cast<size_t>(nnz)};
    idx = cpp11::sexp(alloc_vector(INTSXP, nnz));
    x = cpp11::sexp(alloc_vector(lgl ? LGLSXP : REALSXP, nnz));
  }
  if (nnz > 0) {
    if (nidx == static_cast<size_t>(nnz)) {
      read_numeric(di, ti, Dst::I32, INTEGER(idx), false, ctx);
    } else {
      read_numeric_range(di, ti, Dst::I32, INTEGER(idx), 0, static_cast<hsize_t>(nnz), ctx);
    }
    void* xp = lgl ? static_cast<void*>(LOGICAL(x)) : static_cast<void*>(REAL(x));
    if (ndata == static_cast<size_t>(nnz)) {
      read_numeric(dd, td, lgl ? Dst::LGL : Dst::F64, xp, false, ctx);
    } else {
      read_numeric_range(dd, td, lgl ? Dst::LGL : Dst::F64, xp, 0, static_cast<hsize_t>(nnz), ctx);
    }
    ProfileScope prof_canon{"canon", static_cast<size_t>(nnz)};
    size_t dups = lgl ? canonicalize_sparse<int>(pp, INTEGER(idx), LOGICAL(x), static_cast<size_t>(major), minor, ctx, path)
                      : canonicalize_sparse<double>(pp, INTEGER(idx), REAL(x), static_cast<size_t>(major), minor, ctx, path);
    if (dups > 0) {
      idx = truncate_vector(idx, static_cast<size_t>(nnz) - dups);
      x = truncate_vector(x, static_cast<size_t>(nnz) - dups);
    }
  }
  const char* cls = csr ? (lgl ? "lgRMatrix" : "dgRMatrix") : (lgl ? "lgCMatrix" : "dgCMatrix");
  cpp11::sexp obj(make_s4(cls));
  cpp11::unwind_protect([&] {
    R_do_slot_assign(obj, Rf_install("p"), p);
    R_do_slot_assign(obj, Rf_install(csr ? "j" : "i"), idx);
    R_do_slot_assign(obj, Rf_install("x"), x);
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim)[0] = static_cast<int>(nrow);
    INTEGER(dim)[1] = static_cast<int>(ncol);
    R_do_slot_assign(obj, Rf_install("Dim"), dim);
    UNPROTECT(1);
  });
  ctx.record(path + "/data", "", "", td.dtype);
  ctx.record(path + "/indices", "", "", ti.dtype);
  ctx.record(path + "/indptr", "", "", tp.dtype);
  return obj;
}

// ---------------------------------------------------------------------------
// Categoricals

cpp11::sexp as_character(SEXP x) {
  return cpp11::unwind_protect([&] {
    SEXP call = PROTECT(Rf_lang2(Rf_install("as.character"), x));
    SEXP res = Rf_eval(call, R_BaseEnv);
    UNPROTECT(1);
    return res;
  });
}

cpp11::sexp build_factor(cpp11::sexp codes, cpp11::sexp cats, bool ordered, const std::string& path) {
  cpp11::sexp levels = TYPEOF(cats) == STRSXP ? cats : as_character(cats);
  const R_xlen_t ncat = Rf_xlength(levels);
  int* c = INTEGER(codes);
  const R_xlen_t n = Rf_xlength(codes);
  for (R_xlen_t i = 0; i < n; ++i) {
    int v = c[i];
    if (v == -1) {
      c[i] = NA_INTEGER;
    } else if (v < 0 || v >= ncat) {
      throw std::runtime_error("'" + path + "': categorical code out of range");
    } else {
      c[i] = v + 1;
    }
  }
  cpp11::unwind_protect([&] {
    Rf_setAttrib(codes, R_LevelsSymbol, levels);
    if (ordered) {
      SEXP cl = PROTECT(Rf_allocVector(STRSXP, 2));
      SET_STRING_ELT(cl, 0, Rf_mkChar("ordered"));
      SET_STRING_ELT(cl, 1, Rf_mkChar("factor"));
      Rf_setAttrib(codes, R_ClassSymbol, cl);
      UNPROTECT(1);
    } else {
      Rf_setAttrib(codes, R_ClassSymbol, Rf_mkString("factor"));
    }
  });
  return codes;
}

cpp11::sexp read_codes(hid_t dset, Ctx& ctx) {
  TypeDesc t = describe_type(dataset_type(dset));
  if (t.cls != H5T_INTEGER) throw std::runtime_error("categorical codes must be integers");
  H5Handle space = dataset_space(dset);
  const size_t n = static_cast<size_t>(H5Sget_simple_extent_npoints(space));
  cpp11::sexp codes(alloc_vector(INTSXP, static_cast<R_xlen_t>(n)));
  if (n > 0) read_numeric(dset, t, Dst::I32, INTEGER(codes), false, ctx);
  return codes;
}

cpp11::sexp read_categorical(hid_t g, const std::string& path, Ctx& ctx) {
  H5Handle cd = open_dataset(g, "codes");
  cpp11::sexp codes = read_codes(cd, ctx);
  ctx.record(path + "/codes", attr_string(cd, "encoding-type", ""), "",
             describe_type(dataset_type(cd)).dtype);
  cpp11::sexp cats = read_elem(g, "categories", path + "/categories", ctx, ElemOpts{});
  bool ordered = attr_bool(g, "ordered", false);
  return build_factor(codes, cats, ordered, path);
}

// Columns of 0.1.0 dataframes reference their categories through an
// object reference stored in the "categories" attribute.
cpp11::sexp read_legacy_series(hid_t parent, const std::string& name, const std::string& path,
                               Ctx& ctx) {
  H5Handle d = open_object(parent, name);
  if (object_kind(d) != H5I_DATASET || !has_attr(d, "categories")) {
    return read_elem(parent, name, path, ctx, ElemOpts{});
  }
  H5Handle cats_obj;
  {
    H5Handle attr(H5Aopen(d, "categories", H5P_DEFAULT), H5Aclose, "open categories attribute");
    H5Handle at(H5Aget_type(attr), H5Tclose, "attr type");
#ifndef H5_NO_DEPRECATED_SYMBOLS
    if (H5Tget_class(at) == H5T_REFERENCE) {
      hobj_ref_t ref;
      if (H5Aread(attr, H5T_STD_REF_OBJ, &ref) >= 0) {
        hid_t o = H5Rdereference2(d, H5P_DEFAULT, H5R_OBJECT, &ref);
        if (o >= 0) cats_obj = H5Handle(o, H5Oclose, "dereference");
      }
      H5Eclear2(H5E_DEFAULT);
    }
#endif
  }
  if (cats_obj.get() < 0) {
    std::string alt = "__categories/" + name;
    if (!link_exists(parent, alt)) return read_elem(parent, name, path, ctx, ElemOpts{});
    cats_obj = open_object(parent, alt);
  }
  cpp11::sexp codes = read_codes(d, ctx);
  TypeDesc ct = describe_type(dataset_type(cats_obj));
  cpp11::sexp cats = ct.is_string ? read_strings(cats_obj, ct, true)
                                  : read_numeric_array(cats_obj, ct, true, ctx, path + "/categories");
  bool ordered = attr_bool(cats_obj, "ordered", false);
  ctx.record(path, "categorical", "legacy", describe_type(dataset_type(d)).dtype);
  ctx.record(path + "/categories", "", "", ct.dtype);
  return build_factor(codes, cats, ordered, path);
}

// ---------------------------------------------------------------------------
// Dataframes

// Attaches row names; non-unique or missing names are kept in "h5ad_index".
void attach_row_names(cpp11::sexp df, cpp11::sexp index, size_t nrow) {
  cpp11::sexp idx = TYPEOF(index) == STRSXP ? index : as_character(index);
  bool unique = true;
  {
    std::unordered_set<SEXP> seen;
    seen.reserve(static_cast<size_t>(Rf_xlength(idx)) * 2 + 1);
    for (R_xlen_t i = 0; i < Rf_xlength(idx); ++i) {
      SEXP s = STRING_ELT(idx, i);
      if (s == NA_STRING || !seen.insert(s).second) {
        unique = false;
        break;
      }
    }
  }
  cpp11::unwind_protect([&] {
    if (unique && Rf_xlength(idx) > 0) {
      Rf_setAttrib(df, R_RowNamesSymbol, idx);
    } else {
      SEXP rn = PROTECT(Rf_allocVector(INTSXP, nrow == 0 ? 0 : 2));
      if (nrow > 0) {
        INTEGER(rn)[0] = NA_INTEGER;
        INTEGER(rn)[1] = -static_cast<int>(nrow);
      }
      Rf_setAttrib(df, R_RowNamesSymbol, rn);
      UNPROTECT(1);
      if (Rf_xlength(idx) > 0) Rf_setAttrib(df, Rf_install("h5ad_index"), idx);
    }
  });
}

cpp11::sexp read_dataframe(hid_t g, const std::string& path, const std::string& ver, Ctx& ctx) {
  std::string idx_key = attr_string(g, "_index", "_index");
  std::vector<std::string> cols = attr_strings(g, "column-order");
  const bool legacy = ver == "0.1.0";
  cpp11::sexp index = R_NilValue;
  if (link_exists(g, idx_key)) {
    index = legacy ? read_legacy_series(g, idx_key, join_path(path, idx_key), ctx)
                   : read_elem(g, idx_key, join_path(path, idx_key), ctx, ElemOpts{});
  }
  std::vector<cpp11::sexp> vals;
  std::vector<std::string> names;
  long long nrow = index == R_NilValue ? -1 : static_cast<long long>(Rf_xlength(index));
  for (const auto& c : cols) {
    if (!link_exists(g, c)) {
      ctx.warn("'" + path + "': column '" + c + "' listed in column-order is missing");
      continue;
    }
    cpp11::sexp v = legacy ? read_legacy_series(g, c, join_path(path, c), ctx)
                           : read_elem(g, c, join_path(path, c), ctx, ElemOpts{});
    long long len = Rf_isNull(Rf_getAttrib(v, R_DimSymbol)) ? static_cast<long long>(Rf_xlength(v))
                                                            : INTEGER(Rf_getAttrib(v, R_DimSymbol))[0];
    if (nrow < 0) nrow = len;
    if (len != nrow) {
      throw std::runtime_error("'" + path + "': column '" + c + "' has " + std::to_string(len) +
                               " rows, expected " + std::to_string(nrow));
    }
    vals.push_back(v);
    names.push_back(c);
  }
  if (nrow < 0) nrow = 0;
  cpp11::sexp lst(alloc_vector(VECSXP, static_cast<R_xlen_t>(vals.size())));
  for (size_t i = 0; i < vals.size(); ++i) SET_VECTOR_ELT(lst, static_cast<R_xlen_t>(i), vals[i]);
  cpp11::sexp df = make_data_frame(lst, names, static_cast<size_t>(nrow));
  if (index != R_NilValue) attach_row_names(df, index, static_cast<size_t>(nrow));
  ctx.record(path, "dataframe", ver, "", idx_key);
  return df;
}

// ---------------------------------------------------------------------------
// Nullable arrays

cpp11::sexp read_nullable(hid_t g, const std::string& path, const std::string& enc, Ctx& ctx) {
  H5Handle dv = open_dataset(g, "values");
  H5Handle dm = open_dataset(g, "mask");
  TypeDesc tv = describe_type(dataset_type(dv));
  TypeDesc tm = describe_type(dataset_type(dm));
  cpp11::sexp vals = tv.is_string ? read_strings(dv, tv, true)
                                  : read_numeric_array(dv, tv, true, ctx, path + "/values");
  cpp11::sexp mask = read_numeric_array(dm, tm, true, ctx, path + "/mask");
  if (Rf_xlength(mask) != Rf_xlength(vals)) throw std::runtime_error("'" + path + "': mask length mismatch");
  const R_xlen_t n = Rf_xlength(vals);
  const int* mk = TYPEOF(mask) == LGLSXP ? LOGICAL(mask) : INTEGER(mask);
  switch (TYPEOF(vals)) {
    case INTSXP:
      for (R_xlen_t i = 0; i < n; ++i) if (mk[i]) INTEGER(vals)[i] = NA_INTEGER;
      break;
    case REALSXP:
      for (R_xlen_t i = 0; i < n; ++i) if (mk[i]) REAL(vals)[i] = NA_REAL;
      break;
    case LGLSXP:
      for (R_xlen_t i = 0; i < n; ++i) if (mk[i]) LOGICAL(vals)[i] = NA_LOGICAL;
      break;
    case STRSXP:
      for (R_xlen_t i = 0; i < n; ++i) if (mk[i]) SET_STRING_ELT(vals, i, NA_STRING);
      break;
    default:
      break;
  }
  ctx.record(path + "/values", "", "", tv.dtype);
  std::string extra = enc == "nullable-string-array" ? attr_string(g, "na-value", "NA") : "";
  ctx.record(path, enc, attr_string(g, "encoding-version", "0.1.0"), tv.dtype, extra);
  return vals;
}

// ---------------------------------------------------------------------------
// Generic dispatch

cpp11::sexp named_list(std::vector<std::string>& names, std::vector<cpp11::sexp>& vals) {
  cpp11::sexp lst(alloc_vector(VECSXP, static_cast<R_xlen_t>(vals.size())));
  for (size_t i = 0; i < vals.size(); ++i) SET_VECTOR_ELT(lst, static_cast<R_xlen_t>(i), vals[i]);
  cpp11::unwind_protect([&] {
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(names.size())));
    for (size_t i = 0; i < names.size(); ++i) {
      SET_STRING_ELT(nm, static_cast<R_xlen_t>(i), mkchar(names[i].data(), names[i].size()));
    }
    Rf_setAttrib(lst, R_NamesSymbol, nm);
    UNPROTECT(1);
  });
  return lst;
}

cpp11::sexp read_dict(hid_t g, const std::string& path, Ctx& ctx, const ElemOpts& o,
                      const SlotSel* sel = nullptr) {
  std::vector<std::string> keys = group_children(g);
  std::vector<std::string> names;
  std::vector<cpp11::sexp> vals;
  for (const auto& k : keys) {
    if (sel != nullptr && !sel->wants(k)) continue;
    if (!link_exists(g, k)) {
      ctx.warn("'" + join_path(path, k) + "': link cannot be resolved and was skipped");
      continue;
    }
    names.push_back(k);
    vals.push_back(read_elem(g, k, join_path(path, k), ctx, o));
  }
  return named_list(names, vals);
}

cpp11::sexp read_awkward(hid_t g, const std::string& path, Ctx& ctx) {
  cpp11::sexp buffers = read_dict(g, path, ctx, ElemOpts{});
  std::string form = attr_string(g, "form", "");
  std::vector<int64_t> len = attr_int64s(g, "length");
  cpp11::unwind_protect([&] {
    SEXP f = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(f, 0, mkchar(form.data(), form.size()));
    Rf_setAttrib(buffers, Rf_install("form"), f);
    SEXP l = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(l)[0] = len.empty() ? NA_REAL : static_cast<double>(len[0]);
    Rf_setAttrib(buffers, Rf_install("length"), l);
    Rf_setAttrib(buffers, R_ClassSymbol, Rf_mkString("h5ad_awkward"));
    UNPROTECT(2);
  });
  return buffers;
}

cpp11::sexp read_anndata(hid_t g, const std::string& prefix, Ctx& ctx, bool x_colmajor,
                         const std::vector<SlotSel>* sels, bool raw_only);

cpp11::sexp read_dataset_value(hid_t d, const std::string& path, const std::string& enc,
                               const TypeDesc& t, Ctx& ctx, const ElemOpts& o) {
  if (t.is_string) return read_strings(d, t, o.colmajor);
  if (t.is_complex) return read_complex_array(d, t, o.colmajor);
  if (t.is_compound) return read_compound(d, path, ctx);
  if (t.is_numeric) return read_numeric_array(d, t, o.colmajor, ctx, path);
  ctx.warn("'" + path + "': dtype " + t.dtype + " (encoding '" + enc + "') is not supported; returned NULL");
  return cpp11::sexp(R_NilValue);
}

cpp11::sexp read_elem_impl(hid_t loc, const std::string& name, const std::string& path, Ctx& ctx,
                           const ElemOpts& o);

// Adds the element path to error messages (once, at the innermost level).
cpp11::sexp read_elem(hid_t loc, const std::string& name, const std::string& path, Ctx& ctx,
                      const ElemOpts& o) {
  try {
    return read_elem_impl(loc, name, path, ctx, o);
  } catch (const cpp11::unwind_exception&) {
    throw;
  } catch (const std::exception& e) {
    std::string msg = e.what();
    if (msg.rfind("while reading '", 0) == 0) throw;
    throw std::runtime_error("while reading '" + path + "': " + msg);
  }
}

cpp11::sexp read_elem_impl(hid_t loc, const std::string& name, const std::string& path, Ctx& ctx,
                           const ElemOpts& o) {
  H5Handle obj = open_object(loc, name);
  const H5I_type_t kind = object_kind(obj);
  const std::string enc = attr_string(obj, "encoding-type", "");
  const std::string ver = attr_string(obj, "encoding-version", "");

  if (kind == H5I_GROUP) {
    if (enc == "csr_matrix" || enc == "csc_matrix") {
      ctx.record(path, enc, ver, "");
      return read_sparse(obj, path, enc == "csr_matrix", attr_int64s(obj, "shape"), ctx);
    }
    if (enc.empty() && has_attr(obj, "h5sparse_format")) {
      std::string fmt = attr_string(obj, "h5sparse_format", "csr");
      ctx.record(path, fmt == "csc" ? "csc_matrix" : "csr_matrix", "legacy", "");
      return read_sparse(obj, path, fmt != "csc", attr_int64s(obj, "h5sparse_shape"), ctx);
    }
    if (enc == "categorical") {
      ctx.record(path, enc, ver, "");
      return read_categorical(obj, path, ctx);
    }
    if (enc == "dataframe") return read_dataframe(obj, path, ver, ctx);
    if (enc == "nullable-integer" || enc == "nullable-boolean" || enc == "nullable-string-array") {
      return read_nullable(obj, path, enc, ctx);
    }
    if (enc == "awkward-array") {
      ctx.record(path, enc, ver, "");
      return read_awkward(obj, path, ctx);
    }
    if (enc == "anndata" || enc == "raw") {
      ctx.record(path, enc, ver, "");
      return read_anndata(obj, path, ctx, o.colmajor, nullptr, enc == "raw");
    }
    if (!enc.empty() && enc != "dict") {
      ctx.warn("'" + path + "': unknown encoding '" + enc + "' read as a dict");
    }
    ctx.record(path, enc.empty() ? "" : "dict", ver, "");
    return read_dict(obj, path, ctx, ElemOpts{});
  }

  if (kind == H5I_DATASET) {
    H5Handle ftype = dataset_type(obj);
    TypeDesc t = describe_type(ftype);
    if (enc == "null") {
      ctx.record(path, enc, ver, "");
      return cpp11::sexp(R_NilValue);
    }
    H5Handle space = dataset_space(obj);
    const bool null_space = space_is_null(space);
    std::vector<hsize_t> dims = null_space ? std::vector<hsize_t>{} : space_dims(space);
    const size_t rec = ctx.types.size();
    ctx.record(path, enc, ver, t.dtype);
    if (null_space) return cpp11::sexp(R_NilValue);
    cpp11::sexp v = read_dataset_value(obj, path, enc, t, ctx, o);
    if (enc.empty()) {
      // Map legacy (encoding-less) datasets to what Python anndata returns:
      // 0-d datasets are scalars, compound datasets are rec-arrays, and
      // pre-0.7 files stored scalar strings as length-1 fixed-length arrays.
      std::string& e = ctx.types[rec].encoding;
      if (dims.empty()) {
        e = t.is_string ? "string" : "numeric-scalar";
      } else if (t.is_compound && !t.is_complex) {
        e = "rec-array";
      } else if (t.is_string && !t.is_vlen_string && dims.size() == 1 && dims[0] == 1) {
        e = "string";
      } else if (t.is_string) {
        e = "string-array";
      } else if (t.is_numeric || t.is_complex) {
        e = "array";
      }
    }
    return v;
  }
  ctx.warn("'" + path + "': object is neither a group nor a dataset and was skipped");
  return cpp11::sexp(R_NilValue);
}

// ---------------------------------------------------------------------------
// AnnData containers

const char* const SLOTS[] = {"X", "obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns", "raw"};
constexpr int NSLOTS = 10;

cpp11::sexp read_mapping_slot(hid_t g, const std::string& slot, const std::string& path, Ctx& ctx,
                              const ElemOpts& o, const SlotSel* sel) {
  H5Handle obj = open_object(g, slot);
  if (object_kind(obj) == H5I_DATASET) {
    // Pre-0.7 obsm/varm: one compound dataset whose fields are the entries.
    H5Handle ftype = dataset_type(obj);
    TypeDesc t = describe_type(ftype);
    ctx.record(path, "", "legacy", t.dtype);
    if (t.is_compound) return read_compound(obj, path, ctx);
    return read_dataset_value(obj, path, "", t, ctx, o);
  }
  std::string enc = attr_string(obj, "encoding-type", "");
  ctx.record(path, enc.empty() ? "" : enc, attr_string(obj, "encoding-version", ""), "");
  return read_dict(obj, path, ctx, o, sel);
}

cpp11::sexp read_anndata(hid_t g, const std::string& prefix, Ctx& ctx, bool x_colmajor,
                         const std::vector<SlotSel>* sels, bool raw_only) {
  std::vector<std::string> names;
  std::vector<cpp11::sexp> vals;
  ElemOpts xo{x_colmajor}, plain{};
  for (int s = 0; s < NSLOTS; ++s) {
    const std::string slot = SLOTS[s];
    const SlotSel* sel = sels ? &(*sels)[static_cast<size_t>(s)] : nullptr;
    names.push_back(slot);
    if (raw_only && slot != "X" && slot != "var" && slot != "varm") {
      vals.emplace_back(R_NilValue);
      continue;
    }
    if ((sel && !sel->read) || !link_exists(g, slot)) {
      vals.emplace_back(R_NilValue);
      continue;
    }
    const std::string path = join_path(prefix, slot);
    if (slot == "X") {
      vals.push_back(read_elem(g, slot, path, ctx, xo));
    } else if (slot == "obs" || slot == "var") {
      H5Handle obj = open_object(g, slot);
      if (object_kind(obj) == H5I_DATASET) {
        ctx.record(path, "dataframe", "legacy-recarray", "");
        vals.push_back(read_compound(obj, path, ctx));
      } else {
        vals.push_back(read_elem(g, slot, path, ctx, plain));
      }
    } else if (slot == "layers") {
      vals.push_back(read_mapping_slot(g, slot, path, ctx, xo, sel));
    } else if (slot == "raw") {
      H5Handle obj = open_object(g, slot);
      if (object_kind(obj) != H5I_GROUP) {
        // nested AnnData objects store a missing raw as a "null" dataset
        vals.push_back(read_elem(g, slot, path, ctx, plain));
        continue;
      }
      ctx.record(path, attr_string(obj, "encoding-type", "raw"), attr_string(obj, "encoding-version", ""), "");
      vals.push_back(read_anndata(obj, path, ctx, x_colmajor, nullptr, true));
    } else {
      vals.push_back(read_mapping_slot(g, slot, path, ctx, plain, sel));
    }
  }
  // Legacy raw (anndata < 0.7): raw.X, raw.var, raw.varm at the top level.
  if (!raw_only && prefix.empty() && (!sels || (*sels)[9].read) && !link_exists(g, "raw") &&
      (link_exists(g, "raw.X") || link_exists(g, "raw.var"))) {
    std::vector<std::string> rn{"X", "obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns", "raw"};
    std::vector<cpp11::sexp> rv(10, cpp11::sexp(R_NilValue));
    if (link_exists(g, "raw.X")) rv[0] = read_elem(g, "raw.X", "raw/X", ctx, xo);
    if (link_exists(g, "raw.var")) {
      H5Handle obj = open_object(g, "raw.var");
      if (object_kind(obj) == H5I_DATASET) {
        ctx.record("raw/var", "dataframe", "legacy-recarray", "");
        rv[2] = read_compound(obj, "raw/var", ctx);
      } else {
        rv[2] = read_elem(g, "raw.var", "raw/var", ctx, plain);
      }
    }
    if (link_exists(g, "raw.varm")) rv[4] = read_mapping_slot(g, "raw.varm", "raw/varm", ctx, plain, nullptr);
    vals[9] = named_list(rn, rv);
    ctx.record("raw", "raw", "legacy", "");
  }
  cpp11::sexp out = named_list(names, vals);
  cpp11::unwind_protect([&] { Rf_setAttrib(out, Rf_install("annr_anndata"), Rf_ScalarLogical(TRUE)); });
  return out;
}

SlotSel parse_sel(SEXP v) {
  SlotSel s;
  if (Rf_isNull(v)) return s;
  if (TYPEOF(v) == LGLSXP && Rf_xlength(v) >= 1) {
    s.read = LOGICAL(v)[0] == TRUE;
    s.all = true;
  } else if (TYPEOF(v) == STRSXP) {
    s.read = true;
    s.all = false;
    for (R_xlen_t i = 0; i < Rf_xlength(v); ++i) s.keys.emplace_back(Rf_translateCharUTF8(STRING_ELT(v, i)));
  }
  return s;
}

struct FdGuard {
  int fd = -1;
  ~FdGuard() {
#ifndef _WIN32
    if (fd >= 0) ::close(fd);
#endif
  }
};

cpp11::sexp string_vec(const std::vector<std::string>& v) {
  cpp11::sexp out(alloc_vector(STRSXP, static_cast<R_xlen_t>(v.size())));
  cpp11::unwind_protect([&] {
    for (size_t i = 0; i < v.size(); ++i) SET_STRING_ELT(out, static_cast<R_xlen_t>(i), mkchar(v[i].data(), v[i].size()));
  });
  return out;
}

H5Handle open_file_ro(const std::string& path) {
  H5Handle fapl(H5Pcreate(H5P_FILE_ACCESS), H5Pclose, "H5Pcreate");
#if H5_VERSION_GE(1, 12, 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 10 && H5_VERS_RELEASE >= 7)
  H5Pset_file_locking(fapl, false, true);
#endif
  return H5Handle(H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl), H5Fclose,
                  ("cannot open '" + path + "' as an HDF5 file").c_str());
}

SEXP opt_get(SEXP opts, const char* nm) {
  SEXP names = Rf_getAttrib(opts, R_NamesSymbol);
  for (R_xlen_t i = 0; i < Rf_xlength(opts); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), nm) == 0) return VECTOR_ELT(opts, i);
  }
  return R_NilValue;
}

void configure_ctx(Ctx& ctx, const std::string& path, SEXP opts) {
  ctx.filename = path;
  SEXP th = opt_get(opts, "threads");
  ctx.nthreads = Rf_isNull(th) ? 1 : std::max(1, Rf_asInteger(th));
  SEXP st = opt_get(opts, "strict");
  ctx.strict = !Rf_isNull(st) && Rf_asLogical(st) == TRUE;
}

// Opens the descriptor used by the direct path unless disabled or impossible.
void open_direct(hid_t file, const std::string& path, SEXP opts, Ctx& ctx, FdGuard& fd) {
  SEXP dr = opt_get(opts, "direct");
  if (!(Rf_isNull(dr) || Rf_asLogical(dr) == TRUE)) return;
  hsize_t ub = 0;
  H5Handle fcpl(H5Fget_create_plist(file), H5Pclose, "H5Fget_create_plist");
  H5Pget_userblock(fcpl, &ub);
#ifndef _WIN32
  if (ub == 0) {
    fd.fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    ctx.fd = fd.fd;
  }
#endif
}

void attach_results(cpp11::sexp res, const Ctx& ctx) {
  std::vector<std::string> tp, te, tv, td, tx;
  for (const auto& r : ctx.types) {
    tp.push_back(r.path);
    te.push_back(r.encoding);
    tv.push_back(r.version);
    td.push_back(r.dtype);
    tx.push_back(r.extra);
  }
  std::vector<std::string> tnames{"path", "encoding", "version", "dtype", "extra"};
  std::vector<cpp11::sexp> tcols{string_vec(tp), string_vec(te), string_vec(tv), string_vec(td), string_vec(tx)};
  cpp11::sexp types = named_list(tnames, tcols);
  cpp11::sexp warns = string_vec(ctx.warnings);
  cpp11::sexp notes = string_vec(ctx.fallback_notes);
  cpp11::sexp stats(alloc_vector(REALSXP, 2));
  REAL(stats)[0] = static_cast<double>(ctx.direct_reads);
  REAL(stats)[1] = static_cast<double>(ctx.fallback_reads);
  cpp11::unwind_protect([&] {
    Rf_setAttrib(res, Rf_install("annr_types"), types);
    Rf_setAttrib(res, Rf_install("annr_warnings"), warns);
    Rf_setAttrib(res, Rf_install("annr_fallbacks"), notes);
    Rf_setAttrib(res, Rf_install("annr_stats"), stats);
  });
}

}  // namespace

// Entry point used by api.cpp.
SEXP read_h5ad_impl(const std::string& path, SEXP opts) {
  Ctx ctx;
  configure_ctx(ctx, path, opts);
  auto get = [&](const char* nm) { return opt_get(opts, nm); };
  SEXP xc = get("x_colmajor");
  const bool x_colmajor = Rf_isNull(xc) || Rf_asLogical(xc) == TRUE;
  std::vector<SlotSel> sels(NSLOTS);
  SEXP slots = get("slots");
  if (!Rf_isNull(slots)) {
    SEXP sn = Rf_getAttrib(slots, R_NamesSymbol);
    for (R_xlen_t i = 0; i < Rf_xlength(slots); ++i) {
      std::string nm = CHAR(STRING_ELT(sn, i));
      for (int s = 0; s < NSLOTS; ++s) {
        if (nm == SLOTS[s]) sels[static_cast<size_t>(s)] = parse_sel(VECTOR_ELT(slots, i));
      }
    }
  }

  H5Handle file = open_file_ro(path);
  FdGuard fd;
  open_direct(file, path, opts, ctx, fd);
  cpp11::sexp res = read_anndata(file, "", ctx, x_colmajor, &sels, false);
  attach_results(res, ctx);
  return res;
}

// Reads one element (e.g. "obsm/X_umap") of a file.
SEXP read_elem_impl(const std::string& path, const std::string& elem, SEXP opts) {
  Ctx ctx;
  configure_ctx(ctx, path, opts);
  H5Handle file = open_file_ro(path);
  FdGuard fd;
  open_direct(file, path, opts, ctx, fd);
  if (!link_exists(file, elem)) throw std::runtime_error("element '" + elem + "' not found in '" + path + "'");
  cpp11::sexp res = read_elem(file, elem, elem, ctx, ElemOpts{});
  if (Rf_isNull(res)) return R_NilValue;
  attach_results(res, ctx);
  return res;
}

}  // namespace annr
