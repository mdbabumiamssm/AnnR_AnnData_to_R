// Writer: creates anndata-compatible HDF5 objects from R data.
//
// Values are converted to the target dtype and compressed by worker threads,
// then handed to HDF5 chunk by chunk (H5Dwrite_chunk), or written block-wise
// for contiguous datasets.
#include "common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace annr {

namespace {

// ---------------------------------------------------------------------------
// Output types

enum class Out { F64, F32, F16, I8, I16, I32, I64, U8, U16, U32, U64, BOOL };

size_t out_size(Out o) {
  switch (o) {
    case Out::I8: case Out::U8: case Out::BOOL: return 1;
    case Out::I16: case Out::U16: case Out::F16: return 2;
    case Out::I32: case Out::U32: case Out::F32: return 4;
    default: return 8;
  }
}

const char* out_name(Out o) {
  switch (o) {
    case Out::F64: return "float64";
    case Out::F32: return "float32";
    case Out::F16: return "float16";
    case Out::I8: return "int8";
    case Out::I16: return "int16";
    case Out::I32: return "int32";
    case Out::I64: return "int64";
    case Out::U8: return "uint8";
    case Out::U16: return "uint16";
    case Out::U32: return "uint32";
    case Out::U64: return "uint64";
    case Out::BOOL: return "bool";
  }
  return "?";
}

bool parse_out(const std::string& s_in, Out& o) {
  const std::string s = (!s_in.empty() && s_in[0] == '>') ? s_in.substr(1) : s_in;
  static const std::pair<const char*, Out> table[] = {
      {"float64", Out::F64}, {"float32", Out::F32}, {"float16", Out::F16}, {"int8", Out::I8},
      {"int16", Out::I16},   {"int32", Out::I32},   {"int64", Out::I64},   {"uint8", Out::U8},
      {"uint16", Out::U16},  {"uint32", Out::U32},  {"uint64", Out::U64},  {"bool", Out::BOOL}};
  for (const auto& e : table) {
    if (s == e.first) {
      o = e.second;
      return true;
    }
  }
  return false;
}

void out_range(Out o, double& lo, double& hi) {
  switch (o) {
    case Out::I8: lo = -128; hi = 127; break;
    case Out::I16: lo = -32768; hi = 32767; break;
    case Out::I32: lo = -2147483648.0; hi = 2147483647.0; break;
    case Out::I64: lo = -9223372036854775808.0; hi = 9223372036854774784.0; break;
    case Out::U8: lo = 0; hi = 255; break;
    case Out::U16: lo = 0; hi = 65535; break;
    case Out::U32: lo = 0; hi = 4294967295.0; break;
    case Out::U64: lo = 0; hi = 18446744073709549568.0; break;
    default: lo = -INFINITY; hi = INFINITY; break;
  }
}

H5Handle bool_type() {
  H5Handle t(H5Tenum_create(H5T_NATIVE_INT8), H5Tclose, "H5Tenum_create");
  int8_t f = 0, v = 1;
  H5Tenum_insert(t, "FALSE", &f);
  H5Tenum_insert(t, "TRUE", &v);
  return t;
}

bool hint_big_endian(const std::string& hint) { return !hint.empty() && hint[0] == '>'; }

void byteswap(uint8_t* buf, size_t n, size_t es) {
  if (es <= 1) return;
  for (size_t i = 0; i < n; ++i) std::reverse(buf + i * es, buf + (i + 1) * es);
}

H5Handle file_type(Out o, bool big_endian = false) {
  if (big_endian && o != Out::BOOL) {
    H5Handle t = file_type(o, false);
    if (H5Tset_order(t, H5T_ORDER_BE) < 0) h5_fail("H5Tset_order");
    return t;
  }
  hid_t base;
  switch (o) {
    case Out::F64: base = H5T_IEEE_F64LE; break;
    case Out::F32: base = H5T_IEEE_F32LE; break;
#ifdef H5T_IEEE_F16LE
    case Out::F16: base = H5T_IEEE_F16LE; break;
#else
    case Out::F16: {
      H5Handle t(H5Tcopy(H5T_IEEE_F32LE), H5Tclose, "H5Tcopy");
      H5Tset_fields(t, 15, 10, 5, 0, 10);
      H5Tset_size(t, 2);
      H5Tset_ebias(t, 15);
      return t;
    }
#endif
    case Out::I8: base = H5T_STD_I8LE; break;
    case Out::I16: base = H5T_STD_I16LE; break;
    case Out::I32: base = H5T_STD_I32LE; break;
    case Out::I64: base = H5T_STD_I64LE; break;
    case Out::U8: base = H5T_STD_U8LE; break;
    case Out::U16: base = H5T_STD_U16LE; break;
    case Out::U32: base = H5T_STD_U32LE; break;
    case Out::U64: base = H5T_STD_U64LE; break;
    case Out::BOOL: return bool_type();
  }
  return H5Handle(H5Tcopy(base), H5Tclose, "H5Tcopy");
}

H5Handle utf8_type() {
  H5Handle t(H5Tcopy(H5T_C_S1), H5Tclose, "H5Tcopy");
  H5Tset_size(t, H5T_VARIABLE);
  H5Tset_cset(t, H5T_CSET_UTF8);
  return t;
}

uint16_t double_to_half(double d) {
  float f = static_cast<float>(d);
  uint32_t x;
  std::memcpy(&x, &f, 4);
  uint32_t sign = (x >> 16) & 0x8000u;
  int32_t exp = static_cast<int32_t>((x >> 23) & 0xff) - 127 + 15;
  uint32_t mant = x & 0x7fffffu;
  if (((x >> 23) & 0xff) == 0xff) return static_cast<uint16_t>(sign | 0x7c00u | (mant ? 0x200u : 0));
  if (exp >= 31) return static_cast<uint16_t>(sign | 0x7c00u);
  if (exp <= 0) {
    if (exp < -10) return static_cast<uint16_t>(sign);
    mant |= 0x800000u;
    uint32_t shift = static_cast<uint32_t>(14 - exp);
    uint32_t h = mant >> shift;
    if ((mant >> (shift - 1)) & 1u) h += 1;  // round half up (lossless check guards)
    return static_cast<uint16_t>(sign | h);
  }
  uint32_t h = sign | (static_cast<uint32_t>(exp) << 10) | (mant >> 13);
  if (mant & 0x1000u) h += 1;
  return static_cast<uint16_t>(h);
}

double half_to_double_w(uint16_t h) {
  uint32_t sign = (h >> 15) & 1u, exp = (h >> 10) & 0x1fu, mant = h & 0x3ffu;
  double v;
  if (exp == 0) v = std::ldexp(static_cast<double>(mant), -24);
  else if (exp == 31) v = mant == 0 ? INFINITY : NAN;
  else v = std::ldexp(static_cast<double>(mant | 0x400u), static_cast<int>(exp) - 25);
  return sign ? -v : v;
}

// ---------------------------------------------------------------------------
// Source data (R memory)

enum class Src { DBL, INT, LGL, I64, U64 };

struct Source {
  Src kind;
  const void* data;
  size_t n;
};

// Decides the on-disk type: the hint is honoured only when lossless.
Out resolve_out(const Source& s, const std::string& hint, int nthreads) {
  Out want;
  bool has_hint = parse_out(hint, want);
  if (s.kind == Src::LGL) return Out::BOOL;
  if (s.kind == Src::I64) return Out::I64;
  if (s.kind == Src::U64) return Out::U64;
  const size_t n = s.n;
  const size_t nb = block_count(n, 1 << 16, nthreads);
  if (s.kind == Src::INT) {
    const int* v = static_cast<const int*>(s.data);
    if (!has_hint) return Out::I32;
    if (want == Out::F64) return Out::F64;
    std::atomic<bool> ok{true};
    parallel_for(nb, nthreads, [&](size_t b, int) {
      size_t lo = n * b / nb, hi = n * (b + 1) / nb;
      double rlo, rhi;
      out_range(want, rlo, rhi);
      for (size_t i = lo; i < hi && ok; ++i) {
        int x = v[i];
        if (want == Out::BOOL) {
          if (x != 0 && x != 1) ok = false;
        } else if (want == Out::F32 || want == Out::F16) {
          double lim = want == Out::F32 ? 16777216.0 : 2048.0;
          if (std::fabs(static_cast<double>(x)) > lim) ok = false;
        } else if (x < rlo || x > rhi) {
          ok = false;
        }
      }
    });
    return ok ? want : Out::I32;
  }
  // doubles
  const double* v = static_cast<const double*>(s.data);
  if (!has_hint || want == Out::F64) return Out::F64;
  std::atomic<bool> ok{true};
  parallel_for(nb, nthreads, [&](size_t b, int) {
    size_t lo = n * b / nb, hi = n * (b + 1) / nb;
    double rlo, rhi;
    out_range(want, rlo, rhi);
    for (size_t i = lo; i < hi && ok; ++i) {
      double x = v[i];
      if (want == Out::F32) {
        if (!std::isnan(x) && static_cast<double>(static_cast<float>(x)) != x) ok = false;
      } else if (want == Out::F16) {
        if (!std::isnan(x) && half_to_double_w(double_to_half(x)) != x) ok = false;
      } else if (want == Out::BOOL) {
        if (x != 0.0 && x != 1.0) ok = false;
      } else {
        if (!(x >= rlo && x <= rhi) || std::trunc(x) != x) ok = false;
      }
    }
  });
  return ok ? want : Out::F64;
}

template <typename T>
inline void store(uint8_t* dst, T v) {
  std::memcpy(dst, &v, sizeof(T));
}

// Converts element i of the source to the output type at dst.
inline void put(const Source& s, size_t i, Out o, uint8_t* dst) {
  if (s.kind == Src::I64 || s.kind == Src::U64) {  // exact 64-bit patterns
    std::memcpy(dst, static_cast<const int64_t*>(s.data) + i, 8);
    return;
  }
  double d;
  int64_t k;
  if (s.kind == Src::DBL) {
    d = static_cast<const double*>(s.data)[i];
    k = 0;
  } else {
    k = static_cast<const int*>(s.data)[i];
    d = static_cast<double>(k);
  }
  const bool from_int = s.kind != Src::DBL;
  switch (o) {
    case Out::F64: store<double>(dst, d); break;
    case Out::F32: store<float>(dst, static_cast<float>(d)); break;
    case Out::F16: store<uint16_t>(dst, double_to_half(d)); break;
    case Out::I8: store<int8_t>(dst, static_cast<int8_t>(from_int ? k : static_cast<int64_t>(d))); break;
    case Out::I16: store<int16_t>(dst, static_cast<int16_t>(from_int ? k : static_cast<int64_t>(d))); break;
    case Out::I32: store<int32_t>(dst, static_cast<int32_t>(from_int ? k : static_cast<int64_t>(d))); break;
    case Out::I64: store<int64_t>(dst, from_int ? k : static_cast<int64_t>(d)); break;
    case Out::U8: store<uint8_t>(dst, static_cast<uint8_t>(from_int ? k : static_cast<int64_t>(d))); break;
    case Out::U16: store<uint16_t>(dst, static_cast<uint16_t>(from_int ? k : static_cast<int64_t>(d))); break;
    case Out::U32: store<uint32_t>(dst, static_cast<uint32_t>(from_int ? k : static_cast<int64_t>(d))); break;
    case Out::U64:
      store<uint64_t>(dst, from_int ? static_cast<uint64_t>(k)
                                    : (d >= 9223372036854775808.0 ? static_cast<uint64_t>(d)
                                                                  : static_cast<uint64_t>(static_cast<int64_t>(d))));
      break;
    case Out::BOOL: store<int8_t>(dst, (from_int ? k != 0 : d != 0.0) ? 1 : 0); break;
  }
}

// Fills a C-order block (dims cd, origin org, valid extent ext) of the file
// array from R memory. colmajor: R dims == file shape (axis reversal);
// otherwise R memory is already in C order.
void fill_block(const Source& s, Out o, const std::vector<hsize_t>& shape, bool colmajor,
                const hsize_t* org, const hsize_t* ext, const hsize_t* cd, uint8_t* buf) {
  const size_t rank = shape.size();
  const size_t os = out_size(o);
  size_t total = 1;
  for (size_t k = 0; k < rank; ++k) total *= cd[k];
  bool full = true;
  for (size_t k = 0; k < rank; ++k) full = full && ext[k] == cd[k];
  if (!full) std::memset(buf, 0, total * os);
  if (rank <= 1) {
    const size_t n = rank == 0 ? 1 : ext[0];
    const size_t base = rank == 0 ? 0 : org[0];
    for (size_t i = 0; i < n; ++i) put(s, base + i, o, buf + i * os);
    return;
  }
  // C-order strides of the block and of the file array.
  std::vector<size_t> bstride(rank), cstride(rank), fstride(rank);
  bstride[rank - 1] = 1;
  cstride[rank - 1] = 1;
  for (size_t k = rank - 1; k-- > 0;) {
    bstride[k] = bstride[k + 1] * cd[k + 1];
    cstride[k] = cstride[k + 1] * shape[k + 1];
  }
  fstride[0] = 1;
  for (size_t k = 1; k < rank; ++k) fstride[k] = fstride[k - 1] * shape[k - 1];
  if (colmajor && rank == 2) {
    constexpr size_t T = 64;
    const size_t nrow = shape[0];
    for (size_t r0 = 0; r0 < ext[0]; r0 += T) {
      size_t rn = std::min(T, static_cast<size_t>(ext[0]) - r0);
      for (size_t c = 0; c < ext[1]; ++c) {
        size_t src = (org[1] + c) * nrow + org[0] + r0;
        uint8_t* dst = buf + ((r0)*cd[1] + c) * os;
        for (size_t r = 0; r < rn; ++r) put(s, src + r, o, dst + r * cd[1] * os);
      }
    }
    return;
  }
  std::vector<hsize_t> idx(rank, 0);
  for (;;) {
    size_t b = 0, src = 0;
    for (size_t k = 0; k < rank; ++k) {
      b += idx[k] * bstride[k];
      src += (org[k] + idx[k]) * (colmajor ? fstride[k] : cstride[k]);
    }
    if (!colmajor) {
      // the innermost run is contiguous in both
      for (size_t i = 0; i < ext[rank - 1]; ++i) put(s, src + i, o, buf + (b + i) * os);
      size_t k = rank - 1;
      bool done = true;
      while (k-- > 0) {
        if (++idx[k] < ext[k]) {
          done = false;
          break;
        }
        idx[k] = 0;
      }
      if (done) break;
      continue;
    }
    put(s, src, o, buf + b * os);
    size_t k = rank;
    bool done = true;
    while (k-- > 0) {
      if (++idx[k] < ext[k]) {
        done = false;
        break;
      }
      idx[k] = 0;
    }
    if (done) break;
  }
}

// ---------------------------------------------------------------------------
// Compression settings

struct Comp {
  int method = 0;  // 0 none, 1 gzip, 2 lzf, 3 zstd
  int level = 4;
  bool shuffle = false;
  size_t chunk_bytes = size_t(1) << 20;
};

Comp parse_comp(SEXP opts) {
  Comp c;
  if (Rf_isNull(opts)) return c;
  SEXP names = Rf_getAttrib(opts, R_NamesSymbol);
  for (R_xlen_t i = 0; i < Rf_xlength(opts); ++i) {
    std::string nm = CHAR(STRING_ELT(names, i));
    SEXP v = VECTOR_ELT(opts, i);
    if (nm == "compression") {
      std::string m = CHAR(STRING_ELT(v, 0));
      c.method = m == "gzip" ? 1 : m == "lzf" ? 2 : m == "zstd" ? 3 : 0;
    } else if (nm == "level") {
      c.level = Rf_asInteger(v);
    } else if (nm == "shuffle") {
      c.shuffle = Rf_asLogical(v) == TRUE;
    } else if (nm == "chunk_bytes") {
      c.chunk_bytes = static_cast<size_t>(Rf_asReal(v));
    }
  }
  if (c.method == 3 && !zstd_available()) throw std::runtime_error("AnnR was built without zstd support");
  return c;
}

void set_filters(hid_t dcpl, const Comp& c, size_t chunk_bytes) {
  if (c.shuffle && c.method != 0) H5Pset_shuffle(dcpl);
  if (c.method == 1) {
    if (H5Pset_deflate(dcpl, static_cast<unsigned>(std::clamp(c.level, 0, 9))) < 0) h5_fail("H5Pset_deflate");
  } else if (c.method == 2) {
    unsigned cd[3] = {4, 0x0105, static_cast<unsigned>(chunk_bytes)};
    if (H5Pset_filter(dcpl, FILTER_LZF, H5Z_FLAG_OPTIONAL, 3, cd) < 0) h5_fail("H5Pset_filter(lzf)");
  } else if (c.method == 3) {
    unsigned cd[1] = {static_cast<unsigned>(std::clamp(c.level, 1, 22))};
    if (H5Pset_filter(dcpl, FILTER_ZSTD, H5Z_FLAG_OPTIONAL, 1, cd) < 0) h5_fail("H5Pset_filter(zstd)");
  }
}

std::vector<hsize_t> choose_chunks(const std::vector<hsize_t>& shape, size_t es, size_t target) {
  const size_t rank = shape.size();
  std::vector<hsize_t> c(rank, 1);
  if (rank == 0) return c;
  size_t target_elems = std::max<size_t>(1, target / es);
  // Fill from the innermost dimension outwards.
  size_t remaining = target_elems;
  for (size_t k = rank; k-- > 0;) {
    hsize_t d = std::max<hsize_t>(1, shape[k]);
    hsize_t take = std::min<hsize_t>(d, std::max<size_t>(1, remaining));
    c[k] = take;
    remaining = std::max<size_t>(1, remaining / take);
  }
  return c;
}

// ---------------------------------------------------------------------------
// Writer handle

struct Writer {
  hid_t file = H5I_INVALID_HID;
  int nthreads = 1;
  ~Writer() {
    if (file >= 0) H5Fclose(file);
  }
};

Writer* get_writer(SEXP xp) {
  Writer* w = static_cast<Writer*>(R_ExternalPtrAddr(xp));
  if (w == nullptr || w->file < 0) throw std::runtime_error("the h5ad writer is closed");
  return w;
}

void writer_finalizer(SEXP xp) {
  Writer* w = static_cast<Writer*>(R_ExternalPtrAddr(xp));
  if (w != nullptr) {
    delete w;
    R_ClearExternalPtr(xp);
  }
}

H5Handle make_lcpl() {
  H5Handle lcpl(H5Pcreate(H5P_LINK_CREATE), H5Pclose, "H5Pcreate");
  H5Pset_create_intermediate_group(lcpl, 1);
  H5Pset_char_encoding(lcpl, H5T_CSET_UTF8);
  return lcpl;
}

H5Handle open_any(hid_t file, const std::string& path) {
  if (path.empty() || path == "/") return H5Handle(H5Gopen2(file, "/", H5P_DEFAULT), H5Gclose, "open root");
  return open_object(file, path);
}

// ---------------------------------------------------------------------------
// Attributes

void delete_attr_if_exists(hid_t obj, const char* name) {
  if (has_attr(obj, name)) H5Adelete(obj, name);
}

void attr_write_string(hid_t obj, const char* name, const std::string& value) {
  delete_attr_if_exists(obj, name);
  H5Handle t = utf8_type();
  H5Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "H5Screate");
  H5Handle a(H5Acreate2(obj, name, t, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, "H5Acreate2");
  const char* p = value.c_str();
  if (H5Awrite(a, t, &p) < 0) h5_fail("H5Awrite");
}

void attr_write_strings(hid_t obj, const char* name, const std::vector<std::string>& values) {
  delete_attr_if_exists(obj, name);
  hsize_t n = values.size();
  if (n == 0) {
    // h5py stores an empty list as an empty float64 array
    H5Handle sp(H5Screate_simple(1, &n, nullptr), H5Sclose, "H5Screate_simple");
    H5Handle a(H5Acreate2(obj, name, H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, "H5Acreate2");
    double dummy = 0;
    H5Awrite(a, H5T_NATIVE_DOUBLE, &dummy);
    return;
  }
  H5Handle t = utf8_type();
  H5Handle sp(H5Screate_simple(1, &n, nullptr), H5Sclose, "H5Screate_simple");
  H5Handle a(H5Acreate2(obj, name, t, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, "H5Acreate2");
  std::vector<const char*> ptrs;
  for (const auto& v : values) ptrs.push_back(v.c_str());
  if (H5Awrite(a, t, ptrs.data()) < 0) h5_fail("H5Awrite");
}

void attr_write_int64s(hid_t obj, const char* name, const std::vector<int64_t>& values, bool scalar) {
  delete_attr_if_exists(obj, name);
  hsize_t n = values.size();
  H5Handle sp(scalar ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &n, nullptr), H5Sclose, "H5Screate");
  H5Handle a(H5Acreate2(obj, name, H5T_STD_I64LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, "H5Acreate2");
  if (H5Awrite(a, H5T_NATIVE_INT64, values.data()) < 0) h5_fail("H5Awrite");
}

void attr_write_bool(hid_t obj, const char* name, bool value) {
  delete_attr_if_exists(obj, name);
  H5Handle t = bool_type();
  H5Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "H5Screate");
  H5Handle a(H5Acreate2(obj, name, t, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, "H5Acreate2");
  int8_t v = value ? 1 : 0;
  if (H5Awrite(a, t, &v) < 0) h5_fail("H5Awrite");
}

void set_encoding(hid_t obj, const std::string& enc, const std::string& ver) {
  if (enc.empty()) return;
  attr_write_string(obj, "encoding-type", enc);
  attr_write_string(obj, "encoding-version", ver);
}

std::vector<std::string> strings_from(SEXP x) {
  std::vector<std::string> out;
  if (TYPEOF(x) != STRSXP) return out;
  for (R_xlen_t i = 0; i < Rf_xlength(x); ++i) {
    SEXP s = STRING_ELT(x, i);
    out.emplace_back(s == NA_STRING ? "" : Rf_translateCharUTF8(s));
  }
  return out;
}

// ---------------------------------------------------------------------------
// Numeric datasets

Source source_of(SEXP x) {
  switch (TYPEOF(x)) {
    case REALSXP: return {Src::DBL, REAL(x), static_cast<size_t>(Rf_xlength(x))};
    case INTSXP: return {Src::INT, INTEGER(x), static_cast<size_t>(Rf_xlength(x))};
    case LGLSXP: return {Src::LGL, LOGICAL(x), static_cast<size_t>(Rf_xlength(x))};
    default: throw std::runtime_error("unsupported R type for a numeric dataset");
  }
}

// Uses the exact 64-bit patterns kept by the reader when the doubles still
// match them (i.e. the values were not modified in R).
void exact_int64_source(SEXP x, const std::string& hint, Source& src) {
  if (TYPEOF(x) != REALSXP) return;
  SEXP bits = Rf_getAttrib(x, Rf_install("h5ad_int64"));
  if (TYPEOF(bits) != REALSXP || Rf_xlength(bits) != Rf_xlength(x)) return;
  std::string h = (!hint.empty() && hint[0] == '>') ? hint.substr(1) : hint;
  const bool is_u = h == "uint64";
  if (!is_u && h != "int64") return;
  const int64_t* b = reinterpret_cast<const int64_t*>(REAL(bits));
  const double* v = REAL(x);
  for (R_xlen_t i = 0; i < Rf_xlength(x); ++i) {
    double expect = is_u ? static_cast<double>(static_cast<uint64_t>(b[i])) : static_cast<double>(b[i]);
    if (expect != v[i]) return;
  }
  src = {is_u ? Src::U64 : Src::I64, b, static_cast<size_t>(Rf_xlength(x))};
}

struct DsetOpts {
  bool unlimited = false;  // maxshape=(None,) like anndata's sparse components
  bool force_chunked = false;
};

// Creates and fills a numeric dataset. Returns the dtype written.
std::string write_numeric(Writer* w, const std::string& path, const Source& src,
                          const std::vector<hsize_t>& shape, bool colmajor, const std::string& hint,
                          const Comp& comp, const DsetOpts& dopt) {
  const Out o = resolve_out(src, hint, w->nthreads);
  const size_t os = out_size(o);
  const size_t rank = shape.size();
  // Big-endian output only when the source element was big-endian and the
  // hinted dtype survived the lossless check.
  Out hinted;
  const bool be = hint_big_endian(hint) && parse_out(hint, hinted) && hinted == o && o != Out::BOOL;
  size_t n = 1;
  for (hsize_t d : shape) n *= d;
  H5Handle ftype = file_type(o, be);
  std::vector<hsize_t> maxdims(shape);
  if (dopt.unlimited && rank == 1) maxdims[0] = H5S_UNLIMITED;
  H5Handle space(rank == 0 ? H5Screate(H5S_SCALAR)
                           : H5Screate_simple(static_cast<int>(rank), shape.data(), maxdims.data()),
                 H5Sclose, "H5Screate");
  H5Handle dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, "H5Pcreate");
  const bool chunked = rank > 0 && (comp.method != 0 || dopt.unlimited || dopt.force_chunked);
  std::vector<hsize_t> chunk;
  size_t chunk_elems = 0;
  if (chunked) {
    std::vector<hsize_t> base(shape);
    for (auto& d : base) d = std::max<hsize_t>(d, 1);
    chunk = choose_chunks(base, os, comp.chunk_bytes);
    if (n == 0 && rank == 1) chunk[0] = 1024;
    chunk_elems = 1;
    for (hsize_t c : chunk) chunk_elems *= c;
    if (H5Pset_chunk(dcpl, static_cast<int>(rank), chunk.data()) < 0) h5_fail("H5Pset_chunk");
    set_filters(dcpl, comp, chunk_elems * os);
  }
  H5Handle lcpl = make_lcpl();
  H5Handle dset(H5Dcreate2(w->file, path.c_str(), ftype, space, lcpl, dcpl, H5P_DEFAULT), H5Dclose,
                ("cannot create dataset '" + path + "'").c_str());
  if (n == 0) return out_name(o);

  const std::string dtype_name = (be ? ">" : "") + std::string(out_name(o));
  if (rank == 0) {
    std::vector<uint8_t> buf(os);
    put(src, 0, o, buf.data());
    H5Handle mt = file_type(o);
    if (H5Dwrite(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) < 0) h5_fail("H5Dwrite");
    return dtype_name;
  }
  const int nthreads = std::max(1, w->nthreads);

  if (!chunked) {
    // Contiguous: fill row blocks in parallel, write them as hyperslabs.
    size_t row_elems = 1;
    for (size_t k = 1; k < rank; ++k) row_elems *= shape[k];
    size_t rows_per = std::max<size_t>(1, (size_t(4) << 20) / std::max<size_t>(1, row_elems * os));
    const size_t nrows = shape[0];
    const size_t nblocks = (nrows + rows_per - 1) / rows_per;
    const size_t batch = static_cast<size_t>(nthreads) * 2;
    std::vector<std::vector<uint8_t>> bufs(batch);
    H5Handle mt = file_type(o);
    for (size_t b0 = 0; b0 < nblocks; b0 += batch) {
      size_t bn = std::min(batch, nblocks - b0);
      parallel_for(bn, nthreads, [&](size_t j, int) {
        size_t r0 = (b0 + j) * rows_per;
        size_t rn = std::min(rows_per, nrows - r0);
        std::vector<hsize_t> org(rank, 0), ext(shape);
        org[0] = r0;
        ext[0] = rn;
        bufs[j].resize(rn * row_elems * os);
        fill_block(src, o, shape, colmajor, org.data(), ext.data(), ext.data(), bufs[j].data());
      });
      for (size_t j = 0; j < bn; ++j) {
        size_t r0 = (b0 + j) * rows_per;
        size_t rn = std::min(rows_per, nrows - r0);
        std::vector<hsize_t> start(rank, 0), count(shape);
        start[0] = r0;
        count[0] = rn;
        H5Handle fs = dataset_space(dset);
        if (H5Sselect_hyperslab(fs, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr) < 0) {
          h5_fail("H5Sselect_hyperslab");
        }
        H5Handle ms(H5Screate_simple(static_cast<int>(rank), count.data(), nullptr), H5Sclose, "H5Screate_simple");
        if (H5Dwrite(dset, mt, ms, fs, H5P_DEFAULT, bufs[j].data()) < 0) h5_fail("H5Dwrite");
      }
    }
    return dtype_name;
  }

  // Chunked: fill + filter in parallel, then H5Dwrite_chunk in file order.
  std::vector<size_t> grid(rank);
  size_t nchunks = 1;
  for (size_t k = 0; k < rank; ++k) {
    grid[k] = (shape[k] + chunk[k] - 1) / chunk[k];
    nchunks *= grid[k];
  }
  const size_t raw_bytes = chunk_elems * os;
  size_t cap = raw_bytes;
  if (comp.method == 1) cap = deflate_bound(raw_bytes);
  if (comp.method == 3) cap = zstd_bound(raw_bytes);
  const unsigned comp_index = comp.shuffle && comp.method != 0 ? 1u : 0u;
  struct Out_ {
    std::vector<uint8_t> data;
    size_t size = 0;
    unsigned mask = 0;
    std::vector<hsize_t> org;
  };
  const size_t batch = static_cast<size_t>(nthreads) * 4;
  std::vector<Out_> outs(batch);
  std::vector<void*> states(static_cast<size_t>(nthreads), nullptr);
  struct StateGuard {
    std::vector<void*>& s;
    ~StateGuard() {
      for (void* p : s) free_deflate_state(p);
    }
  } guard{states};
  std::vector<std::vector<uint8_t>> scratch(static_cast<size_t>(nthreads)), scratch2(static_cast<size_t>(nthreads));
  for (size_t c0 = 0; c0 < nchunks; c0 += batch) {
    size_t cn = std::min(batch, nchunks - c0);
    parallel_for(cn, nthreads, [&](size_t j, int tid) {
      size_t lin = c0 + j;
      Out_& out = outs[j];
      out.org.assign(rank, 0);
      std::vector<hsize_t> ext(rank);
      size_t rem = lin;
      for (size_t k = rank; k-- > 0;) {
        out.org[k] = (rem % grid[k]) * chunk[k];
        rem /= grid[k];
      }
      for (size_t k = 0; k < rank; ++k) ext[k] = std::min<hsize_t>(chunk[k], shape[k] - out.org[k]);
      auto& raw = scratch[static_cast<size_t>(tid)];
      raw.resize(raw_bytes);
      fill_block(src, o, shape, colmajor, out.org.data(), ext.data(), chunk.data(), raw.data());
      if (be) byteswap(raw.data(), chunk_elems, os);
      const uint8_t* cur = raw.data();
      if (comp.method != 0 && comp.shuffle) {
        auto& sh = scratch2[static_cast<size_t>(tid)];
        sh.resize(raw_bytes);
        shuffle(raw.data(), sh.data(), raw_bytes, os);
        cur = sh.data();
      }
      out.mask = 0;
      if (comp.method == 0) {
        out.data.assign(cur, cur + raw_bytes);
        out.size = raw_bytes;
        return;
      }
      out.data.resize(cap);
      size_t got = 0;
      if (comp.method == 1) {
        got = deflate_zlib(cur, raw_bytes, out.data.data(), cap, std::clamp(comp.level, 0, 9),
                           &states[static_cast<size_t>(tid)]);
        if (got == 0) throw std::runtime_error("deflate failed");
      } else if (comp.method == 2) {
        got = raw_bytes > 4 ? lzf_compress(cur, raw_bytes, out.data.data(), raw_bytes - 1) : 0;
        if (got == 0) {  // incompressible: store unfiltered, flag the optional filter as skipped
          std::memcpy(out.data.data(), cur, raw_bytes);
          got = raw_bytes;
          out.mask = 1u << comp_index;
        }
      } else {
        got = zstd_compress(cur, raw_bytes, out.data.data(), cap, std::clamp(comp.level, 1, 22));
        if (got == 0) throw std::runtime_error("zstd compression failed");
      }
      out.size = got;
    });
    for (size_t j = 0; j < cn; ++j) {
      if (H5Dwrite_chunk(dset, H5P_DEFAULT, outs[j].mask, outs[j].org.data(), outs[j].size,
                         outs[j].data.data()) < 0) {
        h5_fail("H5Dwrite_chunk failed for '" + path + "'");
      }
    }
  }
  return dtype_name;
}

std::vector<hsize_t> shape_from(SEXP shape) {
  std::vector<hsize_t> s;
  for (R_xlen_t i = 0; i < Rf_xlength(shape); ++i) s.push_back(static_cast<hsize_t>(REAL(shape)[i]));
  return s;
}

// ---------------------------------------------------------------------------
// String datasets

// Returns the fixed width n when hint is "S<n>" and every string fits.
size_t fixed_width_for(SEXP x, const std::string& hint, bool& ascii) {
  ascii = true;
  if (hint.size() < 2 || hint[0] != 'S') return 0;
  size_t width = 0;
  try {
    width = static_cast<size_t>(std::stoul(hint.substr(1)));
  } catch (...) {
    return 0;
  }
  if (width == 0) return 0;
  bool fits = true;
  cpp11::unwind_protect([&] {
    for (R_xlen_t i = 0; i < Rf_xlength(x) && fits; ++i) {
      SEXP s = STRING_ELT(x, i);
      if (s == NA_STRING) {
        fits = false;
        break;
      }
      const char* c = Rf_translateCharUTF8(s);
      size_t len = std::strlen(c);
      for (size_t j = 0; j < len; ++j) {
        if (static_cast<unsigned char>(c[j]) & 0x80) ascii = false;
      }
      if (len > width) fits = false;
    }
  });
  return fits ? width : 0;
}

void write_strings_impl(Writer* w, const std::string& path, SEXP x, const std::vector<hsize_t>& shape,
                        bool colmajor, const Comp& comp, const std::string& hint) {
  const size_t rank = shape.size();
  size_t n = 1;
  for (hsize_t d : shape) n *= d;
  bool ascii = true;
  const size_t width = fixed_width_for(x, hint, ascii);
  H5Handle t = utf8_type();
  if (width > 0) {
    t = H5Handle(H5Tcopy(H5T_C_S1), H5Tclose, "H5Tcopy");
    H5Tset_size(t, width);
    H5Tset_strpad(t, H5T_STR_NULLPAD);
    H5Tset_cset(t, ascii ? H5T_CSET_ASCII : H5T_CSET_UTF8);
  }
  H5Handle space(rank == 0 ? H5Screate(H5S_SCALAR)
                           : H5Screate_simple(static_cast<int>(rank), shape.data(), nullptr),
                 H5Sclose, "H5Screate");
  H5Handle dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, "H5Pcreate");
  if (rank > 0 && comp.method != 0 && n > 0) {
    std::vector<hsize_t> base(shape);
    std::vector<hsize_t> chunk = choose_chunks(base, sizeof(char*) * 2, comp.chunk_bytes);
    H5Pset_chunk(dcpl, static_cast<int>(rank), chunk.data());
    size_t ce = 1;
    for (hsize_t c : chunk) ce *= c;
    set_filters(dcpl, comp, ce * 16);
  }
  H5Handle lcpl = make_lcpl();
  H5Handle dset(H5Dcreate2(w->file, path.c_str(), t, space, lcpl, dcpl, H5P_DEFAULT), H5Dclose,
                ("cannot create dataset '" + path + "'").c_str());
  if (n == 0) return;
  std::vector<const char*> ptrs(n);
  // Rf_translateCharUTF8 may allocate: done on the main thread, under protection.
  cpp11::unwind_protect([&] {
    for (size_t k = 0; k < n; ++k) {
      size_t src = k;
      if (colmajor && rank >= 2) {
        // C-order index k -> column-major position
        size_t rem = k, f = 0, stride = 1;
        std::vector<size_t> idx(rank);
        for (size_t j = rank; j-- > 0;) {
          idx[j] = rem % shape[j];
          rem /= shape[j];
        }
        for (size_t j = 0; j < rank; ++j) {
          f += idx[j] * stride;
          stride *= shape[j];
        }
        src = f;
      }
      SEXP s = STRING_ELT(x, static_cast<R_xlen_t>(src));
      ptrs[k] = s == NA_STRING ? "" : Rf_translateCharUTF8(s);
    }
  });
  if (width > 0) {
    std::vector<char> fixed(n * width, '\0');
    for (size_t k = 0; k < n; ++k) {
      size_t len = std::min(width, std::strlen(ptrs[k]));
      std::memcpy(fixed.data() + k * width, ptrs[k], len);
    }
    if (H5Dwrite(dset, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, fixed.data()) < 0) h5_fail("H5Dwrite strings");
    return;
  }
  if (H5Dwrite(dset, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data()) < 0) h5_fail("H5Dwrite strings");
}

}  // namespace

// ---------------------------------------------------------------------------
// R entry points

SEXP w_open(const std::string& path, int threads) {
  H5Handle fcpl(H5Pcreate(H5P_FILE_CREATE), H5Pclose, "H5Pcreate");
  hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, fcpl, H5P_DEFAULT);
  if (f < 0) h5_fail("cannot create '" + path + "'");
  Writer* w = new Writer();
  w->file = f;
  w->nthreads = std::max(1, threads);
  return cpp11::unwind_protect([&] {
    SEXP xp = PROTECT(R_MakeExternalPtr(w, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(xp, writer_finalizer, TRUE);
    UNPROTECT(1);
    return xp;
  });
}

void w_close(SEXP xp) {
  Writer* w = static_cast<Writer*>(R_ExternalPtrAddr(xp));
  if (w == nullptr) return;
  hid_t f = w->file;
  w->file = H5I_INVALID_HID;
  delete w;
  R_ClearExternalPtr(xp);
  if (f >= 0 && H5Fclose(f) < 0) h5_fail("error while closing the file");
}

void w_group(SEXP xp, const std::string& path, const std::string& enc, const std::string& ver) {
  Writer* w = get_writer(xp);
  if (path.empty() || path == "/") {
    H5Handle root(H5Gopen2(w->file, "/", H5P_DEFAULT), H5Gclose, "open root");
    set_encoding(root, enc, ver);
    return;
  }
  H5Handle lcpl = make_lcpl();
  H5Handle g(H5Gcreate2(w->file, path.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
             ("cannot create group '" + path + "'").c_str());
  set_encoding(g, enc, ver);
}

void w_attr(SEXP xp, const std::string& path, const std::string& name, SEXP value, const std::string& kind) {
  Writer* w = get_writer(xp);
  H5Handle obj = open_any(w->file, path);
  if (kind == "string") {
    std::vector<std::string> v = strings_from(value);
    attr_write_string(obj, name.c_str(), v.empty() ? "" : v[0]);
  } else if (kind == "strings") {
    attr_write_strings(obj, name.c_str(), strings_from(value));
  } else if (kind == "int64s" || kind == "int64") {
    if (TYPEOF(value) != REALSXP) throw std::runtime_error("integer attributes must be passed as doubles");
    std::vector<int64_t> v;
    for (R_xlen_t i = 0; i < Rf_xlength(value); ++i) v.push_back(static_cast<int64_t>(REAL(value)[i]));
    attr_write_int64s(obj, name.c_str(), v, kind == "int64");
  } else if (kind == "bool") {
    attr_write_bool(obj, name.c_str(), Rf_asLogical(value) == TRUE);
  } else {
    throw std::runtime_error("unknown attribute kind " + kind);
  }
}

std::string w_array(SEXP xp, const std::string& path, SEXP x, SEXP shape, bool colmajor,
                    const std::string& hint, const std::string& enc, const std::string& ver,
                    SEXP comp_opts, bool unlimited) {
  Writer* w = get_writer(xp);
  Comp comp = parse_comp(comp_opts);
  std::vector<hsize_t> sh = shape_from(shape);
  std::string dtype;
  if (TYPEOF(x) == STRSXP) {
    write_strings_impl(w, path, x, sh, colmajor, comp, hint);
    dtype = "str";
  } else if (TYPEOF(x) == CPLXSXP) {
    const size_t rank = sh.size();
    size_t n = 1;
    for (hsize_t d : sh) n *= d;
    H5Handle ct(H5Tcreate(H5T_COMPOUND, 16), H5Tclose, "H5Tcreate");
    H5Tinsert(ct, "r", 0, H5T_IEEE_F64LE);
    H5Tinsert(ct, "i", 8, H5T_IEEE_F64LE);
    H5Handle space(rank == 0 ? H5Screate(H5S_SCALAR) : H5Screate_simple(static_cast<int>(rank), sh.data(), nullptr),
                   H5Sclose, "H5Screate");
    H5Handle lcpl = make_lcpl();
    H5Handle dset(H5Dcreate2(w->file, path.c_str(), ct, space, lcpl, H5P_DEFAULT, H5P_DEFAULT), H5Dclose,
                  "cannot create complex dataset");
    std::vector<double> buf(2 * n);
    const Rcomplex* c = COMPLEX(x);
    for (size_t k = 0; k < n; ++k) {
      size_t src = k;
      if (colmajor && rank >= 2) {
        size_t rem = k, f = 0, stride = 1;
        std::vector<size_t> idx(rank);
        for (size_t j = rank; j-- > 0;) {
          idx[j] = rem % sh[j];
          rem /= sh[j];
        }
        for (size_t j = 0; j < rank; ++j) {
          f += idx[j] * stride;
          stride *= sh[j];
        }
        src = f;
      }
      buf[2 * k] = c[src].r;
      buf[2 * k + 1] = c[src].i;
    }
    H5Handle mt(H5Tcreate(H5T_COMPOUND, 16), H5Tclose, "H5Tcreate");
    H5Tinsert(mt, "r", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(mt, "i", 8, H5T_NATIVE_DOUBLE);
    if (n > 0 && H5Dwrite(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) < 0) h5_fail("H5Dwrite complex");
    dtype = "complex128";
  } else {
    DsetOpts dopt;
    dopt.unlimited = unlimited;
    Source src = source_of(x);
    exact_int64_source(x, hint, src);
    dtype = write_numeric(w, path, src, sh, colmajor, hint, comp, dopt);
  }
  if (!enc.empty()) {
    H5Handle obj = open_object(w->file, path);
    set_encoding(obj, enc, ver);
  }
  return dtype;
}

void w_null(SEXP xp, const std::string& path) {
  Writer* w = get_writer(xp);
  H5Handle sp(H5Screate(H5S_NULL), H5Sclose, "H5Screate");
  H5Handle lcpl = make_lcpl();
  H5Handle d(H5Dcreate2(w->file, path.c_str(), H5T_IEEE_F32LE, sp, lcpl, H5P_DEFAULT, H5P_DEFAULT),
             H5Dclose, "cannot create null dataset");
  set_encoding(d, "null", "0.1.0");
}

SEXP w_sparse(SEXP xp, const std::string& path, const std::string& fmt, SEXP shape, SEXP p, SEXP idx,
              SEXP x, SEXP hints, SEXP comp_opts) {
  Writer* w = get_writer(xp);
  Comp comp = parse_comp(comp_opts);
  std::vector<std::string> h = strings_from(hints);
  while (h.size() < 3) h.emplace_back("");
  {
    H5Handle lcpl = make_lcpl();
    H5Handle g(H5Gcreate2(w->file, path.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
               ("cannot create group '" + path + "'").c_str());
    set_encoding(g, fmt == "csr" ? "csr_matrix" : "csc_matrix", "0.1.0");
    std::vector<int64_t> sh{static_cast<int64_t>(REAL(shape)[0]), static_cast<int64_t>(REAL(shape)[1])};
    attr_write_int64s(g, "shape", sh, false);
  }
  DsetOpts dopt;
  dopt.unlimited = true;
  std::vector<hsize_t> nnz{static_cast<hsize_t>(Rf_xlength(x))};
  std::string d1 = write_numeric(w, path + "/data", source_of(x), nnz, false, h[0], comp, dopt);
  std::vector<hsize_t> ni{static_cast<hsize_t>(Rf_xlength(idx))};
  std::string idx_hint = h[1].empty() ? "int32" : h[1];
  std::string d2 = write_numeric(w, path + "/indices", source_of(idx), ni, false, idx_hint, comp, dopt);
  std::vector<hsize_t> np{static_cast<hsize_t>(Rf_xlength(p))};
  std::string ptr_hint = h[2].empty() ? "int32" : h[2];
  std::string d3 = write_numeric(w, path + "/indptr", source_of(p), np, false, ptr_hint, comp, dopt);
  cpp11::writable::strings out({d1, d2, d3});
  return out;
}

// Writes a data frame-like list as a compound (rec-array) dataset.
void w_recarray(SEXP xp, const std::string& path, SEXP cols, SEXP names, SEXP hints) {
  Writer* w = get_writer(xp);
  const R_xlen_t ncol = Rf_xlength(cols);
  std::vector<std::string> nm = strings_from(names);
  std::vector<std::string> hs = strings_from(hints);
  hs.resize(static_cast<size_t>(ncol));
  size_t n = 0;
  if (ncol > 0) {
    SEXP c0 = VECTOR_ELT(cols, 0);
    SEXP d0 = Rf_getAttrib(c0, R_DimSymbol);
    n = Rf_isNull(d0) ? static_cast<size_t>(Rf_xlength(c0)) : static_cast<size_t>(INTEGER(d0)[0]);
  }
  struct Field {
    Out out = Out::F64;
    bool is_str = false;
    size_t k = 1;  // values per record (array fields)
    size_t offset = 0, size = 0;
    H5Handle ftype;
  };
  std::vector<Field> fields(static_cast<size_t>(ncol));
  size_t rec_size = 0;
  for (R_xlen_t j = 0; j < ncol; ++j) {
    SEXP c = VECTOR_ELT(cols, j);
    Field& f = fields[static_cast<size_t>(j)];
    SEXP dim = Rf_getAttrib(c, R_DimSymbol);
    f.k = Rf_isNull(dim) ? 1 : static_cast<size_t>(INTEGER(dim)[1]);
    H5Handle base;
    if (TYPEOF(c) == STRSXP) {
      f.is_str = true;
      base = utf8_type();
    } else {
      f.out = resolve_out(source_of(c), hs[static_cast<size_t>(j)], w->nthreads);
      base = file_type(f.out);
    }
    if (f.k > 1 || !Rf_isNull(dim)) {
      hsize_t ad = f.k;
      f.ftype = H5Handle(H5Tarray_create2(base, 1, &ad), H5Tclose, "H5Tarray_create2");
    } else {
      f.ftype = H5Handle(H5Tcopy(base), H5Tclose, "H5Tcopy");
    }
    f.size = H5Tget_size(f.ftype);
    f.offset = rec_size;
    rec_size += f.size;
  }
  H5Handle ct(H5Tcreate(H5T_COMPOUND, std::max<size_t>(rec_size, 1)), H5Tclose, "H5Tcreate");
  for (R_xlen_t j = 0; j < ncol; ++j) {
    const Field& f = fields[static_cast<size_t>(j)];
    if (H5Tinsert(ct, nm[static_cast<size_t>(j)].c_str(), f.offset, f.ftype) < 0) h5_fail("H5Tinsert");
  }
  std::vector<uint8_t> buf(std::max<size_t>(1, rec_size * n), 0);
  std::vector<std::vector<std::string>> keep;
  for (R_xlen_t j = 0; j < ncol; ++j) {
    SEXP c = VECTOR_ELT(cols, j);
    const Field& f = fields[static_cast<size_t>(j)];
    if (f.is_str) {
      std::vector<std::string> sv = strings_from(c);
      keep.push_back(std::move(sv));
      const auto& ref = keep.back();
      for (size_t r = 0; r < n; ++r)
        for (size_t q = 0; q < f.k; ++q) {
          const char* p = ref[q * n + r].c_str();
          std::memcpy(buf.data() + r * rec_size + f.offset + q * sizeof(char*), &p, sizeof(char*));
        }
    } else {
      Source s = source_of(c);
      size_t es = out_size(f.out);
      for (size_t r = 0; r < n; ++r)
        for (size_t q = 0; q < f.k; ++q) put(s, q * n + r, f.out, buf.data() + r * rec_size + f.offset + q * es);
    }
  }
  hsize_t dims = n;
  H5Handle space(H5Screate_simple(1, &dims, nullptr), H5Sclose, "H5Screate_simple");
  H5Handle lcpl = make_lcpl();
  H5Handle dset(H5Dcreate2(w->file, path.c_str(), ct, space, lcpl, H5P_DEFAULT, H5P_DEFAULT), H5Dclose,
                ("cannot create dataset '" + path + "'").c_str());
  if (n > 0 && H5Dwrite(dset, ct, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()) < 0) h5_fail("H5Dwrite rec-array");
  set_encoding(dset, "rec-array", "0.2.0");
}

}  // namespace annr
