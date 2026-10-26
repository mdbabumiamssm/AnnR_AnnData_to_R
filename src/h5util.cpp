// HDF5 helpers: error reporting, handles, type classification, attributes.
#include "common.h"

namespace annr {

// ---------------------------------------------------------------------------
// Errors

static herr_t collect_error(unsigned n, const H5E_error2_t* err, void* data) {
  auto* out = static_cast<std::vector<std::string>*>(data);
  if (err != nullptr && err->desc != nullptr && err->desc[0] != '\0' && n < 6) {
    out->push_back(err->desc);
  }
  return 0;
}

std::string h5_error_stack() {
  std::vector<std::string> parts;
  H5Ewalk2(H5E_DEFAULT, H5E_WALK_UPWARD, collect_error, &parts);
  H5Eclear2(H5E_DEFAULT);
  std::string s;
  for (const auto& p : parts) {
    if (!s.empty()) s += "; ";
    s += p;
  }
  return s;
}

void h5_fail(const std::string& what) {
  std::string stack = h5_error_stack();
  if (stack.empty()) throw std::runtime_error(what);
  throw std::runtime_error(what + " [HDF5: " + stack + "]");
}

// ---------------------------------------------------------------------------
// Handles

H5Handle open_object(hid_t loc, const std::string& name) {
  return H5Handle(H5Oopen(loc, name.c_str(), H5P_DEFAULT), H5Oclose,
                  ("cannot open '" + name + "'").c_str());
}

H5Handle open_dataset(hid_t loc, const std::string& name) {
  return H5Handle(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose,
                  ("cannot open dataset '" + name + "'").c_str());
}

H5Handle open_group(hid_t loc, const std::string& name) {
  return H5Handle(H5Gopen2(loc, name.c_str(), H5P_DEFAULT), H5Gclose,
                  ("cannot open group '" + name + "'").c_str());
}

H5Handle dataset_type(hid_t dset) {
  return H5Handle(H5Dget_type(dset), H5Tclose, "H5Dget_type failed");
}

H5Handle dataset_space(hid_t dset) {
  return H5Handle(H5Dget_space(dset), H5Sclose, "H5Dget_space failed");
}

bool link_exists(hid_t loc, const std::string& name) {
  htri_t r = H5Lexists(loc, name.c_str(), H5P_DEFAULT);
  if (r < 0) {
    H5Eclear2(H5E_DEFAULT);
    return false;
  }
  if (r == 0) return false;
  // Dangling soft/external links exist but cannot be opened.
  htri_t ok = H5Oexists_by_name(loc, name.c_str(), H5P_DEFAULT);
  if (ok < 0) H5Eclear2(H5E_DEFAULT);
  return ok > 0;
}

H5I_type_t object_kind(hid_t obj) { return H5Iget_type(obj); }

static herr_t collect_link(hid_t, const char* name, const H5L_info_t*, void* data) {
  static_cast<std::vector<std::string>*>(data)->emplace_back(name);
  return 0;
}

std::vector<std::string> group_children(hid_t group) {
  H5_index_t idx = H5_INDEX_NAME;
  hid_t gcpl = H5Gget_create_plist(group);
  if (gcpl >= 0) {
    unsigned flags = 0;
    if (H5Pget_link_creation_order(gcpl, &flags) >= 0 && (flags & H5P_CRT_ORDER_TRACKED)) {
      idx = H5_INDEX_CRT_ORDER;
    }
    H5Pclose(gcpl);
  }
  std::vector<std::string> names;
  hsize_t pos = 0;
  if (H5Literate(group, idx, H5_ITER_INC, &pos, collect_link, &names) < 0) {
    if (idx == H5_INDEX_CRT_ORDER) {
      H5Eclear2(H5E_DEFAULT);
      names.clear();
      pos = 0;
      if (H5Literate(group, H5_INDEX_NAME, H5_ITER_INC, &pos, collect_link, &names) < 0) {
        h5_fail("cannot list group members");
      }
    } else {
      h5_fail("cannot list group members");
    }
  }
  return names;
}

std::vector<hsize_t> space_dims(hid_t space) {
  int rank = H5Sget_simple_extent_ndims(space);
  if (rank < 0) h5_fail("cannot query dataspace rank");
  std::vector<hsize_t> dims(static_cast<size_t>(rank));
  if (rank > 0 && H5Sget_simple_extent_dims(space, dims.data(), nullptr) < 0) {
    h5_fail("cannot query dataspace dims");
  }
  return dims;
}

bool space_is_null(hid_t space) { return H5Sget_simple_extent_type(space) == H5S_NULL; }

// ---------------------------------------------------------------------------
// Types

size_t prim_size(Prim p) {
  switch (p) {
    case Prim::I8: case Prim::U8: case Prim::BOOL: return 1;
    case Prim::I16: case Prim::U16: case Prim::F16: return 2;
    case Prim::I32: case Prim::U32: case Prim::F32: return 4;
    case Prim::I64: case Prim::U64: case Prim::F64: return 8;
    default: return 0;
  }
}

bool prim_is_float(Prim p) { return p == Prim::F16 || p == Prim::F32 || p == Prim::F64; }

static bool is_h5py_bool(hid_t type) {
  if (H5Tget_nmembers(type) != 2) return false;
  hid_t super = H5Tget_super(type);
  if (super < 0) {
    H5Eclear2(H5E_DEFAULT);
    return false;
  }
  size_t ssize = H5Tget_size(super);
  H5Tclose(super);
  if (ssize != 1) return false;
  bool ok = true;
  for (unsigned i = 0; i < 2 && ok; ++i) {
    char* name = H5Tget_member_name(type, i);
    if (name == nullptr) {
      H5Eclear2(H5E_DEFAULT);
      return false;
    }
    int8_t value = 0;
    herr_t rc = H5Tget_member_value(type, i, &value);
    std::string nm(name);
    H5free_memory(name);
    if (rc < 0) ok = false;
    else if (nm == "FALSE") ok = (value == 0);
    else if (nm == "TRUE") ok = (value == 1);
    else ok = false;
  }
  return ok;
}

TypeDesc describe_type(hid_t type) {
  TypeDesc d;
  d.cls = H5Tget_class(type);
  d.size = H5Tget_size(type);
  switch (d.cls) {
    case H5T_INTEGER: {
      d.is_numeric = true;
      d.is_signed = H5Tget_sign(type) == H5T_SGN_2;
      size_t bits = d.size * 8;
      d.dtype = std::string(d.is_signed ? "int" : "uint") + std::to_string(bits);
      bool le = H5Tget_order(type) == H5T_ORDER_LE;
      if (H5Tget_order(type) == H5T_ORDER_BE && d.size > 1) d.dtype = ">" + d.dtype;
      bool std_layout = H5Tget_precision(type) == bits && H5Tget_offset(type) == 0;
      if (le && std_layout) {
        switch (d.size) {
          case 1: d.prim = d.is_signed ? Prim::I8 : Prim::U8; break;
          case 2: d.prim = d.is_signed ? Prim::I16 : Prim::U16; break;
          case 4: d.prim = d.is_signed ? Prim::I32 : Prim::U32; break;
          case 8: d.prim = d.is_signed ? Prim::I64 : Prim::U64; break;
          default: break;
        }
      }
      break;
    }
    case H5T_FLOAT: {
      d.is_numeric = true;
      d.dtype = "float" + std::to_string(d.size * 8);
      if (H5Tget_order(type) == H5T_ORDER_BE) d.dtype = ">" + d.dtype;
      if (H5Tequal(type, H5T_IEEE_F32LE) > 0) {
        d.prim = Prim::F32;
      } else if (H5Tequal(type, H5T_IEEE_F64LE) > 0) {
        d.prim = Prim::F64;
      } else if (d.size == 2 && H5Tget_order(type) == H5T_ORDER_LE) {
        size_t spos, epos, esize, mpos, msize;
        if (H5Tget_fields(type, &spos, &epos, &esize, &mpos, &msize) >= 0 && spos == 15 &&
            epos == 10 && esize == 5 && mpos == 0 && msize == 10 && H5Tget_ebias(type) == 15) {
          d.prim = Prim::F16;
        }
      }
      break;
    }
    case H5T_ENUM: {
      if (is_h5py_bool(type)) {
        d.is_numeric = true;
        d.is_bool = true;
        d.prim = Prim::BOOL;
        d.dtype = "bool";
      } else {
        d.dtype = "enum";
      }
      break;
    }
    case H5T_STRING: {
      d.is_string = true;
      htri_t v = H5Tis_variable_str(type);
      d.is_vlen_string = v > 0;
      d.cset = H5Tget_cset(type);
      d.strpad = H5Tget_strpad(type);
      d.dtype = d.is_vlen_string ? "str" : ("S" + std::to_string(d.size));
      break;
    }
    case H5T_COMPOUND: {
      d.is_compound = true;
      d.dtype = "compound";
      if (H5Tget_nmembers(type) == 2) {
        char* n0 = H5Tget_member_name(type, 0);
        char* n1 = H5Tget_member_name(type, 1);
        bool names_ok = n0 && n1 && std::string(n0) == "r" && std::string(n1) == "i";
        if (n0) H5free_memory(n0);
        if (n1) H5free_memory(n1);
        if (names_ok && H5Tget_member_class(type, 0) == H5T_FLOAT &&
            H5Tget_member_class(type, 1) == H5T_FLOAT) {
          d.is_complex = true;
          d.dtype = d.size == 8 ? "complex64" : "complex128";
        }
      }
      break;
    }
    case H5T_REFERENCE: d.dtype = "reference"; break;
    case H5T_OPAQUE: d.dtype = "opaque"; break;
    case H5T_ARRAY: d.dtype = "array"; break;
    case H5T_VLEN: d.dtype = "vlen"; break;
    case H5T_BITFIELD: d.dtype = "bitfield"; break;
    default: d.dtype = "unknown"; break;
  }
  H5Eclear2(H5E_DEFAULT);
  return d;
}

// ---------------------------------------------------------------------------
// Attributes

bool has_attr(hid_t obj, const char* name) {
  htri_t r = H5Aexists(obj, name);
  if (r < 0) {
    H5Eclear2(H5E_DEFAULT);
    return false;
  }
  return r > 0;
}

namespace {

struct VlenStrings {
  hid_t memtype;
  hid_t space;
  std::vector<char*> ptrs;
  ~VlenStrings() {
    if (!ptrs.empty()) {
#if H5_VERSION_GE(1, 12, 0)
      H5Treclaim(memtype, space, H5P_DEFAULT, ptrs.data());
#else
      H5Dvlen_reclaim(memtype, space, H5P_DEFAULT, ptrs.data());
#endif
    }
  }
};

std::vector<std::string> read_string_attr(hid_t attr, hid_t ftype, hid_t space) {
  std::vector<std::string> out;
  hssize_t n = H5Sget_simple_extent_npoints(space);
  if (n <= 0 || space_is_null(space)) return out;
  if (H5Tis_variable_str(ftype) > 0) {
    H5Handle mt(H5Tcopy(H5T_C_S1), H5Tclose, "H5Tcopy");
    H5Tset_size(mt, H5T_VARIABLE);
    H5Tset_cset(mt, H5Tget_cset(ftype));
    VlenStrings holder{mt, space, std::vector<char*>(static_cast<size_t>(n), nullptr)};
    if (H5Aread(attr, mt, holder.ptrs.data()) < 0) h5_fail("cannot read string attribute");
    for (char* p : holder.ptrs) out.emplace_back(p ? p : "");
  } else {
    size_t sz = H5Tget_size(ftype);
    H5Handle mt(H5Tcopy(ftype), H5Tclose, "H5Tcopy");
    std::vector<char> buf(sz * static_cast<size_t>(n) + 1, '\0');
    if (H5Aread(attr, mt, buf.data()) < 0) h5_fail("cannot read string attribute");
    for (hssize_t i = 0; i < n; ++i) {
      const char* s = buf.data() + i * sz;
      size_t len = strnlen(s, sz);
      out.emplace_back(s, len);
    }
  }
  return out;
}

}  // namespace

std::vector<std::string> attr_strings(hid_t obj, const char* name) {
  std::vector<std::string> out;
  if (!has_attr(obj, name)) return out;
  H5Handle attr(H5Aopen(obj, name, H5P_DEFAULT), H5Aclose, "cannot open attribute");
  H5Handle ftype(H5Aget_type(attr), H5Tclose, "H5Aget_type");
  H5Handle space(H5Aget_space(attr), H5Sclose, "H5Aget_space");
  if (H5Tget_class(ftype) != H5T_STRING) return out;  // e.g. empty float64 array
  return read_string_attr(attr, ftype, space);
}

std::string attr_string(hid_t obj, const char* name, const std::string& dflt) {
  if (!has_attr(obj, name)) return dflt;
  H5Handle attr(H5Aopen(obj, name, H5P_DEFAULT), H5Aclose, "cannot open attribute");
  H5Handle ftype(H5Aget_type(attr), H5Tclose, "H5Aget_type");
  H5Handle space(H5Aget_space(attr), H5Sclose, "H5Aget_space");
  if (H5Tget_class(ftype) != H5T_STRING) return dflt;
  std::vector<std::string> v = read_string_attr(attr, ftype, space);
  return v.empty() ? dflt : v[0];
}

std::vector<int64_t> attr_int64s(hid_t obj, const char* name) {
  std::vector<int64_t> out;
  if (!has_attr(obj, name)) return out;
  H5Handle attr(H5Aopen(obj, name, H5P_DEFAULT), H5Aclose, "cannot open attribute");
  H5Handle ftype(H5Aget_type(attr), H5Tclose, "H5Aget_type");
  H5Handle space(H5Aget_space(attr), H5Sclose, "H5Aget_space");
  H5T_class_t cls = H5Tget_class(ftype);
  if (cls != H5T_INTEGER && cls != H5T_FLOAT) return out;
  hssize_t n = H5Sget_simple_extent_npoints(space);
  if (n <= 0 || space_is_null(space)) return out;
  if (cls == H5T_INTEGER) {
    out.resize(static_cast<size_t>(n));
    if (H5Aread(attr, H5T_NATIVE_INT64, out.data()) < 0) h5_fail("cannot read integer attribute");
  } else {
    std::vector<double> tmp(static_cast<size_t>(n));
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, tmp.data()) < 0) h5_fail("cannot read attribute");
    for (double v : tmp) out.push_back(static_cast<int64_t>(v));
  }
  return out;
}

bool attr_bool(hid_t obj, const char* name, bool dflt) {
  if (!has_attr(obj, name)) return dflt;
  H5Handle attr(H5Aopen(obj, name, H5P_DEFAULT), H5Aclose, "cannot open attribute");
  H5Handle ftype(H5Aget_type(attr), H5Tclose, "H5Aget_type");
  H5Handle space(H5Aget_space(attr), H5Sclose, "H5Aget_space");
  hssize_t n = H5Sget_simple_extent_npoints(space);
  if (n <= 0 || space_is_null(space)) return dflt;
  H5T_class_t cls = H5Tget_class(ftype);
  if (cls == H5T_ENUM) {
    H5Handle mt(H5Tget_native_type(ftype, H5T_DIR_ASCEND), H5Tclose, "native enum");
    std::vector<uint8_t> buf(H5Tget_size(mt) * static_cast<size_t>(n));
    if (H5Aread(attr, mt, buf.data()) < 0) h5_fail("cannot read boolean attribute");
    bool any = false;
    for (size_t i = 0; i < H5Tget_size(mt); ++i) any = any || buf[i] != 0;
    return any;
  }
  if (cls == H5T_INTEGER) {
    std::vector<int64_t> buf(static_cast<size_t>(n));
    if (H5Aread(attr, H5T_NATIVE_INT64, buf.data()) < 0) h5_fail("cannot read attribute");
    return buf[0] != 0;
  }
  if (cls == H5T_FLOAT) {
    std::vector<double> buf(static_cast<size_t>(n));
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, buf.data()) < 0) h5_fail("cannot read attribute");
    return buf[0] != 0;
  }
  if (cls == H5T_STRING) {
    std::string s = attr_string(obj, name, "");
    return s == "True" || s == "true" || s == "TRUE" || s == "1";
  }
  return dflt;
}

// ---------------------------------------------------------------------------
// R helpers

SEXP make_s4(const char* cls) {
  return cpp11::unwind_protect([&] {
    SEXP def = PROTECT(R_do_MAKE_CLASS(cls));
    SEXP obj = R_do_new_object(def);
    UNPROTECT(1);
    return obj;
  });
}

void set_dim(SEXP x, const std::vector<hsize_t>& dims) {
  for (hsize_t d : dims) {
    if (d > static_cast<hsize_t>(INT_MAX)) throw std::runtime_error("dimension exceeds R's integer range");
  }
  cpp11::unwind_protect([&] {
    SEXP d = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(dims.size())));
    for (size_t i = 0; i < dims.size(); ++i) INTEGER(d)[i] = static_cast<int>(dims[i]);
    Rf_setAttrib(x, R_DimSymbol, d);
    UNPROTECT(1);
  });
}

std::string join_path(const std::string& parent, const std::string& name) {
  if (parent.empty() || parent == "/") return name;
  return parent + "/" + name;
}

}  // namespace annr
