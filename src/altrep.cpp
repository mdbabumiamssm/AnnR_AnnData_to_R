// Large numeric vectors backed by malloc'd memory (ALTREP).
//
// R forces a full garbage collection whenever a large allocation exceeds the
// current vector heap, and grows the heap only marginally, so every big
// array of a read paid a full GC (most of the read time once Bioconductor
// packages are loaded). These vectors live outside R's accounted heap, behave
// like ordinary double/integer/logical vectors (DATAPTR, copy on modify,
// standard serialisation) and are freed by a finalizer. The number of live
// bytes is tracked so that read_h5ad() can run a GC itself when earlier
// buffers might be garbage.
#include "common.h"

#include <R_ext/Altrep.h>

#include <cstdlib>

namespace annr {

namespace {

R_altrep_class_t cls_real, cls_int, cls_lgl;
std::atomic<size_t> live_bytes{0};
bool classes_ready = false;

struct ExtBuf {
  void* data;
  R_xlen_t n;
  size_t bytes;
};

void ext_finalize(SEXP xp) {
  ExtBuf* b = static_cast<ExtBuf*>(R_ExternalPtrAddr(xp));
  if (b != nullptr) {
    live_bytes -= b->bytes;
    std::free(b->data);
    delete b;
    R_ClearExternalPtr(xp);
  }
}

inline ExtBuf* buf_of(SEXP x) { return static_cast<ExtBuf*>(R_ExternalPtrAddr(R_altrep_data1(x))); }

R_xlen_t ext_length(SEXP x) { return buf_of(x)->n; }
void* ext_dataptr(SEXP x, Rboolean) { return buf_of(x)->data; }
const void* ext_dataptr_or_null(SEXP x) { return buf_of(x)->data; }

Rboolean ext_inspect(SEXP x, int, int, int, void (*)(SEXP, int, int, int)) {
  Rprintf(" AnnR buffer (%.1f MB)\n", static_cast<double>(buf_of(x)->bytes) / 1048576.0);
  return TRUE;
}

// No Serialized_state method: objects are serialised as standard vectors, so
// saved files never depend on AnnR.
void register_methods(R_altrep_class_t cls) {
  R_set_altrep_Length_method(cls, ext_length);
  R_set_altrep_Inspect_method(cls, ext_inspect);
  R_set_altvec_Dataptr_method(cls, ext_dataptr);
  R_set_altvec_Dataptr_or_null_method(cls, ext_dataptr_or_null);
}

}  // namespace

void register_altrep_classes(DllInfo* dll) {
  cls_real = R_make_altreal_class("annr_real", "AnnR", dll);
  cls_int = R_make_altinteger_class("annr_int", "AnnR", dll);
  cls_lgl = R_make_altlogical_class("annr_lgl", "AnnR", dll);
  register_methods(cls_real);
  register_methods(cls_int);
  register_methods(cls_lgl);
  classes_ready = true;
}

double external_live_bytes() { return static_cast<double>(live_bytes.load()); }

// Threshold above which vectors are allocated outside R's heap.
constexpr R_xlen_t EXT_MIN_LENGTH = R_xlen_t(1) << 16;

SEXP alloc_vector(SEXPTYPE type, R_xlen_t n) {
  if (!classes_ready || n < EXT_MIN_LENGTH || (type != REALSXP && type != INTSXP && type != LGLSXP)) {
    return cpp11::safe[Rf_allocVector](type, n);
  }
  const size_t es = type == REALSXP ? sizeof(double) : sizeof(int);
  const size_t bytes = static_cast<size_t>(n) * es;
  void* mem = std::malloc(bytes);
  if (mem == nullptr) throw std::runtime_error("cannot allocate " + std::to_string(bytes >> 20) + " MB");
  ExtBuf* b = new ExtBuf{mem, n, bytes};
  live_bytes += bytes;
  bool owned = false;  // set once the finalizer owns the buffer
  try {
    return cpp11::unwind_protect([&] {
      SEXP xp = PROTECT(R_MakeExternalPtr(b, R_NilValue, R_NilValue));
      R_RegisterCFinalizerEx(xp, ext_finalize, TRUE);
      owned = true;
      R_altrep_class_t cls = type == REALSXP ? cls_real : (type == INTSXP ? cls_int : cls_lgl);
      SEXP res = R_new_altrep(cls, xp, R_NilValue);
      UNPROTECT(1);
      return res;
    });
  } catch (...) {
    if (!owned) {
      live_bytes -= bytes;
      std::free(mem);
      delete b;
    }
    throw;
  }
}

}  // namespace annr
