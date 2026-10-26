// Shared declarations for the AnnR engine.
//
// Threading rule: worker threads never touch the R API or HDF5. They only
// read the file through pread(), decode, and convert into memory that the
// main thread allocated beforehand.
#pragma once

#include <cpp11.hpp>
#include <hdf5.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if !H5_VERSION_GE(1, 10, 5)
#error "AnnR requires HDF5 >= 1.10.5"
#endif

namespace annr {

// ---------------------------------------------------------------------------
// Optional timing output (ANNR_PROFILE=1), printed from the main thread only.

inline bool profile_enabled() {
  static const bool on = [] {
    const char* v = std::getenv("ANNR_PROFILE");
    return v != nullptr && v[0] == '1';
  }();
  return on;
}

struct ProfileScope {
  const char* what;
  size_t elems;
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  ~ProfileScope() {
    if (!profile_enabled()) return;
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    REprintf("[AnnR] %-10s %12zu elems %9.2f ms\n", what, elems, ms);
  }
};

// ---------------------------------------------------------------------------
// Errors

std::string h5_error_stack();
[[noreturn]] void h5_fail(const std::string& what);

// Raised when the direct (parallel) path cannot handle a dataset; the caller
// then falls back to H5Dread, which is always correct.
struct DirectFallback : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// ---------------------------------------------------------------------------
// RAII owner of an HDF5 identifier.

class H5Handle {
 public:
  using closer_t = herr_t (*)(hid_t);
  H5Handle() = default;
  H5Handle(hid_t id, closer_t closer, const char* what) : id_(id), close_(closer) {
    if (id_ < 0) h5_fail(what);
  }
  ~H5Handle() { reset(); }
  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;
  H5Handle(H5Handle&& o) noexcept : id_(o.id_), close_(o.close_) { o.id_ = H5I_INVALID_HID; }
  H5Handle& operator=(H5Handle&& o) noexcept {
    if (this != &o) {
      reset();
      id_ = o.id_;
      close_ = o.close_;
      o.id_ = H5I_INVALID_HID;
    }
    return *this;
  }
  void reset() {
    if (id_ >= 0 && close_ != nullptr) close_(id_);
    id_ = H5I_INVALID_HID;
  }
  hid_t get() const { return id_; }
  operator hid_t() const { return id_; }

 private:
  hid_t id_ = H5I_INVALID_HID;
  closer_t close_ = nullptr;
};

H5Handle open_object(hid_t loc, const std::string& name);
H5Handle open_dataset(hid_t loc, const std::string& name);
H5Handle open_group(hid_t loc, const std::string& name);
H5Handle dataset_type(hid_t dset);
H5Handle dataset_space(hid_t dset);
bool link_exists(hid_t loc, const std::string& name);
H5I_type_t object_kind(hid_t obj);
std::vector<std::string> group_children(hid_t group);
std::vector<hsize_t> space_dims(hid_t space);
bool space_is_null(hid_t space);

// ---------------------------------------------------------------------------
// Type model

// Primitive element types the direct path can decode without HDF5.
enum class Prim : int { I8, U8, I16, U16, I32, U32, I64, U64, F16, F32, F64, BOOL, NONE };

struct TypeDesc {
  H5T_class_t cls = H5T_NO_CLASS;
  size_t size = 0;
  Prim prim = Prim::NONE;   // set for little-endian standard numeric types
  bool is_numeric = false;  // integer, float or h5py bool
  bool is_bool = false;     // h5py-style enum {FALSE=0, TRUE=1}
  bool is_string = false;
  bool is_vlen_string = false;
  bool is_complex = false;  // h5py complex compound {r, i}
  bool is_compound = false;
  bool is_signed = true;
  H5T_cset_t cset = H5T_CSET_ASCII;
  H5T_str_t strpad = H5T_STR_NULLTERM;
  std::string dtype;        // numpy-style name, e.g. "float32", "str", "S7"
};

TypeDesc describe_type(hid_t type);
size_t prim_size(Prim p);
bool prim_is_float(Prim p);

// ---------------------------------------------------------------------------
// Read context

struct TypeRecord {
  std::string path, encoding, version, dtype, extra;
};

struct Ctx {
  std::string filename;
  int fd = -1;          // descriptor for direct reads, -1 disables the direct path
  int nthreads = 1;
  bool strict = false;  // turn direct-path fallbacks into errors (testing)
  size_t direct_reads = 0;
  size_t fallback_reads = 0;
  std::vector<TypeRecord> types;
  std::vector<std::string> warnings;
  std::vector<std::string> fallback_notes;

  void record(const std::string& path, const std::string& enc, const std::string& ver,
              const std::string& dtype, const std::string& extra = "") {
    types.push_back({path, enc, ver, dtype, extra});
  }
  void warn(const std::string& msg) { warnings.push_back(msg); }
};

// ---------------------------------------------------------------------------
// Attributes

bool has_attr(hid_t obj, const char* name);
std::string attr_string(hid_t obj, const char* name, const std::string& dflt = "");
std::vector<std::string> attr_strings(hid_t obj, const char* name);
std::vector<int64_t> attr_int64s(hid_t obj, const char* name);
bool attr_bool(hid_t obj, const char* name, bool dflt);

// ---------------------------------------------------------------------------
// Parallel loop. body(i, thread_id) must not call R or HDF5.

template <typename Body>
void parallel_for(size_t n, int nthreads, Body&& body) {
  if (n == 0) return;
  if (nthreads <= 1 || n == 1) {
    for (size_t i = 0; i < n; ++i) body(i, 0);
    return;
  }
  std::atomic<size_t> next{0};
  std::atomic<bool> failed{false};
  bool is_fallback = false;
  std::string err;
  std::mutex mtx;
  auto worker = [&](int tid) {
    try {
      for (;;) {
        if (failed.load(std::memory_order_relaxed)) break;
        size_t i = next.fetch_add(1, std::memory_order_relaxed);
        if (i >= n) break;
        body(i, tid);
      }
    } catch (const DirectFallback& e) {
      std::lock_guard<std::mutex> lock(mtx);
      if (!failed.exchange(true)) {
        err = e.what();
        is_fallback = true;
      }
    } catch (const std::exception& e) {
      std::lock_guard<std::mutex> lock(mtx);
      if (!failed.exchange(true)) err = e.what();
    } catch (...) {
      std::lock_guard<std::mutex> lock(mtx);
      if (!failed.exchange(true)) err = "unknown error in worker thread";
    }
  };
  size_t nt = std::min<size_t>(static_cast<size_t>(nthreads), n);
  std::vector<std::thread> pool;
  pool.reserve(nt - 1);
  for (size_t t = 1; t < nt; ++t) {
    try {
      pool.emplace_back(worker, static_cast<int>(t));
    } catch (...) {
      break;  // could not start more threads; continue with what we have
    }
  }
  worker(0);
  for (auto& th : pool) th.join();
  if (failed) {
    if (is_fallback) throw DirectFallback(err);
    throw std::runtime_error(err);
  }
}

// Splits [0, n) into contiguous blocks so that each thread gets a share.
inline size_t block_count(size_t n, size_t min_block, int nthreads) {
  if (n == 0) return 0;
  size_t by_size = (n + min_block - 1) / min_block;
  size_t want = static_cast<size_t>(nthreads) * 4;
  return std::max<size_t>(1, std::min(by_size, want));
}

// ---------------------------------------------------------------------------
// Numeric reading engine (chunkio.cpp)

// Destination element types in R memory.
enum class Dst { F64, I32, LGL, I64 };

// Reads the complete dataset into `out` (n elements of `dst`).
// colmajor=true lays out rank>=2 data column-major with R dims equal to the
// HDF5 dims (an axis reversal); otherwise the C-order element sequence is
// kept, which equals an R array with reversed dims.
void read_numeric(hid_t dset, const TypeDesc& t, Dst dst, void* out, bool colmajor, Ctx& ctx);

// Reads a 1-D hyperslab [start, start+count) of a numeric dataset.
void read_numeric_range(hid_t dset, const TypeDesc& t, Dst dst, void* out, hsize_t start,
                        hsize_t count, Ctx& ctx);

// ---------------------------------------------------------------------------
// Codecs (codecs.cpp)

size_t lzf_decompress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len);
size_t lzf_compress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len);
void unshuffle(const uint8_t* in, uint8_t* out, size_t nbytes, size_t elem);
void shuffle(const uint8_t* in, uint8_t* out, size_t nbytes, size_t elem);
bool inflate_zlib(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len, void** state);
void free_inflate_state(void* state);
size_t deflate_zlib(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_cap, int level,
                    void** state);
size_t deflate_bound(size_t in_len);
void free_deflate_state(void* state);
bool zstd_available();
bool zstd_decompress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len);
size_t zstd_compress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_cap, int level);
size_t zstd_bound(size_t in_len);
void register_filters();

constexpr H5Z_filter_t FILTER_LZF = 32000;
constexpr H5Z_filter_t FILTER_ZSTD = 32015;

// ---------------------------------------------------------------------------
// R helpers

// Allocates an R vector with longjmp protection (throws on failure). Large
// numeric vectors live outside R's heap (altrep.cpp).
SEXP alloc_vector(SEXPTYPE type, R_xlen_t n);
void register_altrep_classes(DllInfo* dll);
double external_live_bytes();
SEXP make_s4(const char* cls);
void set_dim(SEXP x, const std::vector<hsize_t>& dims);
std::string join_path(const std::string& parent, const std::string& name);

}  // namespace annr
