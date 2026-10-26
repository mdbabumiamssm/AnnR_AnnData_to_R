// Parallel direct reader for numeric datasets.
//
// The HDF5 library is used only for metadata (layout, filters, chunk
// addresses). Raw bytes are fetched with pread() by worker threads, the filter
// pipeline is decoded in-thread, and values are converted (and optionally
// transposed) straight into memory owned by an R vector. Anything outside the
// supported envelope falls back to H5Dread, which is always correct.
#include "common.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <type_traits>
#ifndef _WIN32
#include <unistd.h>
#endif

namespace annr {

namespace {

struct Filter {
  H5Z_filter_t id;
  unsigned flags;
  std::vector<unsigned> cd;
};

struct Layout {
  int rank = 0;
  std::vector<hsize_t> dims;
  H5D_layout_t layout = H5D_LAYOUT_ERROR;
  std::vector<hsize_t> chunk;
  std::vector<Filter> filters;
  bool external = false;
};

size_t product(const std::vector<hsize_t>& v) {
  size_t p = 1;
  for (hsize_t x : v) p *= static_cast<size_t>(x);
  return p;
}

Layout describe_layout(hid_t dset) {
  Layout L;
  H5Handle space = dataset_space(dset);
  L.dims = space_dims(space);
  L.rank = static_cast<int>(L.dims.size());
  H5Handle dcpl(H5Dget_create_plist(dset), H5Pclose, "H5Dget_create_plist");
  L.layout = H5Pget_layout(dcpl);
  if (L.layout == H5D_CHUNKED) {
    L.chunk.resize(static_cast<size_t>(L.rank));
    if (H5Pget_chunk(dcpl, L.rank, L.chunk.data()) < 0) h5_fail("H5Pget_chunk");
  }
  int nf = H5Pget_nfilters(dcpl);
  for (int i = 0; i < nf; ++i) {
    unsigned flags = 0, config = 0;
    size_t ncd = 16;
    std::vector<unsigned> cd(16);
    char name[64];
    H5Z_filter_t id = H5Pget_filter2(dcpl, static_cast<unsigned>(i), &flags, &ncd, cd.data(),
                                     sizeof(name), name, &config);
    if (id < 0) h5_fail("H5Pget_filter2");
    cd.resize(std::min<size_t>(ncd, 16));
    L.filters.push_back({id, flags, cd});
  }
  int next = H5Pget_external_count(dcpl);
  L.external = next > 0;
  H5Eclear2(H5E_DEFAULT);
  return L;
}

bool filters_supported(const std::vector<Filter>& fs) {
  for (const auto& f : fs) {
    switch (f.id) {
      case H5Z_FILTER_DEFLATE:
      case H5Z_FILTER_SHUFFLE:
      case FILTER_LZF:
        break;
      case FILTER_ZSTD:
        if (!zstd_available()) return false;
        break;
      default:
        return false;
    }
  }
  return true;
}

// ---------------------------------------------------------------------------
// Element conversion kernels

template <typename S>
inline S ld(const uint8_t* p) {
  S v;
  std::memcpy(&v, p, sizeof(S));
  return v;
}

inline double half_to_double(uint16_t h) {
  uint32_t sign = (h >> 15) & 1u, exp = (h >> 10) & 0x1fu, mant = h & 0x3ffu;
  double v;
  if (exp == 0) {
    v = std::ldexp(static_cast<double>(mant), -24);
  } else if (exp == 31) {
    v = mant == 0 ? INFINITY : NAN;
  } else {
    v = std::ldexp(static_cast<double>(mant | 0x400u), static_cast<int>(exp) - 25);
  }
  return sign ? -v : v;
}

struct F16 {
  uint16_t bits;
};

template <typename S>
inline double to_f64(S v) {
  return static_cast<double>(v);
}
template <>
inline double to_f64<F16>(F16 v) {
  return half_to_double(v.bits);
}

// Out-of-range values map to INT_MIN, which every consumer of I32 rejects or
// treats as missing (sparse indices, categorical codes).
template <typename S>
inline int32_t to_i32(S v) {
  if constexpr (std::is_floating_point_v<S> || std::is_same_v<S, F16>) {
    return INT32_MIN;
  } else if constexpr (std::is_signed_v<S>) {
    if constexpr (sizeof(S) <= 4) {
      return static_cast<int32_t>(v);
    } else {
      return (v > INT32_MAX || v < -INT32_MAX) ? INT32_MIN : static_cast<int32_t>(v);
    }
  } else {
    if constexpr (sizeof(S) < 4) {
      return static_cast<int32_t>(v);
    } else {
      return v > static_cast<S>(INT32_MAX) ? INT32_MIN : static_cast<int32_t>(v);
    }
  }
}

template <typename S>
inline int32_t to_lgl(S v) {
  if constexpr (std::is_same_v<S, F16>) {
    return v.bits != 0;
  } else {
    return v != 0 ? 1 : 0;
  }
}

// I64 keeps the raw 64-bit pattern for unsigned 64-bit input; the caller
// interprets it according to the source type.
template <typename S>
inline int64_t to_i64(S v) {
  if constexpr (std::is_same_v<S, F16> || std::is_floating_point_v<S>) {
    return 0;
  } else if constexpr (std::is_same_v<S, uint64_t>) {
    int64_t out;
    std::memcpy(&out, &v, 8);
    return out;
  } else {
    return static_cast<int64_t>(v);
  }
}

using ConvFn = void (*)(const uint8_t* src, size_t sstride, void* dst, size_t dstride, size_t n);

template <typename D, typename S, D (*F)(S)>
void conv_kernel(const uint8_t* src, size_t sstride, void* dstv, size_t dstride, size_t n) {
  D* dst = static_cast<D*>(dstv);
  if (sstride == sizeof(S) && dstride == 1) {
    for (size_t i = 0; i < n; ++i) dst[i] = F(ld<S>(src + i * sizeof(S)));
  } else {
    for (size_t i = 0; i < n; ++i) dst[i * dstride] = F(ld<S>(src + i * sstride));
  }
}

template <typename S>
ConvFn conv_for_dst(Dst d) {
  switch (d) {
    case Dst::F64: return conv_kernel<double, S, to_f64<S>>;
    case Dst::I32: return conv_kernel<int32_t, S, to_i32<S>>;
    case Dst::LGL: return conv_kernel<int32_t, S, to_lgl<S>>;
    case Dst::I64: return conv_kernel<int64_t, S, to_i64<S>>;
  }
  return nullptr;
}

ConvFn get_conv(Prim sp, Dst d) {
  switch (sp) {
    case Prim::I8: return conv_for_dst<int8_t>(d);
    case Prim::U8: return conv_for_dst<uint8_t>(d);
    case Prim::BOOL: return conv_for_dst<uint8_t>(d);
    case Prim::I16: return conv_for_dst<int16_t>(d);
    case Prim::U16: return conv_for_dst<uint16_t>(d);
    case Prim::I32: return conv_for_dst<int32_t>(d);
    case Prim::U32: return conv_for_dst<uint32_t>(d);
    case Prim::I64: return conv_for_dst<int64_t>(d);
    case Prim::U64: return conv_for_dst<uint64_t>(d);
    case Prim::F16: return d == Dst::F64 || d == Dst::LGL ? conv_for_dst<F16>(d) : nullptr;
    case Prim::F32: return d == Dst::F64 || d == Dst::LGL ? conv_for_dst<float>(d) : nullptr;
    case Prim::F64: return d == Dst::F64 || d == Dst::LGL ? conv_for_dst<double>(d) : nullptr;
    default: return nullptr;
  }
}

size_t dst_size(Dst d) {
  switch (d) {
    case Dst::F64: return 8;
    case Dst::I32: return 4;
    case Dst::LGL: return 4;
    case Dst::I64: return 8;
  }
  return 0;
}

bool is_identity(Prim sp, Dst d) {
  return (d == Dst::F64 && sp == Prim::F64) || (d == Dst::I32 && sp == Prim::I32) ||
         (d == Dst::I64 && (sp == Prim::I64 || sp == Prim::U64));
}

// ---------------------------------------------------------------------------
// Scatter one block (chunk) into the output array.
//
// src holds a C-order block with dims `cd` whose first element sits at
// dataset coordinate `org`; `ext` is the valid extent (clipped at edges).

struct Scatter {
  ConvFn conv;
  size_t es;     // source element size
  size_t ds;     // destination element size
  int rank;
  const hsize_t* D;
  uint8_t* out;
  bool colmajor;

  void operator()(const uint8_t* src, const hsize_t* cd, const hsize_t* org,
                  const hsize_t* ext) const {
    if (rank == 1) {
      conv(src, es, out + org[0] * ds, 1, ext[0]);
      return;
    }
    if (!colmajor) {
      // C-order output: rows of the innermost dimension are contiguous.
      std::vector<size_t> dstride(rank), cstride(rank);
      dstride[rank - 1] = 1;
      cstride[rank - 1] = 1;
      for (int k = rank - 2; k >= 0; --k) {
        dstride[k] = dstride[k + 1] * D[k + 1];
        cstride[k] = cstride[k + 1] * cd[k + 1];
      }
      std::vector<hsize_t> idx(rank, 0);
      size_t run = ext[rank - 1];
      for (;;) {
        size_t so = 0, doff = 0;
        for (int k = 0; k < rank - 1; ++k) {
          so += idx[k] * cstride[k];
          doff += (org[k] + idx[k]) * dstride[k];
        }
        doff += org[rank - 1];
        conv(src + so * es, es, out + doff * ds, 1, run);
        int k = rank - 2;
        while (k >= 0) {
          if (++idx[k] < ext[k]) break;
          idx[k] = 0;
          --k;
        }
        if (k < 0) break;
      }
      return;
    }
    if (rank == 2) {
      // Column-major output (axis reversal) with cache tiling.
      constexpr size_t T = 64;
      const size_t nrow = D[0];
      const size_t sstride = cd[1] * es;
      for (size_t r0 = 0; r0 < ext[0]; r0 += T) {
        size_t rn = std::min(T, static_cast<size_t>(ext[0]) - r0);
        for (size_t c = 0; c < ext[1]; ++c) {
          const uint8_t* s = src + (r0 * cd[1] + c) * es;
          size_t doff = (org[1] + c) * nrow + org[0] + r0;
          conv(s, sstride, out + doff * ds, 1, rn);
        }
      }
      return;
    }
    // Generic N-D axis reversal, one element at a time.
    std::vector<size_t> fstride(rank);
    fstride[0] = 1;
    for (int k = 1; k < rank; ++k) fstride[k] = fstride[k - 1] * D[k - 1];
    std::vector<size_t> cstride(rank);
    cstride[rank - 1] = 1;
    for (int k = rank - 2; k >= 0; --k) cstride[k] = cstride[k + 1] * cd[k + 1];
    std::vector<hsize_t> idx(rank, 0);
    for (;;) {
      size_t so = 0, doff = 0;
      for (int k = 0; k < rank; ++k) {
        so += idx[k] * cstride[k];
        doff += (org[k] + idx[k]) * fstride[k];
      }
      conv(src + so * es, es, out + doff * ds, 1, 1);
      int k = rank - 1;
      while (k >= 0) {
        if (++idx[k] < ext[k]) break;
        idx[k] = 0;
        --k;
      }
      if (k < 0) break;
    }
  }
};

void pread_full(int fd, uint8_t* buf, size_t n, uint64_t off) {
#ifdef _WIN32
  (void)fd; (void)buf; (void)n; (void)off;
  throw DirectFallback("direct reads are not available on Windows");
#else
  while (n > 0) {
    ssize_t r = ::pread(fd, buf, n, static_cast<off_t>(off));
    if (r < 0) {
      if (errno == EINTR) continue;
      throw DirectFallback("pread failed");
    }
    if (r == 0) throw DirectFallback("unexpected end of file");
    buf += r;
    n -= static_cast<size_t>(r);
    off += static_cast<uint64_t>(r);
  }
#endif
}

// ---------------------------------------------------------------------------
// Chunk enumeration

struct ChunkSet {
  std::vector<hsize_t> offsets;  // nchunks * rank, element coordinates
  std::vector<haddr_t> addr;
  std::vector<hsize_t> size;
  std::vector<unsigned> mask;
  size_t count() const { return addr.size(); }
};

#if H5_VERSION_GE(1, 14, 0)
struct IterState {
  ChunkSet* cs;
  int rank;
};

int chunk_iter_cb(const hsize_t* offset, unsigned filter_mask, haddr_t addr, hsize_t size,
                  void* op_data) {
  auto* st = static_cast<IterState*>(op_data);
  for (int k = 0; k < st->rank; ++k) st->cs->offsets.push_back(offset[k]);
  st->cs->addr.push_back(addr);
  st->cs->size.push_back(size);
  st->cs->mask.push_back(filter_mask);
  return H5_ITER_CONT;
}
#endif

ChunkSet enumerate_chunks(hid_t dset, const Layout& L) {
  ChunkSet cs;
  const int rank = L.rank;
  std::vector<size_t> grid(rank);
  size_t expected = 1;
  for (int k = 0; k < rank; ++k) {
    grid[k] = (L.dims[k] + L.chunk[k] - 1) / L.chunk[k];
    expected *= grid[k];
  }
#if H5_VERSION_GE(1, 14, 0)
  IterState st{&cs, rank};
  if (H5Dchunk_iter(dset, H5P_DEFAULT, chunk_iter_cb, &st) < 0) {
    H5Eclear2(H5E_DEFAULT);
    throw DirectFallback("H5Dchunk_iter failed");
  }
  // Some 1.14.0 builds report offsets in chunk units; detect and rescale.
  bool scaled = false;
  for (size_t i = 0; i < cs.count() && !scaled; ++i) {
    for (int k = 0; k < rank; ++k) {
      if (cs.offsets[i * rank + k] % L.chunk[k] != 0) {
        scaled = true;
        break;
      }
    }
  }
  if (scaled) {
    for (size_t i = 0; i < cs.count(); ++i) {
      for (int k = 0; k < rank; ++k) cs.offsets[i * rank + k] *= L.chunk[k];
    }
  }
#else
  std::vector<hsize_t> coord(rank, 0);
  for (size_t c = 0; c < expected; ++c) {
    size_t rem = c;
    for (int k = rank - 1; k >= 0; --k) {
      coord[k] = (rem % grid[k]) * L.chunk[k];
      rem /= grid[k];
    }
    unsigned fmask = 0;
    haddr_t addr = HADDR_UNDEF;
    hsize_t sz = 0;
    if (H5Dget_chunk_info_by_coord(dset, coord.data(), &fmask, &addr, &sz) < 0) {
      H5Eclear2(H5E_DEFAULT);
      throw DirectFallback("H5Dget_chunk_info_by_coord failed");
    }
    if (addr == HADDR_UNDEF) throw DirectFallback("unallocated chunk");
    for (int k = 0; k < rank; ++k) cs.offsets.push_back(coord[k]);
    cs.addr.push_back(addr);
    cs.size.push_back(sz);
    cs.mask.push_back(fmask);
  }
#endif
  if (cs.count() != expected) throw DirectFallback("dataset has unallocated chunks");
  // Every grid cell must be covered exactly once.
  std::vector<uint8_t> seen(expected, 0);
  for (size_t i = 0; i < cs.count(); ++i) {
    size_t lin = 0;
    for (int k = 0; k < rank; ++k) {
      hsize_t o = cs.offsets[i * rank + k];
      if (o % L.chunk[k] != 0 || o >= L.dims[k]) throw DirectFallback("unexpected chunk offset");
      lin = lin * grid[k] + o / L.chunk[k];
    }
    if (seen[lin]) throw DirectFallback("duplicate chunk");
    seen[lin] = 1;
    if (cs.addr[i] == HADDR_UNDEF || cs.size[i] == 0) throw DirectFallback("invalid chunk record");
  }
  return cs;
}

struct ThreadBuf {
  std::vector<uint8_t> raw, a, b;
  void* inflate_state = nullptr;
  ~ThreadBuf() { free_inflate_state(inflate_state); }
};

// ---------------------------------------------------------------------------
// Direct read

void direct_read(hid_t dset, const Layout& L, const TypeDesc& t, Dst dst, void* out,
                 bool colmajor, Ctx& ctx) {
  const size_t es = prim_size(t.prim);
  if (es == 0 || es != t.size) throw DirectFallback("unsupported element size");
  ConvFn conv = get_conv(t.prim, dst);
  if (conv == nullptr) throw DirectFallback("unsupported conversion");
  const size_t ds = dst_size(dst);
  const int rank = L.rank;
  const bool layout_identity = is_identity(t.prim, dst) && (rank == 1 || !colmajor);
  uint8_t* outb = static_cast<uint8_t*>(out);
  Scatter scatter{conv, es, ds, rank, L.dims.data(), outb, colmajor};
  const int nthreads = std::max(1, ctx.nthreads);

  if (L.layout == H5D_CONTIGUOUS) {
    haddr_t base = H5Dget_offset(dset);
    if (base == HADDR_UNDEF) {
      H5Eclear2(H5E_DEFAULT);
      throw DirectFallback("contiguous storage not allocated");
    }
    const size_t n = product(L.dims);
    if (H5Dget_storage_size(dset) != static_cast<hsize_t>(n * es)) {
      throw DirectFallback("unexpected contiguous storage size");
    }
    if (layout_identity) {
      // Bytes on disk are exactly the output bytes.
      const size_t total = n * es;
      const size_t blk = size_t(8) << 20;
      const size_t nb = (total + blk - 1) / blk;
      parallel_for(nb, nthreads, [&](size_t b, int) {
        size_t off = b * blk;
        size_t len = std::min(blk, total - off);
        pread_full(ctx.fd, outb + off, len, base + off);
      });
      return;
    }
    // Virtual chunks of whole rows along the first dimension.
    const size_t row_elems = rank > 1 ? n / L.dims[0] : 1;
    const size_t row_bytes = row_elems * es;
    const size_t target = size_t(4) << 20;
    size_t rows_per = std::max<size_t>(1, target / std::max<size_t>(1, row_bytes));
    if (rank == 1) rows_per = std::max<size_t>(rows_per, 1);
    const size_t nrows = L.dims[0];
    const size_t nb = (nrows + rows_per - 1) / rows_per;
    std::vector<ThreadBuf> tb(static_cast<size_t>(nthreads));
    parallel_for(nb, nthreads, [&](size_t b, int tid) {
      size_t r0 = b * rows_per;
      size_t rn = std::min(rows_per, nrows - r0);
      auto& buf = tb[static_cast<size_t>(tid)].raw;
      buf.resize(rn * row_bytes);
      pread_full(ctx.fd, buf.data(), rn * row_bytes, base + r0 * row_bytes);
      std::vector<hsize_t> cd(L.dims), org(rank, 0), ext(L.dims);
      cd[0] = rows_per;
      org[0] = r0;
      ext[0] = rn;
      scatter(buf.data(), cd.data(), org.data(), ext.data());
    });
    return;
  }

  if (L.layout != H5D_CHUNKED) throw DirectFallback("unsupported layout");
  ChunkSet cs = enumerate_chunks(dset, L);
  const size_t chunk_elems = product(L.chunk);
  const size_t chunk_bytes = chunk_elems * es;
  hsize_t max_stored = 0;
  for (hsize_t s : cs.size) max_stored = std::max(max_stored, s);

  // Visit chunks in file order so the kernel's read-ahead helps.
  std::vector<size_t> order(cs.count());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return cs.addr[a] < cs.addr[b]; });

  const int nf = static_cast<int>(L.filters.size());
  std::vector<ThreadBuf> tb(static_cast<size_t>(nthreads));
  parallel_for(order.size(), nthreads, [&](size_t oi, int tid) {
    const size_t i = order[oi];
    ThreadBuf& T = tb[static_cast<size_t>(tid)];
    const hsize_t* org = &cs.offsets[i * rank];
    std::vector<hsize_t> ext(rank);
    for (int k = 0; k < rank; ++k) ext[k] = std::min<hsize_t>(L.chunk[k], L.dims[k] - org[k]);
    bool any_filter = false;
    for (int f = 0; f < nf; ++f) {
      if (!(cs.mask[i] & (1u << f))) any_filter = true;
    }
    if (!any_filter && cs.size[i] != chunk_bytes) throw DirectFallback("unfiltered chunk size mismatch");
    if (!any_filter && rank == 1 && layout_identity) {
      pread_full(ctx.fd, outb + org[0] * ds, ext[0] * es, cs.addr[i]);
      return;
    }
    T.raw.resize(static_cast<size_t>(max_stored));
    pread_full(ctx.fd, T.raw.data(), static_cast<size_t>(cs.size[i]), cs.addr[i]);
    const uint8_t* cur = T.raw.data();
    size_t cur_len = static_cast<size_t>(cs.size[i]);
    if (any_filter) {
      T.a.resize(chunk_bytes);
      T.b.resize(chunk_bytes);
      for (int f = nf - 1; f >= 0; --f) {
        if (cs.mask[i] & (1u << f)) continue;
        const Filter& F = L.filters[static_cast<size_t>(f)];
        uint8_t* dstbuf = (cur == T.a.data()) ? T.b.data() : T.a.data();
        switch (F.id) {
          case H5Z_FILTER_DEFLATE:
            if (!inflate_zlib(cur, cur_len, dstbuf, chunk_bytes, &T.inflate_state)) {
              throw DirectFallback("inflate failed");
            }
            cur_len = chunk_bytes;
            break;
          case H5Z_FILTER_SHUFFLE: {
            size_t fes = F.cd.empty() ? es : F.cd[0];
            if (cur_len != chunk_bytes) throw DirectFallback("shuffle size mismatch");
            unshuffle(cur, dstbuf, cur_len, fes);
            break;
          }
          case FILTER_LZF: {
            size_t got = lzf_decompress(cur, cur_len, dstbuf, chunk_bytes);
            if (got != chunk_bytes) throw DirectFallback("lzf decode failed");
            cur_len = got;
            break;
          }
          case FILTER_ZSTD:
            if (!zstd_decompress(cur, cur_len, dstbuf, chunk_bytes)) {
              throw DirectFallback("zstd decode failed");
            }
            cur_len = chunk_bytes;
            break;
          default:
            throw DirectFallback("unsupported filter");
        }
        cur = dstbuf;
      }
    }
    if (cur_len != chunk_bytes) throw DirectFallback("decoded chunk size mismatch");
    scatter(cur, L.chunk.data(), org, ext.data());
  });
}

// ---------------------------------------------------------------------------
// HDF5 fallback

struct MemType {
  H5Handle type;
  Prim prim;
};

MemType mem_type_for(hid_t ftype, const TypeDesc& t, Dst dst) {
  switch (dst) {
    case Dst::F64:
      return {H5Handle(H5Tcopy(H5T_NATIVE_DOUBLE), H5Tclose, "H5Tcopy"), Prim::F64};
    case Dst::I32:
      return {H5Handle(H5Tcopy(H5T_NATIVE_INT32), H5Tclose, "H5Tcopy"), Prim::I32};
    case Dst::LGL:
      if (t.is_bool) {
        return {H5Handle(H5Tget_native_type(ftype, H5T_DIR_ASCEND), H5Tclose, "native enum"),
                Prim::U8};
      }
      return {H5Handle(H5Tcopy(H5T_NATIVE_INT32), H5Tclose, "H5Tcopy"), Prim::I32};
    case Dst::I64:
      if (!t.is_signed && t.size == 8) {
        return {H5Handle(H5Tcopy(H5T_NATIVE_UINT64), H5Tclose, "H5Tcopy"), Prim::U64};
      }
      return {H5Handle(H5Tcopy(H5T_NATIVE_INT64), H5Tclose, "H5Tcopy"), Prim::I64};
  }
  throw std::runtime_error("bad destination type");
}

void hdf5_read(hid_t dset, const Layout& L, const TypeDesc& t, Dst dst, void* out,
               bool colmajor) {
  H5Handle ftype = dataset_type(dset);
  MemType mt = mem_type_for(ftype, t, dst);
  const size_t n = product(L.dims);
  const bool direct_ok = (dst != Dst::LGL || !t.is_bool) && (L.rank < 2 || !colmajor);
  if (direct_ok) {
    if (H5Dread(dset, mt.type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out) < 0) h5_fail("H5Dread failed");
    return;
  }
  const size_t msz = H5Tget_size(mt.type);
  std::vector<uint8_t> tmp(std::max<size_t>(1, n * msz));
  if (H5Dread(dset, mt.type, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp.data()) < 0) h5_fail("H5Dread failed");
  ConvFn conv = get_conv(mt.prim, dst);
  if (conv == nullptr) throw std::runtime_error("unsupported conversion in fallback read");
  std::vector<hsize_t> org(L.rank, 0);
  Scatter sc{conv, msz, dst_size(dst), L.rank, L.dims.data(), static_cast<uint8_t*>(out), colmajor};
  if (L.rank == 0) {
    conv(tmp.data(), msz, out, 1, 1);
  } else {
    sc(tmp.data(), L.dims.data(), org.data(), L.dims.data());
  }
}

}  // namespace

void read_numeric(hid_t dset, const TypeDesc& t, Dst dst, void* out, bool colmajor, Ctx& ctx) {
  Layout L = describe_layout(dset);
  ProfileScope prof{"numeric", L.rank == 0 ? 1 : product(L.dims)};
  const size_t n = L.rank == 0 ? 1 : product(L.dims);
  if (n == 0) return;
  const bool can_direct = ctx.fd >= 0 && t.prim != Prim::NONE && !L.external && L.rank >= 1 &&
                          (L.layout == H5D_CONTIGUOUS ||
                           (L.layout == H5D_CHUNKED && filters_supported(L.filters)));
  if (can_direct) {
    try {
      direct_read(dset, L, t, dst, out, colmajor, ctx);
      ctx.direct_reads++;
      return;
    } catch (const DirectFallback& e) {
      if (ctx.strict) throw std::runtime_error(std::string("direct read failed: ") + e.what());
      ctx.fallback_notes.emplace_back(e.what());
    }
  }
  hdf5_read(dset, L, t, dst, out, colmajor);
  ctx.fallback_reads++;
}

void read_numeric_range(hid_t dset, const TypeDesc& t, Dst dst, void* out, hsize_t start,
                        hsize_t count, Ctx&) {
  if (count == 0) return;
  H5Handle ftype = dataset_type(dset);
  MemType mt = mem_type_for(ftype, t, dst);
  H5Handle fspace = dataset_space(dset);
  if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0) {
    h5_fail("H5Sselect_hyperslab");
  }
  H5Handle mspace(H5Screate_simple(1, &count, nullptr), H5Sclose, "H5Screate_simple");
  if (dst == Dst::LGL && t.is_bool) {
    std::vector<uint8_t> tmp(static_cast<size_t>(count));
    if (H5Dread(dset, mt.type, mspace, fspace, H5P_DEFAULT, tmp.data()) < 0) h5_fail("H5Dread");
    int32_t* o = static_cast<int32_t*>(out);
    for (size_t i = 0; i < tmp.size(); ++i) o[i] = tmp[i] != 0;
    return;
  }
  if (H5Dread(dset, mt.type, mspace, fspace, H5P_DEFAULT, out) < 0) h5_fail("H5Dread");
}

}  // namespace annr
