// Chunk codecs: LZF, HDF5 byte-shuffle, deflate (libdeflate or zlib), zstd,
// plus registration of LZF / zstd as HDF5 filters so that H5Dread and
// H5Dcreate work with them in this process.
#include "common.h"

#include <zlib.h>
#ifdef ANNR_HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif
#ifdef ANNR_HAVE_ZSTD
#include <zstd.h>
#endif

namespace annr {

// ---------------------------------------------------------------------------
// LZF (format of Marc Lehmann's liblzf, as used by h5py's filter 32000)

size_t lzf_decompress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len) {
  const uint8_t* ip = in;
  const uint8_t* const in_end = in + in_len;
  uint8_t* op = out;
  uint8_t* const out_end = out + out_len;
  while (ip < in_end) {
    unsigned ctrl = *ip++;
    if (ctrl < 32) {  // literal run of ctrl + 1 bytes
      size_t len = static_cast<size_t>(ctrl) + 1;
      if (static_cast<size_t>(out_end - op) < len || static_cast<size_t>(in_end - ip) < len) return 0;
      std::memcpy(op, ip, len);
      op += len;
      ip += len;
    } else {  // back reference
      size_t len = ctrl >> 5;
      if (len == 7) {
        if (ip >= in_end) return 0;
        len += *ip++;
      }
      if (ip >= in_end) return 0;
      size_t back = (static_cast<size_t>(ctrl & 0x1f) << 8) + *ip++ + 1;
      len += 2;
      if (back > static_cast<size_t>(op - out) || static_cast<size_t>(out_end - op) < len) return 0;
      const uint8_t* ref = op - back;
      if (back >= len) {
        std::memcpy(op, ref, len);
      } else {
        for (size_t k = 0; k < len; ++k) op[k] = ref[k];
      }
      op += len;
    }
  }
  return static_cast<size_t>(op - out);
}

size_t lzf_compress(const uint8_t* in, size_t n, uint8_t* out, size_t cap) {
  constexpr unsigned HLOG = 14;
  constexpr size_t HSIZE = size_t(1) << HLOG;
  constexpr size_t MAX_LIT = 32, MAX_OFF = 8192, MAX_REF = 264;
  if (n == 0 || cap < 2) return 0;
  thread_local std::vector<int64_t> htab;
  htab.assign(HSIZE, -1);
  auto hash = [&](size_t p) -> size_t {
    uint32_t v = (uint32_t(in[p]) << 16) | (uint32_t(in[p + 1]) << 8) | uint32_t(in[p + 2]);
    return ((v * 2654435761u) >> (32 - HLOG)) & (HSIZE - 1);
  };
  size_t ip = 0, op = 0, lit = 0;
  size_t lit_ctrl = op++;
  while (ip < n) {
    if (ip + 2 < n) {
      size_t h = hash(ip);
      int64_t ref = htab[h];
      htab[h] = static_cast<int64_t>(ip);
      if (ref >= 0) {
        size_t r = static_cast<size_t>(ref);
        size_t off = ip - r - 1;
        if (off < MAX_OFF && in[r] == in[ip] && in[r + 1] == in[ip + 1] && in[r + 2] == in[ip + 2]) {
          size_t maxlen = std::min(n - ip, MAX_REF);
          size_t len = 3;
          while (len < maxlen && in[r + len] == in[ip + len]) ++len;
          if (lit) {
            out[lit_ctrl] = static_cast<uint8_t>(lit - 1);
          } else {
            --op;  // drop the unused literal control byte
          }
          size_t code = len - 2;
          if (op + 4 > cap) return 0;
          if (code < 7) {
            out[op++] = static_cast<uint8_t>((off >> 8) + (code << 5));
          } else {
            out[op++] = static_cast<uint8_t>((off >> 8) + (7 << 5));
            out[op++] = static_cast<uint8_t>(code - 7);
          }
          out[op++] = static_cast<uint8_t>(off & 0xff);
          ip += len;
          lit = 0;
          lit_ctrl = op++;
          continue;
        }
      }
    }
    if (op + 2 > cap) return 0;
    out[op++] = in[ip++];
    if (++lit == MAX_LIT) {
      out[lit_ctrl] = static_cast<uint8_t>(MAX_LIT - 1);
      lit = 0;
      lit_ctrl = op++;
    }
  }
  if (lit) {
    out[lit_ctrl] = static_cast<uint8_t>(lit - 1);
  } else {
    --op;
  }
  return op;
}

// ---------------------------------------------------------------------------
// HDF5 shuffle filter: byte b of element i lives at in[b * n + i].

void unshuffle(const uint8_t* in, uint8_t* out, size_t nbytes, size_t es) {
  if (es <= 1 || nbytes < es) {
    std::memcpy(out, in, nbytes);
    return;
  }
  size_t n = nbytes / es;
  for (size_t b = 0; b < es; ++b) {
    const uint8_t* src = in + b * n;
    uint8_t* dst = out + b;
    for (size_t i = 0; i < n; ++i) dst[i * es] = src[i];
  }
  size_t done = n * es;
  if (nbytes > done) std::memcpy(out + done, in + done, nbytes - done);
}

void shuffle(const uint8_t* in, uint8_t* out, size_t nbytes, size_t es) {
  if (es <= 1 || nbytes < es) {
    std::memcpy(out, in, nbytes);
    return;
  }
  size_t n = nbytes / es;
  for (size_t b = 0; b < es; ++b) {
    uint8_t* dst = out + b * n;
    const uint8_t* src = in + b;
    for (size_t i = 0; i < n; ++i) dst[i] = src[i * es];
  }
  size_t done = n * es;
  if (nbytes > done) std::memcpy(out + done, in + done, nbytes - done);
}

// ---------------------------------------------------------------------------
// Deflate (zlib container, as written by HDF5's filter 1)

bool inflate_zlib(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len, void** state) {
#ifdef ANNR_HAVE_LIBDEFLATE
  if (*state == nullptr) {
    *state = libdeflate_alloc_decompressor();
    if (*state == nullptr) return false;
  }
  size_t actual = 0;
  enum libdeflate_result rc = libdeflate_zlib_decompress(
      static_cast<libdeflate_decompressor*>(*state), in, in_len, out, out_len, &actual);
  if (rc == LIBDEFLATE_SUCCESS && actual == out_len) return true;
  // libdeflate rejects some valid-but-unusual streams (e.g. trailing data);
  // retry with zlib before giving up.
#else
  (void)state;
#endif
  uLongf dl = static_cast<uLongf>(out_len);
  int rc2 = uncompress(out, &dl, in, static_cast<uLong>(in_len));
  return rc2 == Z_OK && dl == out_len;
}

void free_inflate_state(void* state) {
#ifdef ANNR_HAVE_LIBDEFLATE
  if (state != nullptr) libdeflate_free_decompressor(static_cast<libdeflate_decompressor*>(state));
#else
  (void)state;
#endif
}

size_t deflate_bound(size_t in_len) {
#ifdef ANNR_HAVE_LIBDEFLATE
  return libdeflate_zlib_compress_bound(nullptr, in_len);
#else
  return static_cast<size_t>(compressBound(static_cast<uLong>(in_len)));
#endif
}

size_t deflate_zlib(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_cap, int level,
                    void** state) {
#ifdef ANNR_HAVE_LIBDEFLATE
  if (*state == nullptr) {
    *state = libdeflate_alloc_compressor(level);
    if (*state == nullptr) return 0;
  }
  return libdeflate_zlib_compress(static_cast<libdeflate_compressor*>(*state), in, in_len, out,
                                  out_cap);
#else
  (void)state;
  uLongf dl = static_cast<uLongf>(out_cap);
  int rc = compress2(out, &dl, in, static_cast<uLong>(in_len), level);
  return rc == Z_OK ? static_cast<size_t>(dl) : 0;
#endif
}

void free_deflate_state(void* state) {
#ifdef ANNR_HAVE_LIBDEFLATE
  if (state != nullptr) libdeflate_free_compressor(static_cast<libdeflate_compressor*>(state));
#else
  (void)state;
#endif
}

// ---------------------------------------------------------------------------
// Zstandard (HDF5 registered filter 32015: one zstd frame per chunk)

bool zstd_available() {
#ifdef ANNR_HAVE_ZSTD
  return true;
#else
  return false;
#endif
}

bool zstd_decompress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_len) {
#ifdef ANNR_HAVE_ZSTD
  size_t r = ZSTD_decompress(out, out_len, in, in_len);
  return !ZSTD_isError(r) && r == out_len;
#else
  (void)in; (void)in_len; (void)out; (void)out_len;
  return false;
#endif
}

size_t zstd_bound(size_t in_len) {
#ifdef ANNR_HAVE_ZSTD
  return ZSTD_compressBound(in_len);
#else
  return in_len;
#endif
}

size_t zstd_compress(const uint8_t* in, size_t in_len, uint8_t* out, size_t out_cap, int level) {
#ifdef ANNR_HAVE_ZSTD
  size_t r = ZSTD_compress(out, out_cap, in, in_len, level);
  return ZSTD_isError(r) ? 0 : r;
#else
  (void)in; (void)in_len; (void)out; (void)out_cap; (void)level;
  return 0;
#endif
}

// ---------------------------------------------------------------------------
// HDF5 filter callbacks

static size_t lzf_h5_filter(unsigned flags, size_t cd_nelmts, const unsigned cd_values[],
                            size_t nbytes, size_t* buf_size, void** buf) {
  if (flags & H5Z_FLAG_REVERSE) {
    size_t outsize = (cd_nelmts >= 3 && cd_values[2] != 0) ? cd_values[2] : (*buf_size) * 2;
    for (int attempt = 0; attempt < 8; ++attempt) {
      void* outbuf = H5allocate_memory(outsize, false);
      if (outbuf == nullptr) return 0;
      size_t got = lzf_decompress(static_cast<const uint8_t*>(*buf), nbytes,
                                  static_cast<uint8_t*>(outbuf), outsize);
      if (got != 0) {
        H5free_memory(*buf);
        *buf = outbuf;
        *buf_size = outsize;
        return got;
      }
      H5free_memory(outbuf);
      if (cd_nelmts >= 3 && cd_values[2] != 0) return 0;  // size was exact: corrupt
      outsize *= 2;
    }
    return 0;
  }
  if (nbytes < 4) return 0;
  void* outbuf = H5allocate_memory(nbytes, false);
  if (outbuf == nullptr) return 0;
  size_t got = lzf_compress(static_cast<const uint8_t*>(*buf), nbytes, static_cast<uint8_t*>(outbuf),
                            nbytes - 1);
  if (got == 0) {  // incompressible: the optional filter is skipped for this chunk
    H5free_memory(outbuf);
    return 0;
  }
  H5free_memory(*buf);
  *buf = outbuf;
  *buf_size = nbytes;
  return got;
}

#ifdef ANNR_HAVE_ZSTD
static size_t zstd_h5_filter(unsigned flags, size_t cd_nelmts, const unsigned cd_values[],
                             size_t nbytes, size_t* buf_size, void** buf) {
  if (flags & H5Z_FLAG_REVERSE) {
    unsigned long long sz = ZSTD_getFrameContentSize(*buf, nbytes);
    if (sz == ZSTD_CONTENTSIZE_ERROR || sz == ZSTD_CONTENTSIZE_UNKNOWN) return 0;
    void* outbuf = H5allocate_memory(static_cast<size_t>(sz), false);
    if (outbuf == nullptr) return 0;
    size_t r = ZSTD_decompress(outbuf, static_cast<size_t>(sz), *buf, nbytes);
    if (ZSTD_isError(r)) {
      H5free_memory(outbuf);
      return 0;
    }
    H5free_memory(*buf);
    *buf = outbuf;
    *buf_size = static_cast<size_t>(sz);
    return r;
  }
  int level = cd_nelmts >= 1 ? static_cast<int>(cd_values[0]) : 3;
  size_t cap = ZSTD_compressBound(nbytes);
  void* outbuf = H5allocate_memory(cap, false);
  if (outbuf == nullptr) return 0;
  size_t r = ZSTD_compress(outbuf, cap, *buf, nbytes, level);
  if (ZSTD_isError(r)) {
    H5free_memory(outbuf);
    return 0;
  }
  H5free_memory(*buf);
  *buf = outbuf;
  *buf_size = cap;
  return r;
}
#endif

void register_filters() {
  if (H5Zfilter_avail(FILTER_LZF) <= 0) {
    static const H5Z_class2_t lzf_class = {H5Z_CLASS_T_VERS, FILTER_LZF, 1, 1, "lzf (AnnR)",
                                           nullptr, nullptr, lzf_h5_filter};
    H5Zregister(&lzf_class);
  }
#ifdef ANNR_HAVE_ZSTD
  if (H5Zfilter_avail(FILTER_ZSTD) <= 0) {
    static const H5Z_class2_t zstd_class = {H5Z_CLASS_T_VERS, FILTER_ZSTD, 1, 1, "zstd (AnnR)",
                                            nullptr, nullptr, zstd_h5_filter};
    H5Zregister(&zstd_class);
  }
#endif
  H5Eclear2(H5E_DEFAULT);
}

}  // namespace annr
