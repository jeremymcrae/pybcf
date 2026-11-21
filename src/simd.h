#ifndef BCF_SIMD_H
#define BCF_SIMD_H

// SIMD kernels for decoding the per-sample FORMAT fields of a BCF record. Each one decodes
// as many leading values as fit its vector width and returns how many it handled; the
// caller finishes the rest with a scalar loop. x86_64 and ARM64 define the same kernel
// names, so the callers in `sample_data.cpp` have no per-architecture code.

#include <cmath>
#include <cstdint>

#include "types.h"

// Choose the SIMD instruction set once here. The kernels further down are written for
// each one, and `BCF_HAVE_SIMD` tells the decoders that kernels exist.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  #include <immintrin.h>
  #define BCF_SIMD_AVX2 1
  #define BCF_HAVE_SIMD 1
#elif defined(__aarch64__)
  #include <arm_neon.h>
  #define BCF_SIMD_NEON 1
  #define BCF_HAVE_SIMD 1
#endif

namespace bcf {

#if defined(BCF_SIMD_AVX2)
// Cache the runtime CPU check once so x86_64 builds remain portable to pre-AVX2 CPUs
// while automatically using 256-bit SIMD on modern x86_64 machines.
static inline bool simd_supported() {
  static const bool supported = __builtin_cpu_supports("avx2");
  return supported;
}

// Convert a 4-lane vector of signed 32-bit integers into 4 float64 values in `dst`.
// BCF reserves the lowest 8 signed integers of each width (e.g. -128..-121 for int8)
// for missing, end-of-vector, and future sentinels, so any value <= `miss_thresh`
// is replaced with quiet NaN.
__attribute__((target("avx2")))
static inline void store_epi32_as_double(__m128i vals, double * dst, __m256d miss_thresh, __m256d qnan) {
  __m256d d = _mm256_cvtepi32_pd(vals);
  __m256d mask = _mm256_cmp_pd(d, miss_thresh, _CMP_LE_OQ);
  _mm256_storeu_pd(dst, _mm256_blendv_pd(d, qnan, mask));
}

// Widen 32 signed int8 lanes to float64: split the 32 bytes into two 16-byte halves,
// shift each half right in 4-byte steps, and sign-extend the bottom 4 int8s to int32
// (`_mm_cvtepi8_epi32`) before converting to float64.
__attribute__((target("avx2")))
static inline void store_epi8_as_double(__m256i data, double * dst, __m256d miss_thresh, __m256d qnan) {
  __m128i low = _mm256_castsi256_si128(data);
  __m128i hi = _mm256_extracti128_si256(data, 1);
  store_epi32_as_double(_mm_cvtepi8_epi32(low), dst, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)), dst + 4, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)), dst + 8, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)), dst + 12, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(hi), dst + 16, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)), dst + 20, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)), dst + 24, miss_thresh, qnan);
  store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)), dst + 28, miss_thresh, qnan);
}

// Decode 32 signed int8 values per iteration.
__attribute__((target("avx2")))
static std::uint64_t simd_int8(const char * src, double * dst, std::uint64_t max_n) {
  const __m256d miss_thresh = _mm256_set1_pd(static_cast<double>(INT8_MIN + 7));
  const __m256d qnan = _mm256_set1_pd(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 32 <= max_n; n += 32) {
    __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n));
    store_epi8_as_double(data, dst + n, miss_thresh, qnan);
  }
  return n;
}

// Decode 16 signed int16 values (32 bytes) per iteration: each 16-byte half holds 8 int16s;
// sign-extend the lower 4 and upper 4 (shifted right by 8 bytes) to int32, then to float64.
__attribute__((target("avx2")))
static std::uint64_t simd_int16(const char * src, double * dst, std::uint64_t max_n) {
  const __m256d miss_thresh = _mm256_set1_pd(static_cast<double>(INT16_MIN + 7));
  const __m256d qnan = _mm256_set1_pd(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 16 <= max_n; n += 16) {
    __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n * 2));
    __m128i low = _mm256_castsi256_si128(data);
    __m128i hi = _mm256_extracti128_si256(data, 1);
    store_epi32_as_double(_mm_cvtepi16_epi32(low), dst + n, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi16_epi32(_mm_bsrli_si128(low, 8)), dst + n + 4, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi16_epi32(hi), dst + n + 8, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi16_epi32(_mm_bsrli_si128(hi, 8)), dst + n + 12, miss_thresh, qnan);
  }
  return n;
}

// Decode 8 signed int32 values (two 4-lane 128-bit loads) per iteration.
__attribute__((target("avx2")))
static std::uint64_t simd_int32(const char * src, double * dst, std::uint64_t max_n) {
  const __m256d miss_thresh = _mm256_set1_pd(static_cast<double>(INT32_MIN + 7));
  const __m256d qnan = _mm256_set1_pd(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 8 <= max_n; n += 8) {
    __m128i low = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + n * 4));
    __m128i hi = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + n * 4 + 16));
    store_epi32_as_double(low, dst + n, miss_thresh, qnan);
    store_epi32_as_double(hi, dst + n + 4, miss_thresh, qnan);
  }
  return n;
}

// Decode 32 8-bit GT allele entries per iteration into float64 allele indices:
// - BCF encodes allele `a` as `((a + 1) << 1) | phased`, missing `.` as 0 (or 1 if phased),
//   and lower-ploidy padding as negative sentinel bytes (0x80..0x87).
// - AVX2 has no per-byte right shift, so clear bit 0 (`0xfe`) first to prevent the low bit
//   of byte `k + 1` from shifting into bit 7 of byte `k` during `_mm256_srli_epi32(..., 1)`.
// - Subtracting 1 turns `(0 >> 1) - 1` into `-1` for `.` and gives `0, 1, 2, ...` for alleles.
// - `_mm256_cmpgt_epi8(zero, initial)` produces `0xff` (`-1`) for any negative sentinel byte,
//   forcing it to `-1` as well so every missing/padding lane becomes NaN.
__attribute__((target("avx2")))
static std::uint64_t simd_gt(const char * src, double * dst, std::uint64_t max_n) {
  const __m256i mask_gt = _mm256_set1_epi8(static_cast<char>(0xfe));
  const __m256i sub = _mm256_set1_epi8(1);
  const __m256i zero = _mm256_setzero_si256();
  const __m256d miss_thresh = _mm256_set1_pd(-1.0);
  const __m256d qnan = _mm256_set1_pd(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 32 <= max_n; n += 32) {
    __m256i initial = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n));
    __m256i gt = _mm256_and_si256(initial, mask_gt);
    gt = _mm256_sub_epi8(_mm256_srli_epi32(gt, 1), sub);
    gt = _mm256_or_si256(gt, _mm256_cmpgt_epi8(zero, initial));
    store_epi8_as_double(gt, dst + n, miss_thresh, qnan);
  }
  return n;
}

// Extract 16 diploid phase flags (0 or 1) per iteration from 32 bytes of 8-bit GT data:
// - Each sample occupies 2 bytes (`allele0`, `allele1`), and phasing is stored in bit 0
//   of `allele1` (the odd-indexed bytes 1, 3, 5, ..., 15 in each 16-byte load).
// - `_mm_shuffle_epi8` gathers the 8 odd bytes from `lo` into lanes 0..7 and the 8 odd
//   bytes from `hi` into lanes 8..15 (lanes with shuffle mask `-1` are zeroed).
// - Negative `allele1` bytes (e.g. `0x81` end-of-vector padding on haploid samples) set
//   `0xff` via `_mm_cmpgt_epi8(zero, second)` so haploid samples report `phased = 1`.
__attribute__((target("avx2")))
static std::uint32_t simd_phase(const char * src, std::uint8_t * dst, std::uint32_t n_samples) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi8(1);
  const __m128i odd_lo = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);
  const __m128i odd_hi = _mm_setr_epi8(-1, -1, -1, -1, -1, -1, -1, -1, 1, 3, 5, 7, 9, 11, 13, 15);
  std::uint32_t s = 0;
  for (; s + 16 <= n_samples; s += 16) {
    __m128i lo = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + s * 2));
    __m128i hi = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + s * 2 + 16));
    __m128i second = _mm_or_si128(_mm_shuffle_epi8(lo, odd_lo), _mm_shuffle_epi8(hi, odd_hi));
    second = _mm_or_si128(second, _mm_cmpgt_epi8(zero, second));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(dst + s), _mm_and_si128(second, one));
  }
  return s;
}

// Copy 8 float32 values per iteration while replacing BCF's signaling-NaN sentinels
// (`0x7f800001..0x7f800007`) with standard IEEE quiet NaN (`0x7fc00000`). Comparing raw
// 32-bit integer patterns avoids touching signaling NaNs in floating-point registers.
__attribute__((target("avx2")))
static std::uint64_t simd_floats(const char * src, float * dst, std::uint64_t max_n) {
  const __m256i miss_lo = _mm256_set1_epi32(static_cast<std::int32_t>(BCF_FLOAT_MISSING - 1));
  const __m256i miss_hi = _mm256_set1_epi32(static_cast<std::int32_t>(BCF_FLOAT_RESERVED_MAX + 1));
  const __m256i qnan = _mm256_set1_epi32(0x7fc00000);
  std::uint64_t n = 0;
  for (; n + 8 <= max_n; n += 8) {
    __m256i bits = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n * 4));
    __m256i is_miss = _mm256_and_si256(_mm256_cmpgt_epi32(bits, miss_lo),
                                       _mm256_cmpgt_epi32(miss_hi, bits));
    bits = _mm256_blendv_epi8(bits, qnan, is_miss);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(dst + n), bits);
  }
  return n;
}
#elif defined(BCF_SIMD_NEON)
// NEON is part of the ARM64 baseline, so no runtime check is needed.
static inline bool simd_supported() {
  return true;
}

// Widen 4 signed int32s to two pairs of float64s and replace BCF missing/reserved
// sentinels (<= `miss_thresh`) with quiet NaN.
static inline void store_s32_as_double(int32x4_t vals, double * dst, float64x2_t miss_thresh,
                                       float64x2_t qnan) {
  float64x2_t d_lo = vcvtq_f64_s64(vmovl_s32(vget_low_s32(vals)));
  float64x2_t d_hi = vcvtq_f64_s64(vmovl_s32(vget_high_s32(vals)));
  vst1q_f64(dst, vbslq_f64(vcleq_f64(d_lo, miss_thresh), qnan, d_lo));
  vst1q_f64(dst + 2, vbslq_f64(vcleq_f64(d_hi, miss_thresh), qnan, d_hi));
}

// Widen 16 int8s -> two 8-lane int16 vectors -> four 4-lane int32 vectors -> float64.
static inline void store_s8_as_double(int8x16_t data, double * dst, float64x2_t miss_thresh,
                                      float64x2_t qnan) {
  int16x8_t wider = vmovl_s8(vget_low_s8(data));
  store_s32_as_double(vmovl_s16(vget_low_s16(wider)), dst, miss_thresh, qnan);
  store_s32_as_double(vmovl_s16(vget_high_s16(wider)), dst + 4, miss_thresh, qnan);
  wider = vmovl_s8(vget_high_s8(data));
  store_s32_as_double(vmovl_s16(vget_low_s16(wider)), dst + 8, miss_thresh, qnan);
  store_s32_as_double(vmovl_s16(vget_high_s16(wider)), dst + 12, miss_thresh, qnan);
}

static std::uint64_t simd_int8(const char * src, double * dst, std::uint64_t max_n) {
  const float64x2_t miss_thresh = vdupq_n_f64(static_cast<double>(INT8_MIN + 7));
  const float64x2_t qnan = vdupq_n_f64(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 16 <= max_n; n += 16) {
    int8x16_t data = vld1q_s8(reinterpret_cast<const std::int8_t *>(src + n));
    store_s8_as_double(data, dst + n, miss_thresh, qnan);
  }
  return n;
}

static std::uint64_t simd_int16(const char * src, double * dst, std::uint64_t max_n) {
  const float64x2_t miss_thresh = vdupq_n_f64(static_cast<double>(INT16_MIN + 7));
  const float64x2_t qnan = vdupq_n_f64(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 8 <= max_n; n += 8) {
    int16x8_t data = vld1q_s16(reinterpret_cast<const std::int16_t *>(src + n * 2));
    store_s32_as_double(vmovl_s16(vget_low_s16(data)), dst + n, miss_thresh, qnan);
    store_s32_as_double(vmovl_s16(vget_high_s16(data)), dst + n + 4, miss_thresh, qnan);
  }
  return n;
}

static std::uint64_t simd_int32(const char * src, double * dst, std::uint64_t max_n) {
  const float64x2_t miss_thresh = vdupq_n_f64(static_cast<double>(INT32_MIN + 7));
  const float64x2_t qnan = vdupq_n_f64(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 4 <= max_n; n += 4) {
    int32x4_t data = vld1q_s32(reinterpret_cast<const std::int32_t *>(src + n * 4));
    store_s32_as_double(data, dst + n, miss_thresh, qnan);
  }
  return n;
}

// NEON has a true unsigned per-byte shift (`vshrq_n_u8`), so `(initial >> 1) - 1`
// directly yields allele indices while `initial >= 0x80` masks negative sentinels to -1.
static std::uint64_t simd_gt(const char * src, double * dst, std::uint64_t max_n) {
  const uint8x16_t sub = vdupq_n_u8(1);
  const uint8x16_t miss_u8 = vdupq_n_u8(0x80);
  const float64x2_t miss_thresh = vdupq_n_f64(-1.0);
  const float64x2_t qnan = vdupq_n_f64(std::nan(""));
  std::uint64_t n = 0;
  for (; n + 16 <= max_n; n += 16) {
    uint8x16_t initial = vld1q_u8(reinterpret_cast<const std::uint8_t *>(src + n));
    int8x16_t gt = vreinterpretq_s8_u8(vsubq_u8(vshrq_n_u8(initial, 1), sub));
    gt = vorrq_s8(gt, vreinterpretq_s8_u8(vcgeq_u8(initial, miss_u8)));
    store_s8_as_double(gt, dst + n, miss_thresh, qnan);
  }
  return n;
}

// Reinterpret 8 diploid pairs as uint16 lanes and shift right by 8 bits to isolate
// each sample's second byte; OR with bit 7 so negative sentinels map to phased (1).
static std::uint32_t simd_phase(const char * src, std::uint8_t * dst, std::uint32_t n_samples) {
  const uint8x8_t one = vdup_n_u8(1);
  std::uint32_t s = 0;
  for (; s + 8 <= n_samples; s += 8) {
    uint8x16_t initial = vld1q_u8(reinterpret_cast<const std::uint8_t *>(src + s * 2));
    uint8x8_t second = vshrn_n_u16(vreinterpretq_u16_u8(initial), 8);
    second = vorr_u8(second, vshr_n_u8(second, 7));
    vst1_u8(dst + s, vand_u8(second, one));
  }
  return s;
}

static std::uint64_t simd_floats(const char * src, float * dst, std::uint64_t max_n) {
  const uint32x4_t miss_lo = vdupq_n_u32(BCF_FLOAT_MISSING);
  const uint32x4_t miss_hi = vdupq_n_u32(BCF_FLOAT_RESERVED_MAX);
  const uint32x4_t qnan_bits = vdupq_n_u32(0x7fc00000u);
  std::uint64_t n = 0;
  for (; n + 4 <= max_n; n += 4) {
    uint32x4_t bits = vld1q_u32(reinterpret_cast<const std::uint32_t *>(src + n * 4));
    uint32x4_t is_miss = vandq_u32(vcgeq_u32(bits, miss_lo), vcleq_u32(bits, miss_hi));
    vst1q_u32(reinterpret_cast<std::uint32_t *>(dst + n), vbslq_u32(is_miss, qnan_bits, bits));
  }
  return n;
}
#endif

} // namespace bcf

#endif
