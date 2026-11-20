#include <climits>
#include <cmath>
#include <cstring>
#include <stdexcept>

#include "sample_data.h"
#include "types.h"

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  #include <immintrin.h>
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
#endif

namespace bcf {

SampleData::SampleData(GzReader & infile, const Header & _header, std::uint32_t len, std::uint32_t _n_fmt, std::uint32_t _n_samples)
    : header(&_header), n_fmt(_n_fmt), buf_len(len), n_samples(_n_samples) {
  if (len > 0) {
    buf.reset(new char[len]);
    infile.read_exact(buf.get(), len, "truncated BCF sample data");
  }
}

void SampleData::parse() {
  if (is_parsed) {
    return;
  }
  if (n_fmt == 0) {
    is_parsed = true;
    return;
  }
  if (header == nullptr) {
    throw std::invalid_argument("missing BCF header for FORMAT");
  }

  const char * raw = buf.get();
  fields.reserve(n_fmt);
  std::uint32_t buf_idx = 0;
  for (std::uint32_t i = 0; i < n_fmt; i++) {
    std::int32_t format_idx = parse_typed_int(raw, buf_len, buf_idx);
    if (format_idx < 0 || static_cast<std::size_t>(format_idx) >= header->format.size() ||
        header->format[format_idx].id.empty()) {
      throw std::invalid_argument("unknown FORMAT index: " + std::to_string(format_idx));
    }
    std::uint32_t uid = static_cast<std::uint32_t>(format_idx);
    const FormatField * fmt_field = &header->format[uid];
    bool is_geno = fmt_field->is_gt;

    Typed type_val(raw, buf_len, buf_idx);
    if (type_val.type == flag && type_val.n_vals != 0) {
      throw std::invalid_argument("invalid BCF FORMAT flag length");
    }
    if (is_geno && !is_int_type(type_val.type)) {
      throw std::invalid_argument("GT field must be encoded as integer");
    }
    std::uint64_t field_bytes = static_cast<std::uint64_t>(type_val.n_vals) *
                                static_cast<std::uint64_t>(type_val.type_size) *
                                static_cast<std::uint64_t>(n_samples);
    if (static_cast<std::uint64_t>(buf_idx) + field_bytes > buf_len) {
      throw std::invalid_argument("FORMAT field exceeds BCF sample buffer");
    }
    FormatType fmt = {type_val.type, type_val.type_size, buf_idx, type_val.n_vals, is_geno};
    bool updated = false;
    for (auto & entry : fields) {
      if (entry.id_idx == uid) {
        entry.format = fmt;
        updated = true;
        break;
      }
    }
    if (!updated) {
      fields.push_back({uid, fmt_field, fmt});
    }
    buf_idx += static_cast<std::uint32_t>(field_bytes);
  }
  is_parsed = true;
}

std::size_t SampleData::size() {
  parse();
  return fields.size();
}

std::vector<NamedFormat> & SampleData::get_fields() {
  parse();
  return fields;
}

bool SampleData::get_type(const char * key, std::size_t key_len, FormatType & out) {
  parse();
  for (const auto & entry : fields) {
    const std::string & id = entry.field->id;
    if (id.size() == key_len && std::memcmp(id.data(), key, key_len) == 0) {
      out = entry.format;
      return true;
    }
  }
  return false;
}

void SampleData::read_ints(const FormatType & type, double * dst) const {
  if (type.is_geno) {
    read_geno(type, dst);
  } else if (type.type_size == 1) {
    read_8bit_ints(type, dst);
  } else if (type.type_size == 2) {
    read_16bit_ints(type, dst);
  } else if (type.type_size == 4) {
    read_32bit_ints(type, dst);
  } else {
    throw std::invalid_argument("unsupported FORMAT integer byte width: " + std::to_string(type.type_size));
  }
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
static inline bool cpu_has_avx2() {
  static const bool supported = __builtin_cpu_supports("avx2");
  return supported;
}

// Convert 4 signed 32-bit integers to float64, replacing BCF missing/end-of-vector
// sentinels (values <= miss_thresh, e.g. -128/-127 for int8) with NaN.
__attribute__((target("avx2")))
static inline void store_epi32_as_double(__m128i vals, double * dst, __m256d miss_thresh, __m256d qnan) {
  __m256d d = _mm256_cvtepi32_pd(vals);
  __m256d mask = _mm256_cmp_pd(d, miss_thresh, _CMP_LE_OQ);
  _mm256_storeu_pd(dst, _mm256_blendv_pd(d, qnan, mask));
}

__attribute__((target("avx2")))
static std::uint64_t parse_8bit_ints_avx2(const char * src, double * dst, std::uint64_t max_n) {
  const __m256d miss_thresh = _mm256_set1_pd(-127.0);
  const __m256d qnan = _mm256_set1_pd(std::nan(""));

  std::uint64_t n = 0;
  for (; n + 32 <= max_n; n += 32) {
    __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n));
    __m128i low = _mm256_castsi256_si128(data);
    __m128i hi = _mm256_extracti128_si256(data, 1);

    store_epi32_as_double(_mm_cvtepi8_epi32(low), dst + n, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)), dst + n + 4, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)), dst + n + 8, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)), dst + n + 12, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(hi), dst + n + 16, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)), dst + n + 20, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)), dst + n + 24, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)), dst + n + 28, miss_thresh, qnan);
  }
  return n;
}

__attribute__((target("avx2")))
static std::uint64_t parse_16bit_ints_avx2(const char * src, double * dst, std::uint64_t max_n) {
  const __m256d miss_thresh = _mm256_set1_pd(-32767.0);
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

__attribute__((target("avx2")))
static std::uint64_t parse_32bit_ints_avx2(const char * src, double * dst, std::uint64_t max_n) {
  const __m256d miss_thresh = _mm256_set1_pd(static_cast<double>(INT32_MIN + 1));
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

// BCF encodes each GT allele as ((allele + 1) << 1) | phased, with 0 for missing ('.')
// and negative values (0x80/0x81) for missing/end-of-vector padding. Shift right by 1,
// subtract 1, and force negative bytes to -1 so both missing encodings map to NaN.
__attribute__((target("avx2")))
static std::uint64_t parse_geno_avx2(const char * src, double * dst, std::uint64_t max_n) {
  const __m256i mask_geno = _mm256_set1_epi8(static_cast<char>(0xfe));
  const __m256i sub = _mm256_set1_epi8(1);
  const __m256i zero = _mm256_setzero_si256();
  const __m256d miss_thresh = _mm256_set1_pd(-1.0);
  const __m256d qnan = _mm256_set1_pd(std::nan(""));

  std::uint64_t n = 0;
  for (; n + 32 <= max_n; n += 32) {
    __m256i initial = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n));
    __m256i geno = _mm256_and_si256(initial, mask_geno);
    geno = _mm256_sub_epi8(_mm256_srli_epi32(geno, 1), sub);
    geno = _mm256_or_si256(geno, _mm256_cmpgt_epi8(zero, initial));

    __m128i low = _mm256_castsi256_si128(geno);
    __m128i hi = _mm256_extracti128_si256(geno, 1);

    store_epi32_as_double(_mm_cvtepi8_epi32(low), dst + n, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)), dst + n + 4, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)), dst + n + 8, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)), dst + n + 12, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(hi), dst + n + 16, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)), dst + n + 20, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)), dst + n + 24, miss_thresh, qnan);
    store_epi32_as_double(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)), dst + n + 28, miss_thresh, qnan);
  }
  return n;
}

// For diploid 8-bit GT (2 bytes per sample), phasing is stored in the low bit of
// the second allele byte (odd byte offsets 1, 3, 5, ...).
__attribute__((target("avx2")))
static std::uint32_t parse_phase_avx2(const char * src, std::uint8_t * dst, std::uint32_t n_samples) {
  const __m128i one = _mm_set1_epi8(1);
  const __m128i odd_lo = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);
  const __m128i odd_hi = _mm_setr_epi8(-1, -1, -1, -1, -1, -1, -1, -1, 1, 3, 5, 7, 9, 11, 13, 15);

  std::uint32_t s = 0;
  for (; s + 16 <= n_samples; s += 16) {
    __m128i lo = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + s * 2));
    __m128i hi = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + s * 2 + 16));
    __m128i second = _mm_or_si128(_mm_shuffle_epi8(lo, odd_lo), _mm_shuffle_epi8(hi, odd_hi));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(dst + s), _mm_and_si128(second, one));
  }
  return s;
}

__attribute__((target("avx2")))
static std::uint64_t parse_floats_avx2(const char * src, float * dst, std::uint64_t max_n) {
  const __m256i miss1 = _mm256_set1_epi32(static_cast<std::int32_t>(BCF_FLOAT_MISSING));
  const __m256i miss2 = _mm256_set1_epi32(static_cast<std::int32_t>(BCF_FLOAT_VECTOR_END));
  const __m256i qnan = _mm256_set1_epi32(0x7fc00000);

  std::uint64_t n = 0;
  for (; n + 8 <= max_n; n += 8) {
    __m256i bits = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n * 4));
    __m256i is_miss = _mm256_or_si256(_mm256_cmpeq_epi32(bits, miss1),
                                      _mm256_cmpeq_epi32(bits, miss2));
    bits = _mm256_blendv_epi8(bits, qnan, is_miss);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(dst + n), bits);
  }
  return n;
}
#elif defined(__aarch64__)
static inline void store_s32_as_double_neon(int32x4_t vals, double * dst, float64x2_t miss_thresh, float64x2_t qnan) {
  float64x2_t d_lo = vcvtq_f64_s64(vmovl_s32(vget_low_s32(vals)));
  float64x2_t d_hi = vcvtq_f64_s64(vmovl_s32(vget_high_s32(vals)));
  vst1q_f64(dst, vbslq_f64(vcleq_f64(d_lo, miss_thresh), qnan, d_lo));
  vst1q_f64(dst + 2, vbslq_f64(vcleq_f64(d_hi, miss_thresh), qnan, d_hi));
}
#endif

void SampleData::read_8bit_ints(const FormatType & type, double * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n > buf_len) {
    throw std::invalid_argument("FORMAT int8 data exceeds buffer");
  }
  const char * raw = buf.get();
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (cpu_has_avx2()) {
    n = parse_8bit_ints_avx2(raw + offset, dst, max_n);
  }
#elif defined(__aarch64__)
  const float64x2_t miss_thresh = vdupq_n_f64(-127.0);
  const float64x2_t qnan_vec = vdupq_n_f64(std::nan(""));
  for (; n + 16 <= max_n; n += 16) {
    int8x16_t data = vld1q_s8(reinterpret_cast<const std::int8_t *>(raw + offset + n));
    int16x8_t wider = vmovl_s8(vget_low_s8(data));
    store_s32_as_double_neon(vmovl_s16(vget_low_s16(wider)), dst + n, miss_thresh, qnan_vec);
    store_s32_as_double_neon(vmovl_s16(vget_high_s16(wider)), dst + n + 4, miss_thresh, qnan_vec);
    wider = vmovl_s8(vget_high_s8(data));
    store_s32_as_double_neon(vmovl_s16(vget_low_s16(wider)), dst + n + 8, miss_thresh, qnan_vec);
    store_s32_as_double_neon(vmovl_s16(vget_high_s16(wider)), dst + n + 12, miss_thresh, qnan_vec);
  }
#endif
  offset += static_cast<std::uint32_t>(n);
  const double qnan = std::nan("");
  for (; n < max_n; n++) {
    std::int32_t v = parse_int(raw, buf_len, offset, 1);
    dst[n] = (v == INT32_MIN) ? qnan : static_cast<double>(v);
  }
}

void SampleData::read_16bit_ints(const FormatType & type, double * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * 2 > buf_len) {
    throw std::invalid_argument("FORMAT int16 data exceeds buffer");
  }
  const char * raw = buf.get();
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (cpu_has_avx2()) {
    n = parse_16bit_ints_avx2(raw + offset, dst, max_n);
  }
#elif defined(__aarch64__)
  const float64x2_t miss_thresh = vdupq_n_f64(-32767.0);
  const float64x2_t qnan_vec = vdupq_n_f64(std::nan(""));
  for (; n + 8 <= max_n; n += 8) {
    int16x8_t data = vld1q_s16(reinterpret_cast<const std::int16_t *>(raw + offset + n * 2));
    store_s32_as_double_neon(vmovl_s16(vget_low_s16(data)), dst + n, miss_thresh, qnan_vec);
    store_s32_as_double_neon(vmovl_s16(vget_high_s16(data)), dst + n + 4, miss_thresh, qnan_vec);
  }
#endif
  offset += static_cast<std::uint32_t>(n * 2);
  const double qnan = std::nan("");
  for (; n < max_n; n++) {
    std::int32_t v = parse_int(raw, buf_len, offset, 2);
    dst[n] = (v == INT32_MIN) ? qnan : static_cast<double>(v);
  }
}

void SampleData::read_32bit_ints(const FormatType & type, double * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * 4 > buf_len) {
    throw std::invalid_argument("FORMAT int32 data exceeds buffer");
  }
  const char * raw = buf.get();
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (cpu_has_avx2()) {
    n = parse_32bit_ints_avx2(raw + offset, dst, max_n);
  }
#elif defined(__aarch64__)
  const float64x2_t miss_thresh = vdupq_n_f64(static_cast<double>(INT32_MIN + 1));
  const float64x2_t qnan_vec = vdupq_n_f64(std::nan(""));
  for (; n + 4 <= max_n; n += 4) {
    int32x4_t data = vld1q_s32(reinterpret_cast<const std::int32_t *>(raw + offset + n * 4));
    store_s32_as_double_neon(data, dst + n, miss_thresh, qnan_vec);
  }
#endif
  offset += static_cast<std::uint32_t>(n * 4);
  const double qnan = std::nan("");
  for (; n < max_n; n++) {
    std::int32_t v = parse_int(raw, buf_len, offset, 4);
    dst[n] = (v == INT32_MIN) ? qnan : static_cast<double>(v);
  }
}

void SampleData::read_geno(const FormatType & type, double * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * type.type_size > buf_len) {
    throw std::invalid_argument("FORMAT GT data exceeds buffer");
  }
  const char * raw = buf.get();
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (cpu_has_avx2() && type.type_size == 1) {
    n = parse_geno_avx2(raw + offset, dst, max_n);
  }
#elif defined(__aarch64__)
  if (type.type_size == 1) {
    const uint8x16_t sub = vdupq_n_u8(1);
    const uint8x16_t miss_u8 = vdupq_n_u8(0x80);
    const float64x2_t miss_thresh = vdupq_n_f64(-1.0);
    const float64x2_t qnan_vec = vdupq_n_f64(std::nan(""));

    for (; n + 16 <= max_n; n += 16) {
      uint8x16_t initial = vld1q_u8(reinterpret_cast<const std::uint8_t *>(raw + offset + n));
      int8x16_t geno = vreinterpretq_s8_u8(vsubq_u8(vshrq_n_u8(initial, 1), sub));
      geno = vorrq_s8(geno, vreinterpretq_s8_u8(vcgeq_u8(initial, miss_u8)));

      int16x8_t wider = vmovl_s8(vget_low_s8(geno));
      store_s32_as_double_neon(vmovl_s16(vget_low_s16(wider)), dst + n, miss_thresh, qnan_vec);
      store_s32_as_double_neon(vmovl_s16(vget_high_s16(wider)), dst + n + 4, miss_thresh, qnan_vec);
      wider = vmovl_s8(vget_high_s8(geno));
      store_s32_as_double_neon(vmovl_s16(vget_low_s16(wider)), dst + n + 8, miss_thresh, qnan_vec);
      store_s32_as_double_neon(vmovl_s16(vget_high_s16(wider)), dst + n + 12, miss_thresh, qnan_vec);
    }
  }
#endif
  offset += static_cast<std::uint32_t>(n * type.type_size);
  const double qnan = std::nan("");
  for (; n < max_n; n++) {
    std::int32_t v = parse_int(raw, buf_len, offset, type.type_size);
    std::int32_t allele = (v < 0) ? -1 : ((v >> 1) - 1);
    dst[n] = (allele < 0) ? qnan : static_cast<double>(allele);
  }
}

void SampleData::read_phase(const FormatType & type, std::uint8_t * dst) const {
  if (n_samples == 0) {
    return;
  }
  if (type.n_vals == 0) {
    std::memset(dst, 0, n_samples);
    return;
  }
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  if (static_cast<std::uint64_t>(type.offset) + max_n * type.type_size > buf_len) {
    throw std::invalid_argument("FORMAT GT data exceeds buffer");
  }
  const char * raw = buf.get();
  std::uint32_t s = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (cpu_has_avx2() && type.n_vals == 2 && type.type_size == 1) {
    s = parse_phase_avx2(raw + type.offset, dst, n_samples);
  }
#elif defined(__aarch64__)
  if (type.n_vals == 2 && type.type_size == 1) {
    const uint8x8_t one = vdup_n_u8(1);
    for (; s + 8 <= n_samples; s += 8) {
      uint8x16_t initial = vld1q_u8(reinterpret_cast<const std::uint8_t *>(raw + type.offset + s * 2));
      uint8x8_t second = vshrn_n_u16(vreinterpretq_u16_u8(initial), 8);
      vst1_u8(dst + s, vand_u8(second, one));
    }
  }
#endif
  std::uint32_t offset = type.offset + s * type.n_vals * type.type_size;
  for (; s < n_samples; s++) {
    std::uint8_t is_phased = 1;
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      std::int32_t v = parse_int(raw, buf_len, offset, type.type_size);
      if (v >= 0 && i > 0) {
        is_phased &= static_cast<std::uint8_t>(v & 1);
      }
    }
    dst[s] = is_phased;
  }
}

void SampleData::read_floats(const FormatType & type, float * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * sizeof(float) > buf_len) {
    throw std::invalid_argument("FORMAT float data exceeds buffer");
  }
  const char * src = buf.get() + type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (cpu_has_avx2()) {
    n = parse_floats_avx2(src, dst, max_n);
  }
#elif defined(__aarch64__)
  const uint32x4_t miss1 = vdupq_n_u32(BCF_FLOAT_MISSING);
  const uint32x4_t miss2 = vdupq_n_u32(BCF_FLOAT_VECTOR_END);
  const uint32x4_t qnan_bits = vdupq_n_u32(0x7fc00000u);
  for (; n + 4 <= max_n; n += 4) {
    uint32x4_t bits = vld1q_u32(reinterpret_cast<const std::uint32_t *>(src + n * 4));
    uint32x4_t is_miss = vorrq_u32(vceqq_u32(bits, miss1), vceqq_u32(bits, miss2));
    vst1q_u32(reinterpret_cast<std::uint32_t *>(dst + n), vbslq_u32(is_miss, qnan_bits, bits));
  }
#endif
  const float qnan = std::nanf("");
  for (; n < max_n; n++) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, src + n * sizeof(float), sizeof(bits));
    if (is_bcf_float_missing(bits)) {
      dst[n] = qnan;
    } else {
      std::memcpy(dst + n, &bits, sizeof(float));
    }
  }
}

std::vector<std::string> SampleData::get_strings(const FormatType & type) const {
  std::vector<std::string> vals(n_samples);
  if (n_samples == 0 || type.n_vals == 0) {
    return vals;
  }
  std::uint32_t offset = type.offset;
  for (std::uint32_t n = 0; n < n_samples; n++) {
    vals[n] = parse_string(buf.get(), buf_len, offset, type.n_vals);
  }
  return vals;
}

} // namespace bcf
