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

SampleData::SampleData(GzReader & infile, const Header & _header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples) {
  n_samples = _n_samples;
  if (len == 0) {
    return;
  }

  // read the sample data into a buffer, but don't parse until required
  buf.resize(len);
  infile.read_exact(buf.data(), len, "truncated BCF sample data");

  // read the available keys
  std::uint32_t buf_idx = 0;
  Typed type_val;
  for (std::uint32_t i = 0; i < n_fmt; i++) {
    type_val = Typed(buf.data(), buf.size(), buf_idx);
    std::int32_t format_idx = parse_int(buf.data(), buf.size(), buf_idx, type_val.type_size);
    auto fmt_it = (format_idx >= 0) ? _header.format.find(static_cast<std::uint32_t>(format_idx))
                                    : _header.format.end();
    if (fmt_it == _header.format.end()) {
      throw std::invalid_argument("unknown FORMAT index: " + std::to_string(format_idx));
    }
    const std::string & key = fmt_it->second;
    bool is_geno = key == "GT";

    type_val = Typed(buf.data(), buf.size(), buf_idx);
    std::uint64_t field_bytes = static_cast<std::uint64_t>(type_val.n_vals) *
                                static_cast<std::uint64_t>(type_val.type_size) *
                                static_cast<std::uint64_t>(n_samples);
    if (static_cast<std::uint64_t>(buf_idx) + field_bytes > buf.size()) {
      throw std::invalid_argument("FORMAT field exceeds BCF sample buffer");
    }
    keys[key] = {type_val.type, type_val.type_size, buf_idx, type_val.n_vals, is_geno};
    buf_idx += static_cast<std::uint32_t>(field_bytes);
  }
}

/// @brief get key names as strings
/// @return vector of sample format keys
std::vector<std::string> SampleData::get_keys() {
  std::vector<std::string> key_names;
  for (auto & x : keys) {
    key_names.push_back(x.first);
  }
  return key_names;
}

/// @brief look up the data type for a given format key
/// @param key string for the format key to check
/// @param out struct populated with type, number of values, and buffer offset when found
/// @return true if key is present in this variant's sample data
bool SampleData::get_type(const std::string & key, FormatType & out) const {
  auto it = keys.find(key);
  if (it == keys.end()) {
    return false;
  }
  out = it->second;
  return true;
}

/// @brief get vector of ints for a given format key
/// @param type struct containing the type size, number of values etc
/// @return vector of 32-bit ints
std::vector<std::int32_t> SampleData::get_ints(FormatType & type) {
  if (type.is_geno) {
    return get_geno(type);
  } else if (type.type_size == 1) {
    return parse_8bit_ints(type);
  } else if (type.type_size == 2) {
    return parse_16bit_ints(type);
  } else if (type.type_size == 4) {
    return parse_32bit_ints(type);
  }
  throw std::invalid_argument("unsupported FORMAT integer byte width: " + std::to_string(type.type_size));
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
/// @brief replace sign-extended narrow missing sentinel with INT32_MIN
__attribute__((target("avx2")))
static inline __m128i replace_missing_epi32(__m128i data, std::int32_t old_missing) {
  const __m128i missing_old = _mm_set1_epi32(old_missing);
  const __m128i missing_32bit = _mm_set1_epi32(INT32_MIN);
  __m128i mask = _mm_cmpeq_epi32(data, missing_old);
  return _mm_blendv_epi8(data, missing_32bit, mask);
}

__attribute__((target("avx2")))
static std::uint64_t parse_8bit_ints_avx2(const char * src, std::int32_t * dst, std::uint64_t max_n) {
  const __m256i missing = _mm256_set1_epi8(static_cast<char>(0x80));
  const __m256i not_recorded = _mm256_set1_epi8(static_cast<char>(0x81));

  std::uint64_t n = 0;
  for (; n + 32 <= max_n; n += 32) {
    __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n));
    data = _mm256_blendv_epi8(data, missing, _mm256_cmpeq_epi8(data, not_recorded));

    __m128i low = _mm256_castsi256_si128(data);
    __m128i hi = _mm256_extracti128_si256(data, 1);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n]), replace_missing_epi32(_mm_cvtepi8_epi32(low), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 4]), replace_missing_epi32(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 8]), replace_missing_epi32(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 12]), replace_missing_epi32(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 16]), replace_missing_epi32(_mm_cvtepi8_epi32(hi), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 20]), replace_missing_epi32(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 24]), replace_missing_epi32(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)), -128));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 28]), replace_missing_epi32(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)), -128));
  }
  return n;
}

__attribute__((target("avx2")))
static std::uint64_t parse_16bit_ints_avx2(const char * src, std::int32_t * dst, std::uint64_t max_n) {
  const __m256i missing = _mm256_set1_epi16(static_cast<short>(0x8000));
  const __m256i not_recorded = _mm256_set1_epi16(static_cast<short>(0x8001));

  std::uint64_t n = 0;
  for (; n + 16 <= max_n; n += 16) {
    __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n * 2));
    data = _mm256_blendv_epi8(data, missing, _mm256_cmpeq_epi16(data, not_recorded));

    __m128i low = _mm256_castsi256_si128(data);
    __m128i hi = _mm256_extracti128_si256(data, 1);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n]), replace_missing_epi32(_mm_cvtepi16_epi32(low), -32768));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 4]), replace_missing_epi32(_mm_cvtepi16_epi32(_mm_bsrli_si128(low, 8)), -32768));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 8]), replace_missing_epi32(_mm_cvtepi16_epi32(hi), -32768));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 12]), replace_missing_epi32(_mm_cvtepi16_epi32(_mm_bsrli_si128(hi, 8)), -32768));
  }
  return n;
}

__attribute__((target("avx2")))
static std::uint64_t parse_32bit_ints_avx2(const char * src, std::int32_t * dst, std::uint64_t max_n) {
  const __m256i missing = _mm256_set1_epi32(INT32_MIN);
  const __m256i not_recorded = _mm256_set1_epi32(INT32_MIN + 1);

  std::uint64_t n = 0;
  for (; n + 8 <= max_n; n += 8) {
    __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n * 4));
    data = _mm256_blendv_epi8(data, missing, _mm256_cmpeq_epi32(data, not_recorded));
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(&dst[n]), data);
  }
  return n;
}

__attribute__((target("avx2")))
static std::uint64_t parse_geno_avx2(const char * src, std::int32_t * dst, std::uint8_t * phase, std::uint64_t max_n) {
  const __m256i mask_geno = _mm256_set1_epi8(static_cast<char>(0xfe));
  const __m256i sub = _mm256_set1_epi8(1);
  const __m256i miss_byte = _mm256_set1_epi8(static_cast<char>(0x80));
  const __m256i eov_byte = _mm256_set1_epi8(static_cast<char>(0x81));
  const __m128i one = _mm_set1_epi8(1);
  const __m128i odd_lo = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);
  const __m128i odd_hi = _mm_setr_epi8(-1, -1, -1, -1, -1, -1, -1, -1, 1, 3, 5, 7, 9, 11, 13, 15);

  std::uint64_t n = 0;
  for (; n + 32 <= max_n; n += 32) {
    __m256i initial = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + n));

    __m256i geno = _mm256_and_si256(initial, mask_geno);
    geno = _mm256_sub_epi8(_mm256_srli_epi32(geno, 1), sub);

    // set BCF missing (0x80) and end-of-vector (0x81) entries to -1 (0xff)
    __m256i missed = _mm256_or_si256(_mm256_cmpeq_epi8(initial, miss_byte),
                                     _mm256_cmpeq_epi8(initial, eov_byte));
    geno = _mm256_or_si256(geno, missed);

    __m128i low = _mm256_castsi256_si128(geno);
    __m128i hi = _mm256_extracti128_si256(geno, 1);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n]), _mm_cvtepi8_epi32(low));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 4]), _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 8]), _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 12]), _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 16]), _mm_cvtepi8_epi32(hi));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 20]), _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 24]), _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&dst[n + 28]), _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)));

    // sample phasing from each sample's second allele (odd bytes of initial)
    __m128i init_lo = _mm256_castsi256_si128(initial);
    __m128i init_hi = _mm256_extracti128_si256(initial, 1);
    __m128i second_alleles = _mm_or_si128(_mm_shuffle_epi8(init_lo, odd_lo),
                                          _mm_shuffle_epi8(init_hi, odd_hi));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(&phase[n >> 1]), _mm_and_si128(second_alleles, one));
  }
  return n;
}
#elif defined(__aarch64__)
/// @brief replace sign-extended narrow missing sentinel with INT32_MIN
static inline int32x4_t replace_missing_s32(int32x4_t data, std::int32_t old_missing) {
  const int32x4_t missing_old = vdupq_n_s32(old_missing);
  const int32x4_t missing_32bit = vdupq_n_s32(INT32_MIN);
  return vbslq_s32(vceqq_s32(data, missing_old), missing_32bit, data);
}
#endif

/// @brief parse 8-bit ints from the buffer
std::vector<std::int32_t> SampleData::parse_8bit_ints(FormatType & type) {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  std::vector<std::int32_t> vals(max_n);
  if (max_n == 0) {
    return vals;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n > buf.size()) {
    throw std::invalid_argument("FORMAT int8 data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (__builtin_cpu_supports("avx2")) {
    n = parse_8bit_ints_avx2(buf.data() + offset, vals.data(), max_n);
  }
#elif defined(__aarch64__)
  const int8x16_t missing = vdupq_n_s8(INT8_MIN);
  const int8x16_t not_recorded = vdupq_n_s8(INT8_MIN + 1);

  for (; n + 16 <= max_n; n += 16) {
    int8x16_t data = vld1q_s8(reinterpret_cast<const std::int8_t *>(buf.data() + offset + n));
    data = vbslq_s8(vceqq_s8(data, not_recorded), missing, data);

    int16x8_t wider = vmovl_s8(vget_low_s8(data));
    vst1q_s32(&vals[n], replace_missing_s32(vmovl_s16(vget_low_s16(wider)), -128));
    vst1q_s32(&vals[n + 4], replace_missing_s32(vmovl_s16(vget_high_s16(wider)), -128));

    wider = vmovl_s8(vget_high_s8(data));
    vst1q_s32(&vals[n + 8], replace_missing_s32(vmovl_s16(vget_low_s16(wider)), -128));
    vst1q_s32(&vals[n + 12], replace_missing_s32(vmovl_s16(vget_high_s16(wider)), -128));
  }
#endif
  offset += static_cast<std::uint32_t>(n);
  for (; n < max_n; n++) {
    vals[n] = parse_int(buf.data(), buf.size(), offset, type.type_size);
  }
  return vals;
}

/// @brief parse 16-bit ints from the buffer
std::vector<std::int32_t> SampleData::parse_16bit_ints(FormatType & type) {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  std::vector<std::int32_t> vals(max_n);
  if (max_n == 0) {
    return vals;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * 2 > buf.size()) {
    throw std::invalid_argument("FORMAT int16 data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (__builtin_cpu_supports("avx2")) {
    n = parse_16bit_ints_avx2(buf.data() + offset, vals.data(), max_n);
  }
#elif defined(__aarch64__)
  const int16x8_t missing = vdupq_n_s16(INT16_MIN);
  const int16x8_t not_recorded = vdupq_n_s16(INT16_MIN + 1);

  for (; n + 8 <= max_n; n += 8) {
    int16x8_t data = vld1q_s16(reinterpret_cast<const std::int16_t *>(buf.data() + offset + n * 2));
    data = vbslq_s16(vceqq_s16(data, not_recorded), missing, data);

    vst1q_s32(&vals[n], replace_missing_s32(vmovl_s16(vget_low_s16(data)), -32768));
    vst1q_s32(&vals[n + 4], replace_missing_s32(vmovl_s16(vget_high_s16(data)), -32768));
  }
#endif
  offset += static_cast<std::uint32_t>(n * 2);
  for (; n < max_n; n++) {
    vals[n] = parse_int(buf.data(), buf.size(), offset, type.type_size);
  }
  return vals;
}

/// @brief parse 32-bit ints from the buffer
std::vector<std::int32_t> SampleData::parse_32bit_ints(FormatType & type) {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  std::vector<std::int32_t> vals(max_n);
  if (max_n == 0) {
    return vals;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * 4 > buf.size()) {
    throw std::invalid_argument("FORMAT int32 data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (__builtin_cpu_supports("avx2")) {
    n = parse_32bit_ints_avx2(buf.data() + offset, vals.data(), max_n);
  }
#elif defined(__aarch64__)
  const int32x4_t missing = vdupq_n_s32(INT32_MIN);
  const int32x4_t not_recorded = vdupq_n_s32(INT32_MIN + 1);

  for (; n + 4 <= max_n; n += 4) {
    int32x4_t data = vld1q_s32(reinterpret_cast<const std::int32_t *>(buf.data() + offset + n * 4));
    data = vbslq_s32(vceqq_s32(data, not_recorded), missing, data);
    vst1q_s32(&vals[n], data);
  }
#endif
  offset += static_cast<std::uint32_t>(n * 4);
  for (; n < max_n; n++) {
    vals[n] = parse_int(buf.data(), buf.size(), offset, type.type_size);
  }
  return vals;
}

/// @brief optimized parsing of genotype data
///
/// @param type struct containing number of values per person, and where in the
///             buffer to start parsing
/// @return vector of genotype values as 32-bit ints (missing=-1)
std::vector<std::int32_t> SampleData::get_geno(FormatType & type) {
  // confirm we checked sample phasing if we look at the genotype data
  phase_checked = true;
  phase.assign(n_samples, 0);

  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  std::vector<std::int32_t> vals(max_n);
  if (max_n == 0 || type.n_vals == 0) {
    return vals;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * type.type_size > buf.size()) {
    throw std::invalid_argument("FORMAT GT data exceeds buffer");
  }

  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  if (__builtin_cpu_supports("avx2") && (type.n_vals == 2) && (type.type_size == 1)) {
    n = parse_geno_avx2(buf.data() + offset, vals.data(), phase.data(), max_n);
  }
#elif defined(__aarch64__)
  if ((type.type_size == 1) && (type.n_vals == 2)) {
    const uint8x16_t sub = vdupq_n_u8(1);
    const uint8x16_t miss_thresh = vdupq_n_u8(0x80);
    const uint8x8_t one = vdup_n_u8(1);

    for (; n + 16 <= max_n; n += 16) {
      uint8x16_t initial = vld1q_u8(reinterpret_cast<const std::uint8_t *>(buf.data() + offset + n));
      int8x16_t geno = vreinterpretq_s8_u8(vsubq_u8(vshrq_n_u8(initial, 1), sub));

      // set BCF missing (0x80) and end-of-vector (0x81) entries to -1 (0xff)
      uint8x16_t missed = vcgeq_u8(initial, miss_thresh);
      geno = vorrq_s8(geno, vreinterpretq_s8_u8(missed));

      int16x8_t wider = vmovl_s8(vget_low_s8(geno));
      vst1q_s32(&vals[n], vmovl_s16(vget_low_s16(wider)));
      vst1q_s32(&vals[n + 4], vmovl_s16(vget_high_s16(wider)));

      wider = vmovl_s8(vget_high_s8(geno));
      vst1q_s32(&vals[n + 8], vmovl_s16(vget_low_s16(wider)));
      vst1q_s32(&vals[n + 12], vmovl_s16(vget_high_s16(wider)));

      // sample phasing from odd bytes
      uint8x8_t second_alleles = vshrn_n_u16(vreinterpretq_u16_u8(initial), 8);
      vst1_u8(&phase[n >> 1], vand_u8(second_alleles, one));
    }
  }
#endif

  offset += static_cast<std::uint32_t>(n * type.type_size);
  std::uint64_t idx = n;
  for (std::uint64_t s = n / type.n_vals; s < n_samples; s++) {
    std::uint8_t is_phased = 1;
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      std::int32_t raw = parse_int(buf.data(), buf.size(), offset, type.type_size);
      if (raw < 0) {
        vals[idx] = -1;
      } else {
        if (i > 0) {
          is_phased &= static_cast<std::uint8_t>(raw & 1);
        }
        vals[idx] = (raw >> 1) - 1;
      }
      idx++;
    }
    phase[s] = is_phased;
  }
  return vals;
}

/// @brief parse float values from the buffer
/// @param type struct with number of values per person, and buffer offset
/// @return vector of floats
std::vector<float> SampleData::get_floats(FormatType & type) {
  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  std::vector<float> vals(max_n);
  if (max_n == 0) {
    return vals;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * sizeof(float) > buf.size()) {
    throw std::invalid_argument("FORMAT float data exceeds buffer");
  }
  const char * src = buf.data() + type.offset;
  std::memcpy(vals.data(), src, max_n * sizeof(float));
  const float qnan = std::nanf("");
  for (std::uint64_t n = 0; n < max_n; n++) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, src + n * sizeof(float), sizeof(bits));
    if (is_bcf_float_missing(bits)) {
      vals[n] = qnan;
    }
  }
  return vals;
}

/// @brief parse strings from the buffer
/// @param type struct with number of bytes per person, and buffer offset
/// @return vector of strings
std::vector<std::string> SampleData::get_strings(FormatType & type) {
  std::vector<std::string> vals(n_samples);
  if (n_samples == 0 || type.n_vals == 0) {
    return vals;
  }
  std::uint32_t offset = type.offset;
  for (std::uint32_t n = 0; n < n_samples; n++) {
    vals[n] = parse_string(buf.data(), buf.size(), offset, type.n_vals);
  }
  return vals;
}

}
