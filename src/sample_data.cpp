#include <climits>
#include <cstring>
#include <stdexcept>

#include "sample_data.h"
#include "types.h"

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
#endif

namespace bcf {

SampleData::SampleData(GzReader & infile, Header & _header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples) {
  n_samples = _n_samples;
  header = &_header;
  if (len == 0) {
    return;
  }
  phase.assign(n_samples, 0);
  missing.assign(n_samples, 0);

  // read the sample data into a buffer, but don't parse until required
  buf.resize(len);
  infile.read_exact(buf.data(), len, "truncated BCF sample data");

  // read the available keys
  std::uint32_t buf_idx = 0;
  Typed type_val;
  for (std::uint32_t i = 0; i < n_fmt; i++) {
    type_val = Typed(buf.data(), buf.size(), buf_idx);
    std::int32_t format_idx = parse_int(buf.data(), buf.size(), buf_idx, type_val.type_size);
    if (format_idx < 0 || header->format.count(static_cast<std::uint32_t>(format_idx)) == 0) {
      throw std::invalid_argument("unknown FORMAT index: " + std::to_string(format_idx));
    }
    const std::string & key = header->format.at(static_cast<std::uint32_t>(format_idx)).id;
    bool is_geno = key == "GT";

    type_val = Typed(buf.data(), buf.size(), buf_idx);
    std::uint64_t field_bytes = static_cast<std::uint64_t>(type_val.n_vals) *
                                static_cast<std::uint64_t>(type_val.type_size) *
                                static_cast<std::uint64_t>(n_samples);
    if (static_cast<std::uint64_t>(buf_idx) + field_bytes > buf.size()) {
      throw std::out_of_range("FORMAT field exceeds BCF sample buffer");
    }
    keys[key] = {static_cast<std::uint8_t>(type_val.type), type_val.type_size, buf_idx,
                 type_val.n_vals, is_geno};
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

/// @brief figure out the data type for a given format key
/// @param key string for the format key to check
/// @return struct containing details of the type, number of values, and buffer offset
FormatType SampleData::get_type(std::string &key) {
  if (keys.count(key) == 0) {
    throw std::invalid_argument("no entries for " + key + " in data");
  }
  return keys[key];
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

#if defined(__x86_64__)
/// @brief convert vectorized 8-bit missing values to the 32-bit form
/// @param data 128-bit register containing 4 32-bit values
/// @return register where the missing values have been converted from their
///         8-bit encoding to a standard 32-bit encoding
static inline __m128i missing_8bit_to_32bit(__m128i data) {
  // 8-bit missing values were 0x80 (-128), sign-extended to 0xffffff80
  const __m128i missing_8bit = _mm_set1_epi32(-128);
  const __m128i missing_32bit = _mm_set1_epi32(INT32_MIN);

  __m128i mask = _mm_cmpeq_epi32(data, missing_8bit);
  data = _mm_andnot_si128(mask, data);
  return _mm_or_si128(data, _mm_and_si128(missing_32bit, mask));
}
#elif defined(__aarch64__)
/// @brief convert vectorized 8-bit missing values to the 32-bit form
static inline int32x4_t missing_8bit_to_32bit(int32x4_t data) {
  const int32x4_t missing_8bit = vdupq_n_s32(-128);
  const int32x4_t missing_32bit = vdupq_n_s32(INT32_MIN);

  uint32x4_t mask = vceqq_s32(data, missing_8bit);
  return vbslq_s32(mask, missing_32bit, data);
}
#endif

#if defined(__x86_64__)
/// @brief convert vectorized 16-bit missing values to the 32-bit form
static inline __m128i missing_16bit_to_32bit(__m128i data) {
  // 16-bit missing values were 0x8000 (-32768), sign-extended to 0xffff8000
  const __m128i missing_16bit = _mm_set1_epi32(-32768);
  const __m128i missing_32bit = _mm_set1_epi32(INT32_MIN);

  __m128i mask = _mm_cmpeq_epi32(data, missing_16bit);
  data = _mm_andnot_si128(mask, data);
  return _mm_or_si128(data, _mm_and_si128(missing_32bit, mask));
}
#elif defined(__aarch64__)
/// @brief convert vectorized 16-bit missing values to the 32-bit form
static inline int32x4_t missing_16bit_to_32bit(int32x4_t data) {
  const int32x4_t missing_16bit = vdupq_n_s32(-32768);
  const int32x4_t missing_32bit = vdupq_n_s32(INT32_MIN);

  uint32x4_t mask = vceqq_s32(data, missing_16bit);
  return vbslq_s32(mask, missing_32bit, data);
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
    throw std::out_of_range("FORMAT int8 data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    const __m256i missing = _mm256_set1_epi8(static_cast<char>(0x80));
    const __m256i not_recorded = _mm256_set1_epi8(static_cast<char>(0x81));

    for (; n + 32 <= max_n; n += 32) {
      __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(buf.data() + offset + n));

      // replace the not_recorded values with standard missing values
      __m256i mask = _mm256_cmpeq_epi8(data, not_recorded);
      data = _mm256_or_si256(_mm256_andnot_si256(mask, data), _mm256_and_si256(missing, mask));

      // widen 8-bit ints to 32-bit, standardize missing value, and store
      __m128i low = _mm256_castsi256_si128(data);
      __m128i hi = _mm256_extracti128_si256(data, 1);
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(low)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 4]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4))));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 8]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8))));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 12]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12))));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 16]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(hi)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 20]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4))));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 24]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8))));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 28]), missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12))));
    }
  }
#elif defined(__aarch64__)
  const int8x16_t missing = vdupq_n_s8(INT8_MIN);
  const int8x16_t not_recorded = vdupq_n_s8(INT8_MIN + 1);

  for (; n + 16 <= max_n; n += 16) {
    int8x16_t data = vld1q_s8(reinterpret_cast<const std::int8_t *>(buf.data() + offset + n));
    uint8x16_t mask = vceqq_s8(data, not_recorded);
    data = vbslq_s8(mask, missing, data);

    int16x8_t wider = vmovl_s8(vget_low_s8(data));
    vst1q_s32(&vals[n], missing_8bit_to_32bit(vmovl_s16(vget_low_s16(wider))));
    vst1q_s32(&vals[n + 4], missing_8bit_to_32bit(vmovl_s16(vget_high_s16(wider))));

    wider = vmovl_s8(vget_high_s8(data));
    vst1q_s32(&vals[n + 8], missing_8bit_to_32bit(vmovl_s16(vget_low_s16(wider))));
    vst1q_s32(&vals[n + 12], missing_8bit_to_32bit(vmovl_s16(vget_high_s16(wider))));
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
    throw std::out_of_range("FORMAT int16 data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    const __m256i missing = _mm256_set1_epi16(static_cast<short>(0x8000));
    const __m256i not_recorded = _mm256_set1_epi16(static_cast<short>(0x8001));

    for (; n + 16 <= max_n; n += 16) {
      __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(buf.data() + offset + n * 2));

      __m256i mask = _mm256_cmpeq_epi16(data, not_recorded);
      data = _mm256_or_si256(_mm256_andnot_si256(mask, data), _mm256_and_si256(missing, mask));

      __m128i low = _mm256_castsi256_si128(data);
      __m128i hi = _mm256_extracti128_si256(data, 1);
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n]), missing_16bit_to_32bit(_mm_cvtepi16_epi32(low)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 4]), missing_16bit_to_32bit(_mm_cvtepi16_epi32(_mm_bsrli_si128(low, 8))));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 8]), missing_16bit_to_32bit(_mm_cvtepi16_epi32(hi)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 12]), missing_16bit_to_32bit(_mm_cvtepi16_epi32(_mm_bsrli_si128(hi, 8))));
    }
  }
#elif defined(__aarch64__)
  const int16x8_t missing = vdupq_n_s16(INT16_MIN);
  const int16x8_t not_recorded = vdupq_n_s16(INT16_MIN + 1);

  for (; n + 8 <= max_n; n += 8) {
    int16x8_t data = vld1q_s16(reinterpret_cast<const std::int16_t *>(buf.data() + offset + n * 2));
    uint16x8_t mask = vceqq_s16(data, not_recorded);
    data = vbslq_s16(mask, missing, data);

    vst1q_s32(&vals[n], missing_16bit_to_32bit(vmovl_s16(vget_low_s16(data))));
    vst1q_s32(&vals[n + 4], missing_16bit_to_32bit(vmovl_s16(vget_high_s16(data))));
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
    throw std::out_of_range("FORMAT int32 data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    const __m256i missing = _mm256_set1_epi32(INT32_MIN);
    const __m256i not_recorded = _mm256_set1_epi32(INT32_MIN + 1);
    for (; n + 8 <= max_n; n += 8) {
      __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(buf.data() + offset + n * 4));

      __m256i mask = _mm256_cmpeq_epi32(data, not_recorded);
      data = _mm256_or_si256(_mm256_andnot_si256(mask, data), _mm256_and_si256(missing, mask));

      _mm256_storeu_si256(reinterpret_cast<__m256i *>(&vals[n]), data);
    }
  }
#elif defined(__aarch64__)
  const int32x4_t missing = vdupq_n_s32(INT32_MIN);
  const int32x4_t not_recorded = vdupq_n_s32(INT32_MIN + 1);

  for (; n + 4 <= max_n; n += 4) {
    int32x4_t data = vld1q_s32(reinterpret_cast<const std::int32_t *>(buf.data() + offset + n * 4));
    uint32x4_t mask = vceqq_s32(data, not_recorded);
    data = vbslq_s32(mask, missing, data);

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

  std::uint64_t max_n = static_cast<std::uint64_t>(type.n_vals) * n_samples;
  std::vector<std::int32_t> vals(max_n);
  if (max_n == 0 || type.n_vals == 0) {
    return vals;
  }
  if (static_cast<std::uint64_t>(type.offset) + max_n * type.type_size > buf.size()) {
    throw std::out_of_range("FORMAT GT data exceeds buffer");
  }

  std::uint32_t offset = type.offset;
  std::uint64_t n = 0;
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2") && (type.n_vals == 2) && (type.type_size == 1)) {
    const __m256i mask_geno = _mm256_set1_epi8(static_cast<char>(0xfe));
    const __m256i sub = _mm256_set1_epi8(1);
    const __m256i miss_byte = _mm256_set1_epi8(static_cast<char>(0x80));
    const __m256i eov_byte = _mm256_set1_epi8(static_cast<char>(0x81));
    const __m128i neg_one = _mm_set1_epi8(-1);
    const __m128i one = _mm_set1_epi8(1);
    const __m128i even_lo = _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i even_hi = _mm_setr_epi8(-1, -1, -1, -1, -1, -1, -1, -1, 0, 2, 4, 6, 8, 10, 12, 14);
    const __m128i odd_lo = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i odd_hi = _mm_setr_epi8(-1, -1, -1, -1, -1, -1, -1, -1, 1, 3, 5, 7, 9, 11, 13, 15);

    for (; n + 32 <= max_n; n += 32) {
      __m256i initial = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(buf.data() + offset + n));

      __m256i geno = _mm256_and_si256(initial, mask_geno);
      geno = _mm256_sub_epi8(_mm256_srli_epi32(geno, 1), sub);

      // set BCF missing (0x80) and end-of-vector (0x81) entries to -1 (0xff)
      __m256i missed = _mm256_or_si256(_mm256_cmpeq_epi8(initial, miss_byte),
                                       _mm256_cmpeq_epi8(initial, eov_byte));
      geno = _mm256_or_si256(geno, missed);

      __m128i low = _mm256_castsi256_si128(geno);
      __m128i hi = _mm256_extracti128_si256(geno, 1);
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n]), _mm_cvtepi8_epi32(low));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 4]), _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 8]), _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 12]), _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 16]), _mm_cvtepi8_epi32(hi));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 20]), _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 24]), _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)));
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&vals[n + 28]), _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)));

      // sample missingness from each sample's first allele (even bytes)
      __m128i first_alleles = _mm_or_si128(_mm_shuffle_epi8(low, even_lo),
                                           _mm_shuffle_epi8(hi, even_hi));
      __m128i is_miss = _mm_and_si128(_mm_cmpeq_epi8(first_alleles, neg_one), one);
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&missing[n >> 1]), is_miss);

      // sample phasing from each sample's second allele (odd bytes of initial)
      __m128i init_lo = _mm256_castsi256_si128(initial);
      __m128i init_hi = _mm256_extracti128_si256(initial, 1);
      __m128i second_alleles = _mm_or_si128(_mm_shuffle_epi8(init_lo, odd_lo),
                                            _mm_shuffle_epi8(init_hi, odd_hi));
      __m128i phase128 = _mm_and_si128(second_alleles, one);
      _mm_storeu_si128(reinterpret_cast<__m128i *>(&phase[n >> 1]), phase128);
    }
  }
#elif defined(__aarch64__)
  if ((type.type_size == 1) && (type.n_vals == 2)) {
    const uint8x16_t sub = vdupq_n_u8(1);
    const uint8x16_t miss_thresh = vdupq_n_u8(0x80);
    const uint8x8_t all_ones = vdup_n_u8(0xff);
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

      // sample missingness from even bytes
      uint8x8_t first_alleles = vmovn_u16(vreinterpretq_u16_s8(geno));
      uint8x8_t is_miss = vand_u8(vceq_u8(first_alleles, all_ones), one);
      vst1_u8(&missing[n >> 1], is_miss);

      // sample phasing from odd bytes
      uint8x8_t second_alleles = vshrn_n_u16(vreinterpretq_u16_u8(initial), 8);
      uint8x8_t is_phased = vand_u8(second_alleles, one);
      vst1_u8(&phase[n >> 1], is_phased);
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
      if (i == 0) {
        missing[s] = static_cast<std::uint8_t>(vals[idx] == -1);
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
  if (static_cast<std::uint64_t>(type.offset) + max_n * 4 > buf.size()) {
    throw std::out_of_range("FORMAT float data exceeds buffer");
  }
  std::uint32_t offset = type.offset;
  for (std::uint64_t n = 0; n < max_n; n++) {
    vals[n] = parse_float(buf.data(), buf.size(), offset);
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
