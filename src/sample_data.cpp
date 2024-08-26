
#include <cstring>

#include "sample_data.h"
#include "types.h"

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
#endif

namespace bcf {

SampleData::SampleData(igzstream & infile, Header & _header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples) {
  n_samples = _n_samples;
  header = &_header;
  if (len == 0) {
    return;
  }
  phase.resize(n_samples);
  missing.resize(n_samples);
  
  // read the sample data into a buffer, but don't parse until required
  buf.resize(len);
  infile.read(reinterpret_cast<char *>(&buf[0]), len);
  
  // read the available keys
  std::uint32_t buf_idx=0;
  std::uint32_t format_idx=0;
  std::string key;
  Typed type_val;
  bool is_geno;
  for (std::uint32_t i = 0; i < n_fmt; i++ ){
    type_val = {&buf[0], buf_idx};
    format_idx = parse_int(&buf[0], buf_idx, type_val.type_size);
    key = header->format[format_idx].id;
    is_geno = key == "GT";

    type_val = {&buf[0], buf_idx};
    keys[key] = {(std::uint8_t) type_val.type, type_val.type_size, buf_idx, 
                 type_val.n_vals, is_geno};
    buf_idx += (type_val.n_vals * type_val.type_size * n_samples);
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
/// @return 
std::vector<std::int32_t> SampleData::get_ints(FormatType & type) {
  if (type.is_geno) {
    return get_geno(type);
  } else if (type.type_size == 1) {
    return parse_8bit_ints(type);
  } else if (type.type_size == 2) {
    parse_16bit_ints(type);
  } else if (type.type_size == 4) {
    parse_32bit_ints(type);
  }
  
  std::vector<std::int32_t> vals;
  vals.resize(type.n_vals * n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t n=0;
  for (; n < (n_samples * type.n_vals); n++) {
    vals[n] = parse_int(&buf[0], offset, type.type_size);
  }
  return vals;
}

#if defined(__x86_64__)
/// @brief convert vectorized 8-bit missing values to the 32-bit form
/// @param data 128-bit register containing 4 32-bit values
/// @return register where the missing values have been converted from their
///         8-bit encoding to a standard 32-bit encoding
__m128i missing_8bit_to_32bit(__m128i data) {
  __m128 mask;
  // 8-bit missing values were 0x80, but after extended to 32-bit became 0xffffff80
  // we need to convert these to 0x80000000, which is the missing value for 32-bit values
  const __m128i missing_8bit = _mm_set_epi32(0xffffff80, 0xffffff80, 0xffffff80, 0xffffff80);
  const __m128i missing_32bit = _mm_set_epi32(0x80000000, 0x80000000, 0x80000000, 0x80000000);
  
  // find which entries have missing values
  mask = (__m128) _mm_cmpeq_epi32(data, missing_8bit);
  
  // erase original missing values
  data = (__m128i) _mm_andnot_ps((__m128) mask, (__m128) data);
  
  // replace 8-bit missing value with 32-bit missing value
  return (__m128i) _mm_or_ps((__m128) data, _mm_and_ps((__m128) missing_32bit, (__m128) mask));
}
#elif defined(__aarch64__)
/// @brief convert vectorized 8-bit missing values to the 32-bit form
/// @param data 128-bit register containing 4 32-bit values
/// @return register where the missing values have been converted from their
///         8-bit encoding to a standard 32-bit encoding
int32x4_t missing_8bit_to_32bit(int32x4_t data) {
  int32x4_t mask;
  // 8-bit missing values were 0x80, but after extended to 32-bit became 0xffffff80
  // we need to convert these to 0x80000000, which is the missing value for 32-bit values
  const int32x4_t missing_8bit = vdupq_n_s32(0xffffff80);
  const int32x4_t missing_32bit = vdupq_n_s32(0x80000000);
  
  // find which entries have missing values
  mask = vceqq_s32(data, missing_8bit);
  
  // erase original missing values
  data = vandq_s32(data, vmvnq_s32(mask));
  
  // replace 8-bit missing value with 32-bit missing value
  return vorrq_s32(data, vandq_s32(missing_32bit, mask));
}
#endif

#if defined(__x86_64__)
/// @brief convert vectorized 16-bit missing values to the 32-bit form
/// @param data 128-bit register containing 4 32-bit values
/// @return register where the missing values have been converted from their
///         16-bit encoding to a standard 32-bit encoding
__m128i missing_16bit_to_32bit(__m128i data) {
  __m128 mask;
  // 16-bit missing values were 0x8000, but after extended to 32-bit became 0xffff8000
  // we need to convert these to 0x80000000, which is the missing value for 32-bit values
  const __m128i missing_16bit = _mm_set_epi32(0xffff8000, 0xffff8000, 0xffff8000, 0xffff8000);
  const __m128i missing_32bit = _mm_set_epi32(0x80000000, 0x80000000, 0x80000000, 0x80000000);
  
  // find which entries have missing values
  mask = (__m128) _mm_cmpeq_epi32(data, missing_16bit);
  
  // erase original missing values
  data = (__m128i) _mm_andnot_ps((__m128) mask, (__m128) data);
  
  // replace 16-bit missing value with 32-bit missing value
  return (__m128i) _mm_or_ps((__m128) data, _mm_and_ps((__m128) missing_32bit, (__m128) mask));
}
#elif defined(__aarch64__)
/// @brief convert vectorized 16-bit missing values to the 32-bit form
/// @param data 128-bit register containing 4 32-bit values
/// @return register where the missing values have been converted from their
///         16-bit encoding to a standard 32-bit encoding
int32x4_t missing_16bit_to_32bit(int32x4_t data) {
  int32x4_t mask;
  // 16-bit missing values were 0x8000, but after extended to 32-bit became 0xffff8000
  // we need to convert these to 0x80000000, which is the missing value for 32-bit values
  const int32x4_t missing_16bit = vdupq_n_s32(0xffff8000);
  const int32x4_t missing_32bit = vdupq_n_s32(0x80000000);
  
  // find which entries have missing values
  mask = vceqq_s32(data, missing_16bit);
  
  // erase original missing values
  data = vandq_s32(data, vmvnq_s32(mask));
  
  // replace 16-bit missing value with 32-bit missing value
  return vorrq_s32(data, vandq_s32(missing_32bit, mask));
}
#endif

/// @brief parse 8-bit ints from the buffer
///
/// This uses vectorized operations if available on x86_64 and aarch64
///
/// @param type struct containing number of values per sample (so we can determine 
///             how many values to parse), and offset within the buffer (so we
///             know where to start).
/// @return vector of data values as 32-bit ints
std::vector<std::int32_t> SampleData::parse_8bit_ints(FormatType & type) {
  std::vector<std::int32_t> vals;
  std::uint64_t max_n = type.n_vals * n_samples;
  vals.resize(max_n);
  std::uint32_t offset = type.offset;
  std::uint32_t n=0;
  
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    __m256i data, mask;
    __m128i low, hi;
    // There are two types of missing data values, one for actual missing data,
    // and one for data not recorded. We collapse these into a single missing
    // value, since we later convert all the values to float (outisde the c++ 
    // code).
    // std::int32_t missing - 8bit: 0x80, 16bit: 0x8000 etc
    // std::int32_t not_recorded - missing | 0x1 8bit: 0x81, 16bit: 0x8001, 32bit: 0x80000001
    __m256i missing = _mm256_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x80808080,
                                      0x80808080, 0x80808080, 0x80808080, 0x80808080);
    __m256i not_recorded = _mm256_set_epi32(0x81818181, 0x81818181, 0x81818181, 0x81818181,
                                           0x81818181, 0x81818181, 0x81818181, 0x81818181);
    
    for (; n < (max_n - (max_n % 32)); n += 32) {
      data = _mm256_loadu_si256((__m256i *) &buf[offset + n]);
      
      // replace the not_recorded values with standard missing values
      mask = _mm256_cmpeq_epi8(data, not_recorded);
      data = (__m256i) _mm256_andnot_ps((__m256) mask, (__m256) data);  // erase original not recorded values
      data = (__m256i) _mm256_or_ps((__m256) data, _mm256_and_ps((__m256) missing, (__m256) mask));  // replace not recorded values with standard missing value
      
      // widen 8-bit ints to 16 then 32bit, standardize missing value, and store
      low = _mm256_extractf128_si256(data, 0);
      hi = _mm256_extractf128_si256(data, 1);
      _mm_storeu_ps((float *) &vals[n], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(low)));
      _mm_storeu_ps((float *) &vals[n + 4], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4))));
      _mm_storeu_ps((float *) &vals[n + 8], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8))));
      _mm_storeu_ps((float *) &vals[n + 12], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12))));
      _mm_storeu_ps((float *) &vals[n + 16], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(hi)));
      _mm_storeu_ps((float *) &vals[n + 20], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4))));
      _mm_storeu_ps((float *) &vals[n + 24], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8))));
      _mm_storeu_ps((float *) &vals[n + 28], (__m128) missing_8bit_to_32bit(_mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12))));
    }
  }
#elif defined(__aarch64__)
  int8x16_t data, mask;
  int16x8_t wider;
  int8x16_t missing = vdupq_n_s8(0x80);
  int8x16_t not_recorded = vdupq_n_s8(0x81);
  
  for (; n < (max_n - (max_n % 16)); n += 16) {
    // load data from the array into SIMD registers.
    data = vld1q_s8((std::int8_t *)&buf[offset + n]);

    // replace the not_recorded values with standard missing values
    mask = vceqq_s8(data, not_recorded);            // find and set missing values
    data = vandq_s8(data, vmvnq_s8(mask));          // erase original missing values
    data = vorrq_s8(data, vandq_s8(missing, mask)); // swap in new missing values

    // store genotypes as 32-bit ints, have to expand all values in turn
    wider = vmovl_s8(vget_low_s8(data));
    vst1q_s32(&vals[n], missing_8bit_to_32bit(vmovl_s16(vget_low_s16(wider))));
    vst1q_s32(&vals[n + 4], missing_8bit_to_32bit(vmovl_s16(vget_high_s16(wider))));

    wider = vmovl_s8(vget_high_s8(data));
    vst1q_s32(&vals[n + 8], missing_8bit_to_32bit(vmovl_s16(vget_low_s16(wider))));
    vst1q_s32(&vals[n + 12], missing_8bit_to_32bit(vmovl_s16(vget_high_s16(wider))));
  }
#endif
  offset += n;
  for (; n < (n_samples * type.n_vals); n++) {
    vals[n] = parse_int(&buf[0], offset, type.type_size);
  }
  return vals;
}

/// @brief parse 16-bit ints from the buffer
///
/// This uses vectorized operations if available on x86_64 and aarch64
///
/// @param type struct containing number of values per sample (so we can determine
///             how many values to parse), and offset within the buffer (so we
///             know where to start).
/// @return vector of data values as 32-bit ints
std::vector<std::int32_t> SampleData::parse_16bit_ints(FormatType & type) {
  std::vector<std::int32_t> vals;
  std::uint64_t max_n = type.n_vals * n_samples;
  vals.resize(max_n);
  std::uint32_t offset = type.offset;
  std::uint32_t n=0;
  
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    __m256i data, mask;
    __m128i low, hi;
    // std::int32_t missing - 8bit: 0x80, 16bit: 0x8000 etc
    // std::int32_t not_recorded - missing | 0x1 8bit: 0x81, 16bit: 0x8001, 32bit: 0x80000001
    __m256i missing = _mm256_set_epi32(0x80008000, 0x80808080, 0x80808080, 0x80808080,
                                      0x80808080, 0x80808080, 0x80808080, 0x80808080);
    __m256i not_recorded = _mm256_set_epi32(0x80018001, 0x80018001, 0x80018001, 0x80018001,
                                           0x80018001, 0x80018001, 0x80018001, 0x80018001);
    
    for (; n < (max_n - (max_n % 16)); n += 16) {
      data = _mm256_loadu_si256((__m256i *) &buf[offset + n]);

      // replace the not_recorded values with standard missing values
      mask = _mm256_cmpeq_epi16(data, not_recorded);
      data = (__m256i) _mm256_andnot_ps((__m256) mask, (__m256) data);  // erase original not recorded values
      data = (__m256i) _mm256_or_ps((__m256) data, _mm256_and_ps((__m256) missing, (__m256) mask));  // replace not recorded values with standard missing value
      
      // widen 8-bit ints to 16 then 32bit, standardize missing value, and store
      low = _mm256_extractf128_si256(data, 0);
      hi = _mm256_extractf128_si256(data, 1);
      _mm_storeu_ps((float *) &vals[n], (__m128) missing_16bit_to_32bit(_mm_cvtepi16_epi32(low)));
      _mm_storeu_ps((float *) &vals[n + 4], (__m128) missing_16bit_to_32bit(_mm_cvtepi16_epi32(_mm_bsrli_si128(low, 8))));
      _mm_storeu_ps((float *) &vals[n + 8], (__m128) missing_16bit_to_32bit(_mm_cvtepi16_epi32(hi)));
      _mm_storeu_ps((float *) &vals[n + 12], (__m128) missing_16bit_to_32bit(_mm_cvtepi16_epi32(_mm_bsrli_si128(hi, 8))));
    }
  }
#elif defined(__aarch64__)
  int16x8_t data, mask;
  int16x8_t missing = vdupq_n_s16(0x8000);
  int16x8_t not_recorded = vdupq_n_s16(0x8001);
  
  for (; n < (max_n - (max_n % 8)); n += 8) {
    // load data from the array into SIMD registers.
    data = vld1q_s16((std::int16_t *)&buf[offset + n]);

    // replace the not_recorded values with standard missing values
    mask = vceqq_s16(data, not_recorded);            // find and set missing values
    data = vandq_s16(data, vmvnq_s16(mask));          // erase original missing values
    data = vorrq_s16(data, vandq_s16(missing, mask)); // swap in new missing values

    // store genotypes as 32-bit ints, have to expand all values in turn
    vst1q_s32(&vals[n], missing_16bit_to_32bit(vmovl_s16(vget_low_s16(data))));
    vst1q_s32(&vals[n + 4], missing_16bit_to_32bit(vmovl_s16(vget_high_s16(data))));
  }
#endif
  offset += n;
  for (; n < (n_samples * type.n_vals); n++) {
    vals[n] = parse_int(&buf[0], offset, type.type_size);
  }
  return vals;
}

/// @brief parse 32-bit ints from the buffer
///
/// This uses vectorized operations if available on x86_64 and aarch64
///
/// @param type struct containing number of values per sample (so we can determine
///             how many values to parse), and offset within the buffer (so we
///             know where to start).
/// @return vector of data values as 32-bit ints
std::vector<std::int32_t> SampleData::parse_32bit_ints(FormatType & type) {
  std::vector<std::int32_t> vals;
  std::uint64_t max_n = type.n_vals * n_samples;
  vals.resize(max_n);
  std::uint32_t offset = type.offset;
  std::uint32_t n=0;
  
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    __m256i data, mask;
    // std::int32_t missing - 8bit: 0x80, 16bit: 0x8000 etc
    // std::int32_t not_recorded - missing | 0x1 8bit: 0x81, 16bit: 0x8001, 32bit: 0x80000001
    __m256i missing = _mm256_set_epi32(0x80000000, 0x80000000, 0x80000000, 0x80000000,
                                      0x80000000, 0x80000000, 0x80000000, 0x80000000);
    __m256i not_recorded = _mm256_set_epi32(0x80000001, 0x80000001, 0x80000001, 0x80000001,
                                           0x80000001, 0x80000001, 0x80000001, 0x80000001);
    for (; n < (max_n - (max_n % 8)); n += 8) {
      data = _mm256_loadu_si256((__m256i *) &buf[offset + n]);

      // replace the not_recorded values with standard missing values
      mask = _mm256_cmpeq_epi16(data, not_recorded);
      data = (__m256i) _mm256_andnot_ps((__m256) mask, (__m256) data);  // erase original not recorded values
      data = (__m256i) _mm256_or_ps((__m256) data, _mm256_and_ps((__m256) missing, (__m256) mask));  // replace not recorded values with standard missing value
      
      // store data
      _mm256_storeu_ps((float *) &vals[n], (__m256) data);
    }
  }
#elif defined(__aarch64__)
  int32x4_t data, mask;
  int32x4_t missing = vdupq_n_s32(0x8000000);
  int32x4_t not_recorded = vdupq_n_s32(0x81);
  
  for (; n < (max_n - (max_n % 4)); n += 4) {
    // load data from the array into SIMD registers.
    data = vld1q_s32((std::int32_t *)&buf[offset + n]);

    // replace the not_recorded values with standard missing values
    mask = vceqq_s32(data, not_recorded);            // find and set missing values
    data = vandq_s32(data, vmvnq_s32(mask));          // erase original missing values
    data = vorrq_s32(data, vandq_s32(missing, mask)); // swap in new missing values
    
    vst1q_s32(&vals[n], data);
  }
#endif
  offset += n;
  for (; n < (n_samples * type.n_vals); n++) {
    vals[n] = parse_int(&buf[0], offset, type.type_size);
  }
  return vals;
}

/// @brief optimized parsing of genotype data
///
/// The genotypes are stored as ints, so we could use the other int parsing
/// functions to extract the values, but the genotypes also store phase info,
/// are offset by 1, and have a genotype-specific missing value, so it's easier
/// to have a function dedictated to parsing genotypes only. This is optimized
/// for parsing 8-bit genotypes on x84_64 and aarch64.
///
/// @param type struct containing number of values per person, and where in the
///             buffer to start parsing
/// @return vector of genotype values as 32-bit ints (missing=0x80000000)
std::vector<std::int32_t> SampleData::get_geno(FormatType & type) {
  // confirm we checked sample phasing if we look at the genotype data
  phase_checked = true;
  
  std::vector<std::int32_t> vals;
  std::uint64_t max_n = type.n_vals * n_samples;
  vals.resize(max_n);
  std::uint32_t offset = type.offset;
  std::uint32_t n=0;
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2") && (type.n_vals == 2) && (type.type_size == 1)) {
    __m256i initial, geno, phase_vec, missed;
    __m128i low, hi, phase128;
    __m256i mask_phase = _mm256_set_epi32(0x01000100, 0x01000100, 0x01000100, 0x01000100,
                                       0x01000100, 0x01000100, 0x01000100, 0x01000100);
    __m256i mask_geno = _mm256_set_epi32(0xfefefefe, 0xfefefefe, 0xfefefefe, 0xfefefefe,
                                      0xfefefefe, 0xfefefefe, 0xfefefefe, 0xfefefefe);
    __m256i sub = _mm256_set_epi64x(0x0101010101010101, 0x0101010101010101, 0x0101010101010101, 0x0101010101010101);
    __m128i missing_mask = _mm_set_epi32(0x00ff00ff, 0x00ff00ff, 0x00ff00ff, 0x00ff00ff);
    __m128i missing_geno = _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    __m256i not_recorded = _mm256_set_epi32(0x81818181, 0x81818181, 0x81818181, 0x81818181,
                                                 0x81818181, 0x81818181, 0x81818181, 0x81818181);
    __m128i shuffle = _mm_set_epi8(0, 2, 4, 6, 8, 10, 12, 14, 17, 3, 5, 7, 9, 11, 13, 15);
    
    for (; n < (max_n - (max_n % 32)); n += 32) {
      initial = _mm256_loadu_si256((__m256i *) &buf[offset + n]);

      geno = _mm256_and_si256(initial, mask_geno);
      geno = _mm256_sub_epi8(_mm256_srli_epi32(geno, 1), sub);
      
      // account for missing values (due to different ploidy between samples)
      missed = _mm256_cmpeq_epi8(initial, not_recorded);           // find missing values
      geno = (__m256i) _mm256_andnot_ps((__m256)missed, (__m256)geno);  // erase original missing values
      geno = (__m256i) _mm256_or_ps((__m256)geno, (__m256)missed);      // swap in new missing values

      // expand the first 8 values to 32-bits, and store
      low = _mm256_extractf128_si256(geno, 0);
      hi = _mm256_extractf128_si256(geno, 1);
      _mm_storeu_ps((float *) &vals[n], (__m128) _mm_cvtepi8_epi32(low));
      _mm_storeu_ps((float *) &vals[n + 4], (__m128) _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 4)));
      _mm_storeu_ps((float *) &vals[n + 8], (__m128) _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 8)));
      _mm_storeu_ps((float *) &vals[n + 12], (__m128) _mm_cvtepi8_epi32(_mm_bsrli_si128(low, 12)));
      _mm_storeu_ps((float *) &vals[n + 16], (__m128) _mm_cvtepi8_epi32(hi));
      _mm_storeu_ps((float *) &vals[n + 20], (__m128) _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 4)));
      _mm_storeu_ps((float *) &vals[n + 24], (__m128) _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 8)));
      _mm_storeu_ps((float *) &vals[n + 28], (__m128) _mm_cvtepi8_epi32(_mm_bsrli_si128(hi, 12)));
      
      // check for missing genotypes. We can check if every second genotype is
      // -1, since that is the missing genotype indicator, then shuffle the data
      // to remove the interspersed bytes.
      low = _mm_or_si128(low, missing_mask);
      hi = _mm_or_si128(hi, missing_mask);
      low = _mm_or_si128(low, _mm_bsrli_si128(hi, 1));
      low = _mm_abs_epi8(_mm_and_si128(low, missing_geno));
      low = _mm_shuffle_epi8(low, shuffle);
      _mm_storeu_ps((float *) &missing[n >> 1], (__m128) low);

      // reorganize the phase data into correctly sorted form. Phase data is
      // initially every second byte across the m256 register. First convert to
      // two, m128 registers, interleave those, then shuffle to correct order.
      phase_vec = _mm256_and_si256(initial, mask_phase);
      low = _mm256_extractf128_si256(phase_vec, 0);
      hi = _mm256_extractf128_si256(phase_vec, 1);
      
      phase128 = _mm_or_si128(_mm_bsrli_si128(low, 1), hi);
      phase128 = _mm_shuffle_epi8(phase128, shuffle);
      _mm_storeu_ps((float *) &phase[n >> 1], (__m128)phase128);
    }
  }
#elif defined(__aarch64__)
  if ((type.type_size == 1) && (type.n_vals == 2)) {

    int8x16_t initial, geno, missed;
    uint16x8_t wider;
    int8x8_t shrunk;

    uint8x16_t missing_mask = vdupq_n_u64(0x00ff00ff00ff00ff);
    uint8x16_t mask_phase = vdupq_n_u64(0x0001000100010001);
    uint8x16_t mask_geno = vdupq_n_u8(0xfe);
    uint8x16_t sub = vdupq_n_u8(0x01);
    int8x8_t missing_geno = vdup_n_s8(-1);
    int8x16_t not_recorded = vdupq_n_s8(0x80);
    
    for (; n < (max_n - (max_n % 16)); n += 16) {
      // load data from the array into SIMD registers.
      initial = vld1q_s8((std::int8_t *)&buf[offset + n]);

      geno = vandq_s8(initial, mask_geno);
      geno = vsubq_s8(vshrq_n_s8(geno, 1), sub); // shift right to remove phase bit,
                                                 // and subtract 1 to get allele

      // account for missing values (due to different ploidy between samples)
      missed = vceqq_s8(initial, not_recorded);  // find and set missing values
      geno = vandq_s8(geno, vmvnq_s8(missed));   // erase original missing values
      geno = vorrq_s8(geno, missed);             // swap in new missing values

      // store genotypes as 32-bit ints, have to expand all values in turn
      wider = vmovl_s8(vget_low_s8(geno));
      vst1q_s32(&vals[n], vmovl_s16(vget_low_s16(wider)));
      vst1q_s32(&vals[n + 4], vmovl_s16(vget_high_s16(wider)));

      wider = vmovl_s8(vget_high_s8(geno));
      vst1q_s32(&vals[n + 8], vmovl_s16(vget_low_s16(wider)));
      vst1q_s32(&vals[n + 12], vmovl_s16(vget_high_s16(wider)));

      // check for missing genotypes
      geno = vandq_s8(geno, missing_mask);  // mask out every second value
      shrunk = vmovn_s16(geno);  // narrow to remove interspersed bytes
      shrunk = vabs_s8(vand_s8(shrunk, missing_geno));  // check if value == -1
      vst1_u8(&missing[n >> 1], shrunk);

      // check if each sample has phased data
      initial = vandq_s8(initial, mask_phase); // keep one mask bit per sample
      shrunk = vmovn_s16(initial); // narrow to remove interspersed bytes
      vst1_u8(&phase[n >> 1], shrunk);
    }
  }
#endif
  
  std::int32_t missing_indicator = 1 << ((8 * type.type_size) - 1);
  offset += n;
  std::uint32_t idx=n;
  n = n / type.n_vals;
  for (; n < n_samples; n++) {
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      vals[idx] = parse_int(&buf[0], offset, type.type_size);
      if (vals[idx] == missing_indicator) {
        vals[idx] = 0;  // convert missing values to missing genotypes
      }
      phase[n] = vals[idx] & 0x00000001;
      vals[idx] = (vals[idx] >> 1) - 1;
      // this only checks on genotype status, but this should apply to other
      // fields too (AD, DP etc), as if a sample lacks gt, other fields 
      // should also be absent
      missing[n] = vals[idx] == -1;
      idx++;
    }
  }
  return vals;
}

/// @brief parse float values from the buffer
/// @param type struct with number of samples per person, and buffer offset
/// @return vector of floats
std::vector<float> SampleData::get_floats(FormatType & type) {
  std::vector<float> vals;
  vals.resize(type.n_vals * n_samples);
  std::memcpy(&vals[0], &buf[type.offset], type.n_vals * n_samples * 4);
  return vals;
}

/// @brief parse strings from the buffer
/// @param type struct with number of samples per person, and buffer offset
/// @return vector of strings
std::vector<std::string> SampleData::get_strings(FormatType & type) {
  std::vector<std::string> vals;
  vals.resize(n_samples);
  std::uint32_t offset = type.offset;
  for (std::uint32_t n=0; n < n_samples; n++) {
    vals[n] = parse_string(&buf[0], offset, type.n_vals);
  }
  return vals;
}

}