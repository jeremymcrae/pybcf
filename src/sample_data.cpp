
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

std::vector<std::string> SampleData::get_keys() {
  std::vector<std::string> key_names;
  for (auto & x : keys) {
    key_names.push_back(x.first);
  }
  return key_names;
}

FormatType SampleData::get_type(std::string &key) {
  if (keys.count(key) == 0) {
    throw std::invalid_argument("no entries for " + key + " in data");
  }
  return keys[key];
}

std::vector<std::int32_t> SampleData::get_ints(FormatType & type) {
  if (type.is_geno) {
    return get_geno(type);
  }
  std::vector<std::int32_t> vals;
  vals.resize(type.n_vals * n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t idx=0;
  for (std::uint32_t n=0; n < n_samples; n++) {
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      vals[idx] = parse_int(&buf[0], offset, type.type_size);
      idx++;
    }
  }
  return vals;
}

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
    __m256i missing_indicator = _mm256_set_epi32(0x81818181, 0x81818181, 0x81818181, 0x81818181,
                                                 0x81818181, 0x81818181, 0x81818181, 0x81818181);
    __m128i shuffle = _mm_set_epi8(0, 2, 4, 6, 8, 10, 12, 14, 17, 3, 5, 7, 9, 11, 13, 15);
    
    for (; n < (max_n - (max_n % 32)); n += 32) {
      initial = _mm256_loadu_si256((__m256i *) &buf[offset + n]);

      geno = _mm256_and_si256(initial, mask_geno);
      geno = _mm256_sub_epi8(_mm256_srli_epi32(geno, 1), sub);
      
      // account for missing values (due to different ploidy between samples)
      missed = _mm256_cmpeq_epi8(initial, missing_indicator);           // find missing values
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
    int8x16_t missing_indicator = vdupq_n_s8(0x80);
    
    for (; n < (max_n - (max_n % 16)); n += 16) {
      // load data from the array into SIMD registers.
      initial = vld1q_s8((std::int8_t *)&buf[offset + n]);

      geno = vandq_s8(initial, mask_geno);
      geno = vsubq_s8(vshrq_n_s8(geno, 1), sub); // shift right to remove phase bit,
                                                 // and subtract 1 to get allele

      // account for missing values (due to different ploidy between samples)
      missed = vceqq_s8(initial, missing_indicator);  // find and set missing values
      geno = vandq_s8(geno, vmvnq_s8(missed));         // erase original missing values
      geno = vorrq_s8(geno, missed);                   // swap in new missing values

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
  
  std::uint32_t missing_indicator = 1 << ((8 * type.type_size) - 1);
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

std::vector<float> SampleData::get_floats(FormatType & type) {
  std::vector<float> vals;
  vals.resize(type.n_vals * n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t idx=0;
  for (std::uint32_t n=0; n < n_samples; n++) {
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      vals[idx] = parse_float(&buf[0], offset);
      idx++;
    }
  }
  return vals;
}

std::vector<std::string> SampleData::get_strings(FormatType & type) {
  std::vector<std::string> vals;
  vals.resize(n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t idx=0;
  for (std::uint32_t n=0; n < n_samples; n++) {
    vals[n] = parse_string(&buf[0], offset, type.n_vals);
  }
  return vals;
}

}