
#include <iostream>

#include <cstring>

#include "sample_data.h"
#include "types.h"

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

FormatType SampleData::get_type(std::string &key) {
  if (keys.count(key) == 0) {
    throw std::invalid_argument("no entries for " + key + " in data");
  }
  return keys[key];
}

std::vector<std::int32_t> SampleData::get_ints(FormatType & type) {
  if (type.is_geno) {
    // confirm we checked sample phasing if we look at the genotype data
    phase_checked = true;
  }
  std::vector<std::int32_t> vals;
  vals.resize(type.n_vals * n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t idx=0;
  for (std::uint32_t n=0; n < n_samples; n++) {
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      vals[idx] = parse_int(&buf[0], offset, type.type_size);
      if (type.is_geno) {
        phase[n] = vals[idx] & 0x00000001;
        vals[idx] = (vals[idx] >> 1) - 1;
        // this only checks on genotype status, but this should apply to other
        // fields too (AD, DP etc), as if a sample lacks gt, other fields 
        // should also be absent
        missing[n] = vals[idx] == -1;
      }
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
  vals.resize(type.n_vals * n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t idx=0;
  for (std::uint32_t n=0; n < n_samples; n++) {
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      vals[idx] = parse_string(&buf[0], offset, type.type_size);
      idx++;
    }
  }
  return vals;
}

}