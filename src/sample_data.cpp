
#include <iostream>

#include "sample_data.h"
#include "types.h"

namespace bcf {

static std::int32_t get_int(char * buf, std::uint32_t & idx, std::uint8_t type_size) {
  std::int32_t val=0;
  if (type_size == 1) {
    val = *reinterpret_cast<std::int8_t *>(&buf[idx]) & 0x000000FF;
    if (val == 0x80) { val = 0x80000000; }  // handle missing data value
  } else if (type_size == 2) {
    val = *reinterpret_cast<std::int16_t *>(&buf[idx]) & 0x0000FFFF;
    if (val == 0x8000) { val = 0x80000000; }  // handle missing data value
  } else {
    val = *reinterpret_cast<std::int32_t *>(&buf[idx]);
  }
  idx += type_size;
  return val;
}

static float get_float(char * buf, std::uint32_t & idx) {
  float val = *reinterpret_cast<float *>(&buf[idx]);
  idx += 4;
  return val;
}

static std::string get_string(const char * buf, std::uint32_t & idx, std::uint32_t size) {
  std::string val;
  val.resize(size);
  std::memcpy(&val[0], &buf[idx], size);
  idx += size;
  return val;
}

SampleData::SampleData(igzstream & infile, Header & _header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples) {
  n_samples = _n_samples;
  header = &_header;
  if (len == 0) {
    return;
  }
  
  // read the sample data into a buffer, but don't parse until required
  buf = new char[len];
  infile.read(reinterpret_cast<char *>(&buf[0]), len);
  
  // read the available keys
  std::uint32_t buf_idx=0;
  std::uint32_t format_idx=0;
  std::string key;
  Typed type_val;
  bool is_geno;
  for (std::uint32_t i = 0; i < n_fmt; i++ ){
    type_val = {buf, buf_idx};
    format_idx = get_int(buf, buf_idx, type_val.type_size);
    key = header->format[format_idx].id;
    is_geno = key == "GT";

    type_val = {buf, buf_idx};
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
  std::vector<std::int32_t> vals;
  vals.resize(type.n_vals * n_samples);
  std::uint32_t offset = type.offset;
  std::uint32_t idx=0;
  for (std::uint32_t n=0; n < n_samples; n++) {
    for (std::uint32_t i = 0; i < type.n_vals; i++) {
      vals[idx] = get_int(buf, offset, type.type_size);
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
      vals[idx] = get_float(buf, offset);
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
      vals[idx] = get_string(buf, offset, type.type_size);
      idx++;
    }
  }
  return vals;
}

}