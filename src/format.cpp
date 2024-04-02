
#include <iostream>

#include "format.h"
#include "types.h"

namespace bcf {

static std::int32_t get_int(char * buf, std::uint32_t & idx, std::uint8_t size) {
  std::int32_t val=0;
  if (size == 1) {
    val = *reinterpret_cast<std::int8_t *>(&buf[idx]);
  } else if (size == 2) {
    val = *reinterpret_cast<std::int16_t *>(&buf[idx]);
  } else {
    val = *reinterpret_cast<std::int32_t *>(&buf[idx]);
  }
  idx += size;
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

Format::Format(igzstream & infile, Header & _header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples) {
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
  for (std::uint32_t i = 0; i < n_fmt; i++ ){
    type_val = {buf, buf_idx};
    format_idx = get_int(buf, buf_idx, type_val.type_size);
    key = header->format[format_idx].id;

    type_val = {buf, buf_idx};
    keys[key] = {type_val.type, type_val.type_size, buf_idx, type_val.n_vals};
    buf_idx += (type_val.n_vals * type_val.type_size * n_samples);
  }
}

void Format::parse(std::string & key) {
  if (keys.count(key) == 0) {
    throw std::invalid_argument("no entries for " + key + " in data");
  }
  
  FormatType type = keys[key];
  std::uint32_t idx = type.offset;
  if ((type.type == int8) || (type.type == int16) || (type.type == int32)) {
    std::vector<std::int32_t> vals;
    vals.resize(type.len);
    for (std::uint32_t i = 0; i < type.len; i++) {
      vals[i] = get_int(buf, idx, type.size);
    }
  } else if (type.type == float_) {
    std::vector<float> vals;
    vals.resize(type.len);
    for (std::uint32_t i = 0; i < type.len; i++) {
      vals[i] = get_int(buf, idx, type.size);
    }
  } else {
    std::vector<std::string> vals;
    // FIXME: this isn't right, needs to account for the number of characters
    // FIXME: in each vector, rather than char size
    vals.resize(type.len);
    for (std::uint32_t i = 0; i < type.len; i++) {
      vals[i] = get_string(buf, idx, type.size);
    }
  }
}


}