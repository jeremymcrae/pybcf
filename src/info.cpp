
#include <iostream>
#include <bitset>

#include "info.h"
#include "types.h"

namespace bcf {

static std::int32_t parse_int(char * buf, std::uint32_t & idx, std::uint8_t type_size) {
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

static float parse_float(char * buf, std::uint32_t & idx) {
  float val = *reinterpret_cast<float *>(&buf[idx]);
  idx += 4;
  return val;
}

static std::string parse_string(const char * buf, std::uint32_t & idx, std::uint32_t size) {
  std::string val;
  val.resize(size);
  std::memcpy(&val[0], &buf[idx], size);
  idx += size;
  return val;
}

Info::Info(igzstream & infile, Header & _header, std::uint32_t n_info) {
  header = &_header;
  
  // read the sample data into a buffer, but don't parse until required
  buf.resize(n_info);
  infile.read(reinterpret_cast<char *>(&buf[0]), n_info);
}

void Info::parse() {
  std::uint32_t buf_idx = 0;
  std::uint8_t typing;
  Typed type_val;

  // read the info fields. TODO - find out a way to skip this if not required
  std::uint32_t id_idx;
  std::int8_t info_idx;
  std::uint32_t idx;
  std::string key;

  float f_val;
  std::int32_t i_val;
  std::string s_val;

  for (std::uint32_t i = 0; i < buf.size(); i++) {
    type_val = {&buf[0], buf_idx};
    id_idx = parse_int(&buf[0], buf_idx, type_val.type_size);
    key = header->info[id_idx].id;

    // now parse the value
    type_val = {&buf[0], buf_idx};
    
    // figure out whuch datastore to keep values in
    switch (type_val.type) {
      case flag:
        info_idx = -1;
        idx = 0;
        break;
      case float_:
        info_idx = 0;
        break;
      case char_:
        info_idx = 2;
        break;
      default:
        info_idx = 1;
        break;
    }
    
    if (type_val.type != flag) {
      if ((type_val.n_vals > 1) && (type_val.type != char_)) {
        // increase offset for vector values
        info_idx += 3;
      }
      if (type_val.type == float_) {
        if (type_val.n_vals == 1) { 
          idx = scalar_floats.size();
        } else {
          idx = vector_floats.size();
          vector_floats.push_back({});
        }
      } else if (type_val.type == char_) {
        idx = strings.size();
      } else {
        if (type_val.n_vals == 1) {
          idx = scalar_ints.size();
        } else {
          idx = vector_ints.size();
          vector_ints.push_back({});
        }
      }
    }
    
    keys[key] = {info_idx, idx};
    
    if (type_val.type == char_) {
      s_val = parse_string(&buf[0], buf_idx, type_val.n_vals);
      strings.push_back(s_val);
    } else {
      for (std::uint32_t i=0; i < type_val.n_vals; i++) {
        switch(type_val.type) {
          case flag:
            break;
          case float_:
            f_val = parse_float(&buf[0], buf_idx);
            if (type_val.n_vals == 1) {
              scalar_floats.push_back(f_val);
            } else {
              vector_floats[idx].push_back(f_val);
            }
            break;
          default:
            i_val = parse_int(&buf[0], buf_idx, type_val.type_size);
            if (type_val.n_vals == 1) {
              scalar_ints.push_back(i_val);
            } else {
              vector_ints[idx].push_back(i_val);
            }
            break;
        }
      }
    }
  }
}

InfoType Info::get_type(std::string &key) {
  if (!is_parsed) {
    parse();
  }
  
  if (keys.count(key) == 0) {
    throw std::invalid_argument("unknown info field: " + key);
  }
  return keys[key];
}

}