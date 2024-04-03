
#include <iostream>
#include <bitset>

#include "info.h"
#include "types.h"

namespace bcf {

static std::int32_t parse_int(igzstream & infile, std::uint8_t size) {
  std::int32_t val=0;
  infile.read(reinterpret_cast<char *>(&val), size);
  return val;
}

static float parse_float(igzstream & infile) {
  float val;
  infile.read(reinterpret_cast<char *>(&val), 4);
  return val;
}

static std::string parse_string(igzstream & infile, std::uint32_t size) {
  std::string val(size, ' ');
  infile.read(reinterpret_cast<char *>(&val[0]), size);
  return val;
}

Info::Info(igzstream & infile, Header & header, std::uint32_t n_info) {

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

  for (std::uint32_t i = 0; i < n_info; i++) {
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    id_idx = 0;
    infile.read(reinterpret_cast<char *>(&id_idx), type_val.type_size);
    key = header.info[id_idx].id;

    // std::cout << "info idx: " << id_idx << ", info key: " << key << "(" << std::bitset<8>(typing) << ")" << std::endl;

    // now parse the value
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    
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
      if (type_val.n_vals > 1) {
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
      s_val = parse_string(infile, type_val.n_vals);
      strings.push_back(s_val);
    } else {
      for (std::uint32_t i=0; i < type_val.n_vals; i++) {
        switch(type_val.type) {
          case flag:
            break;
          case float_:
            f_val = parse_float(infile);
            if (type_val.n_vals == 1) {
              scalar_floats.push_back(f_val);
            } else {
              vector_floats[idx].push_back(f_val);
            }
            break;
          default:
            i_val = parse_int(infile, type_val.type_size);
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
  if (keys.count(key) == 0) {
    throw std::invalid_argument("unknown info field: " + key);
  }
  return keys[key];
}

}