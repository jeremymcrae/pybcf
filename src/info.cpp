
#include <iostream>
#include <bitset>
#include <vector>

#include "info.h"
#include "types.h"

namespace bcf {

std::vector<std::string> Info::get_keys() {
  if (!is_parsed) {
    parse();
    is_parsed = true;
  }
  std::vector<std::string> key_vals;
  for (auto & x : keys) {
    key_vals.push_back(x.first);
  }
  return key_vals;
}

/// @brief parse the info data for keys and values
void Info::parse() {
  std::uint32_t buf_idx = 0;
  Typed type_val;

  // read the info fields. TODO - find out a way to skip this if not required
  std::uint32_t id_idx;
  std::int8_t info_idx;
  std::uint32_t idx;
  std::string key;
  std::string number;
  bool has_multiple = false; 

  float f_val;
  std::int32_t i_val;
  std::string s_val;

  for (std::uint32_t i = 0; i < n_info; i++) {
    type_val = {buf, buf_idx};
    id_idx = parse_int(&buf[0], buf_idx, type_val.type_size);
    InfoField & field = header->info[id_idx];
    key = field.id;
    number = field.number;
    has_multiple = number == "A" || number == "R" || number == "G";
    if ((number.find_first_not_of("0123456789") == std::string::npos)) {
      has_multiple |= std::stol(number) > 1;
    }

    // now parse the value
    type_val = {&buf[0], buf_idx};
    
    // figure out which datastore to keep values in
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
      if ((type_val.n_vals > 1 || has_multiple) && (type_val.type != char_)) {
        // increase offset for vector values
        info_idx += 3;
      }
      if (type_val.type == float_) {
        if (type_val.n_vals == 1 && !has_multiple) { 
          idx = scalar_floats.size();
        } else {
          idx = vector_floats.size();
          vector_floats.push_back({});
        }
      } else if (type_val.type == char_) {
        idx = strings.size();
      } else {
        if (type_val.n_vals == 1 && !has_multiple) {
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
            if (type_val.n_vals == 1 && !has_multiple) {
              scalar_floats.push_back(f_val);
            } else {
              vector_floats[idx].push_back(f_val);
            }
            break;
          default:
            i_val = parse_int(&buf[0], buf_idx, type_val.type_size);
            if (type_val.n_vals == 1  && !has_multiple) {
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
    is_parsed = true;
  }
  
  if (keys.count(key) == 0) {
    throw std::invalid_argument("unknown info field: " + key);
  }
  return keys[key];
}

}