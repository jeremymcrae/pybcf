#include <stdexcept>
#include <utility>
#include <vector>

#include "info.h"
#include "types.h"

namespace bcf {

std::vector<std::string> Info::get_keys() {
  parse();
  std::vector<std::string> key_vals;
  key_vals.reserve(keys.size());
  for (const auto & x : keys) {
    key_vals.push_back(x.first);
  }
  return key_vals;
}

/// @brief parse the info data for keys and values
void Info::parse() {
  if (is_parsed) {
    return;
  }
  if (n_info == 0) {
    is_parsed = true;
    return;
  }
  if (header == nullptr) {
    throw std::invalid_argument("missing BCF header for INFO");
  }

  std::uint32_t buf_idx = 0;
  for (std::uint32_t i = 0; i < n_info; i++) {
    Typed id_type(buf.data(), buf.size(), buf_idx);
    std::int32_t id_idx = parse_int(buf.data(), buf.size(), buf_idx, id_type.type_size);
    auto field_it = (id_idx >= 0) ? header->info.find(static_cast<std::uint32_t>(id_idx))
                                  : header->info.end();
    if (field_it == header->info.end()) {
      throw std::invalid_argument("unknown INFO ID index: " + std::to_string(id_idx));
    }
    const InfoField & field = field_it->second;

    Typed type_val(buf.data(), buf.size(), buf_idx);
    bool is_vector = (type_val.n_vals > 1) || field.is_list;

    switch (type_val.type) {
      case flag:
        keys[field.id] = {info_flag, 0};
        break;
      case char_:
        keys[field.id] = {field.is_list ? info_vector_string : info_string,
                          static_cast<std::uint32_t>(strings.size())};
        strings.push_back(parse_string(buf.data(), buf.size(), buf_idx, type_val.n_vals));
        break;
      case float_:
        if (!is_vector) {
          keys[field.id] = {info_scalar_float, static_cast<std::uint32_t>(scalar_floats.size())};
          scalar_floats.push_back(parse_float(buf.data(), buf.size(), buf_idx));
        } else {
          std::vector<float> vals(type_val.n_vals);
          for (std::uint32_t j = 0; j < type_val.n_vals; j++) {
            vals[j] = parse_float(buf.data(), buf.size(), buf_idx);
          }
          keys[field.id] = {info_vector_float, static_cast<std::uint32_t>(vector_floats.size())};
          vector_floats.push_back(std::move(vals));
        }
        break;
      default:
        if (!is_vector) {
          keys[field.id] = {info_scalar_int, static_cast<std::uint32_t>(scalar_ints.size())};
          scalar_ints.push_back(parse_int(buf.data(), buf.size(), buf_idx, type_val.type_size));
        } else {
          std::vector<std::int32_t> vals(type_val.n_vals);
          for (std::uint32_t j = 0; j < type_val.n_vals; j++) {
            vals[j] = parse_int(buf.data(), buf.size(), buf_idx, type_val.type_size);
          }
          keys[field.id] = {info_vector_int, static_cast<std::uint32_t>(vector_ints.size())};
          vector_ints.push_back(std::move(vals));
        }
        break;
    }
  }

  std::vector<char>().swap(buf);
  is_parsed = true;
}

bool Info::get_type(const std::string & key, InfoType & out) {
  parse();
  auto it = keys.find(key);
  if (it == keys.end()) {
    return false;
  }
  out = it->second;
  return true;
}

}
