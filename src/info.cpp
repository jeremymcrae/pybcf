#include <climits>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "info.h"
#include "types.h"

namespace bcf {

std::vector<std::string> Info::get_keys() {
  parse();
  std::vector<std::string> key_names;
  key_names.reserve(fields.size());
  for (const auto & entry : fields) {
    key_names.push_back(entry.key);
  }
  return key_names;
}

std::size_t Info::size() {
  parse();
  return fields.size();
}

/// @brief scan INFO field descriptors without eagerly decoding values
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

  fields.reserve(n_info);
  std::uint32_t buf_idx = 0;
  for (std::uint32_t i = 0; i < n_info; i++) {
    std::int32_t id_idx = parse_typed_int(buf.data(), buf.size(), buf_idx);
    auto field_it = (id_idx >= 0) ? header->info.find(static_cast<std::uint32_t>(id_idx))
                                  : header->info.end();
    if (field_it == header->info.end()) {
      throw std::invalid_argument("unknown INFO ID index: " + std::to_string(id_idx));
    }
    const InfoField & field = field_it->second;

    Typed type_val(buf.data(), buf.size(), buf_idx);
    std::uint64_t field_bytes = static_cast<std::uint64_t>(type_val.n_vals) * type_val.type_size;
    if (static_cast<std::uint64_t>(buf_idx) + field_bytes > buf.size()) {
      throw std::invalid_argument("INFO field exceeds BCF shared buffer");
    }

    if (type_val.type == flag && type_val.n_vals != 0) {
      throw std::invalid_argument("invalid BCF INFO flag length");
    }

    InfoKind kind;
    if (field.is_flag) {
      kind = info_flag;
    } else if (type_val.n_vals == 0) {
      kind = field.is_list ? info_vector_int : info_missing;
    } else if (type_val.type == char_) {
      kind = field.is_list ? info_vector_string : info_string;
    } else if (type_val.type == float_) {
      kind = ((type_val.n_vals > 1) || field.is_list) ? info_vector_float : info_scalar_float;
    } else {
      kind = ((type_val.n_vals > 1) || field.is_list) ? info_vector_int : info_scalar_int;
    }

    InfoType info_type = {kind, type_val.type_size, buf_idx, type_val.n_vals};
    bool updated = false;
    for (auto & entry : fields) {
      if (entry.key == field.id) {
        entry.info = info_type;
        updated = true;
        break;
      }
    }
    if (!updated) {
      fields.push_back({field.id, info_type});
    }
    buf_idx += static_cast<std::uint32_t>(field_bytes);
  }

  is_parsed = true;
}

bool Info::get_type(const std::string & key, InfoType & out) {
  parse();
  for (const auto & entry : fields) {
    if (entry.key == key) {
      out = entry.info;
      return true;
    }
  }
  return false;
}

std::int32_t Info::get_int(const InfoType & type) const {
  std::uint32_t idx = type.offset;
  return parse_int(buf.data(), buf.size(), idx, type.type_size);
}

float Info::get_float(const InfoType & type) const {
  std::uint32_t idx = type.offset;
  return parse_float(buf.data(), buf.size(), idx);
}

std::string Info::get_string(const InfoType & type) const {
  std::uint32_t idx = type.offset;
  return parse_string(buf.data(), buf.size(), idx, type.n_vals);
}

std::vector<std::int32_t> Info::get_ints(const InfoType & type) const {
  std::vector<std::int32_t> vals;
  vals.reserve(type.n_vals);
  std::uint32_t idx = type.offset;
  for (std::uint32_t j = 0; j < type.n_vals; j++) {
    std::int32_t v = parse_raw_int(buf.data(), buf.size(), idx, type.type_size);
    if (v == INT32_MIN + 1) {
      break;
    }
    vals.push_back(v);
  }
  return vals;
}

std::vector<float> Info::get_floats(const InfoType & type) const {
  std::vector<float> vals;
  vals.reserve(type.n_vals);
  std::uint32_t idx = type.offset;
  for (std::uint32_t j = 0; j < type.n_vals; j++) {
    if (static_cast<std::size_t>(idx) + 4 > buf.size()) {
      throw std::invalid_argument("buffer overflow while reading BCF float");
    }
    std::uint32_t bits = 0;
    std::memcpy(&bits, buf.data() + idx, 4);
    idx += 4;
    if (bits == BCF_FLOAT_VECTOR_END) {
      break;
    }
    float val = std::nanf("");
    if (bits != BCF_FLOAT_MISSING) {
      std::memcpy(&val, &bits, 4);
    }
    vals.push_back(val);
  }
  return vals;
}

} // namespace bcf
