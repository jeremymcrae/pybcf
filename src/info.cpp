#include <climits>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "info.h"
#include "types.h"

namespace bcf {

std::size_t Info::size() {
  parse();
  return fields.size();
}

std::vector<NamedInfo> & Info::get_fields() {
  parse();
  return fields;
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

  const char * raw = buf.get();
  fields.reserve(n_info);
  std::uint32_t buf_idx = info_offset;
  for (std::uint32_t i = 0; i < n_info; i++) {
    std::int32_t id_idx = parse_typed_int(raw, buf_len, buf_idx);
    if (id_idx < 0 || static_cast<std::size_t>(id_idx) >= header->info.size() ||
        header->info[id_idx].id.empty()) {
      throw std::invalid_argument("unknown INFO ID index: " + std::to_string(id_idx));
    }
    std::uint32_t uid = static_cast<std::uint32_t>(id_idx);
    const InfoField * field = &header->info[uid];

    Typed type_val(raw, buf_len, buf_idx);
    std::uint64_t field_bytes = static_cast<std::uint64_t>(type_val.n_vals) * type_val.type_size;
    if (static_cast<std::uint64_t>(buf_idx) + field_bytes > buf_len) {
      throw std::invalid_argument("INFO field exceeds BCF shared buffer");
    }

    if (type_val.type == flag && type_val.n_vals != 0) {
      throw std::invalid_argument("invalid BCF INFO flag length");
    }

    // Record each INFO field's byte offset and shape now; decode actual values only when accessed.
    InfoKind kind;
    if (field->is_flag) {
      kind = info_flag;
    } else if (type_val.n_vals == 0) {
      // Zero-length non-flag values represent empty tuples for list fields or None for scalars.
      kind = field->is_list ? info_vector_int : info_missing;
    } else if (type_val.type == char_) {
      kind = field->is_list ? info_vector_string : info_string;
    } else if (type_val.type == float_) {
      kind = ((type_val.n_vals > 1) || field->is_list) ? info_vector_float : info_scalar_float;
    } else {
      kind = ((type_val.n_vals > 1) || field->is_list) ? info_vector_int : info_scalar_int;
    }

    // If a duplicate INFO tag appears in the same record, keep the last occurrence.
    InfoType info_type = {kind, type_val.type_size, buf_idx, type_val.n_vals};
    bool updated = false;
    for (auto & entry : fields) {
      if (entry.id_idx == uid) {
        entry.info = info_type;
        updated = true;
        break;
      }
    }
    if (!updated) {
      fields.push_back({uid, field, info_type});
    }
    buf_idx += static_cast<std::uint32_t>(field_bytes);
  }

  is_parsed = true;
}

bool Info::get_type(const char * key, std::size_t key_len, InfoType & out) {
  parse();
  for (const auto & entry : fields) {
    const std::string & id = entry.field->id;
    if (id.size() == key_len && std::memcmp(id.data(), key, key_len) == 0) {
      out = entry.info;
      return true;
    }
  }
  return false;
}

std::int32_t Info::get_int(const InfoType & type) const {
  std::uint32_t idx = type.offset;
  return parse_int(buf.get(), buf_len, idx, type.type_size);
}

std::uint32_t Info::get_float_bits(const InfoType & type) const {
  std::uint32_t idx = type.offset;
  std::uint32_t bits = parse_raw_float(buf.get(), buf_len, idx);
  return (bits == BCF_FLOAT_VECTOR_END) ? BCF_FLOAT_MISSING : bits;
}

std::string Info::get_string(const InfoType & type) const {
  std::uint32_t idx = type.offset;
  return parse_string(buf.get(), buf_len, idx, type.n_vals);
}

std::vector<std::int32_t> Info::get_ints(const InfoType & type) const {
  std::vector<std::int32_t> vals;
  vals.reserve(type.n_vals);
  const char * raw = buf.get();
  std::uint32_t idx = type.offset;
  for (std::uint32_t j = 0; j < type.n_vals; j++) {
    std::int32_t v = parse_raw_int(raw, buf_len, idx, type.type_size);
    if (v == INT32_MIN + 1) {
      break;
    }
    vals.push_back(v);
  }
  return vals;
}

std::vector<std::uint32_t> Info::get_floats_bits(const InfoType & type) const {
  std::vector<std::uint32_t> vals;
  vals.reserve(type.n_vals);
  const char * raw = buf.get();
  std::uint32_t idx = type.offset;
  for (std::uint32_t j = 0; j < type.n_vals; j++) {
    std::uint32_t bits = parse_raw_float(raw, buf_len, idx);
    if (bits == BCF_FLOAT_VECTOR_END) {
      break;
    }
    vals.push_back(bits);
  }
  return vals;
}

} // namespace bcf
