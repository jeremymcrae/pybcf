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

// Scan the `n_info` key/value pairs in `buf` starting at `info_offset`.
// Each entry begins with a typed integer dictionary ID followed by a typed value
// descriptor and `n_vals * type_size` payload bytes. Values themselves are not
// decoded here — only their offsets and types are recorded.
void Info::parse() {
  if (is_parsed) {
    return;
  }
  if (n_info == 0) {
    is_parsed = true;
    return;
  }
  const char * raw = buf.get();
  fields.reserve(n_info);
  std::uint32_t buf_idx = info_offset;
  for (std::uint32_t i = 0; i < n_info; i++) {
    // Decode the INFO tag's header dictionary index and validate it against the header.
    std::int32_t id_idx = parse_typed_int(raw, buf_len, buf_idx);
    if (id_idx < 0 || static_cast<std::size_t>(id_idx) >= header->info.size() ||
        header->info[id_idx].id.empty()) {
      throw std::invalid_argument("unknown INFO ID index: " + std::to_string(id_idx));
    }
    std::uint32_t uid = static_cast<std::uint32_t>(id_idx);
    const InfoField * field = &header->info[uid];

    // Read the value's type/length descriptor and ensure its payload fits within `buf`.
    TypeDescriptor desc(raw, buf_len, buf_idx);
    std::uint64_t field_bytes = static_cast<std::uint64_t>(desc.n_vals) * desc.type_size;
    if (static_cast<std::uint64_t>(buf_idx) + field_bytes > buf_len) {
      throw std::invalid_argument("INFO field exceeds BCF shared buffer");
    }

    if (desc.type == flag && desc.n_vals != 0) {
      throw std::invalid_argument("invalid BCF INFO flag length");
    }

    // Combine the header's `Number`/`Type` metadata with the record's element count:
    // list fields (`Number` != 0, 1) always decode as Python tuples even when `n_vals <= 1`.
    InfoKind kind;
    if (field->is_flag) {
      kind = info_flag;
    } else if (desc.n_vals == 0) {
      // Zero-length non-flag values represent empty tuples `()` for list fields or `None` for scalars.
      kind = field->is_list ? info_vector_int : info_missing;
    } else if (desc.type == char_) {
      // BCF encodes string lists as a single comma-separated `char_` string.
      kind = field->is_list ? info_vector_string : info_string;
    } else if (desc.type == float_) {
      kind = ((desc.n_vals > 1) || field->is_list) ? info_vector_float : info_scalar_float;
    } else {
      kind = ((desc.n_vals > 1) || field->is_list) ? info_vector_int : info_scalar_int;
    }

    fields.push_back({uid, field, {kind, desc.type_size, buf_idx, desc.n_vals}});
    buf_idx += static_cast<std::uint32_t>(field_bytes);
  }

  is_parsed = true;
}

// Look up an INFO field by name without allocating a temporary std::string. If a malformed
// record repeats a key, the first copy wins, matching htslib.
bool Info::find(const char * key, std::size_t key_len, InfoLayout & out) {
  parse();
  for (const auto & entry : fields) {
    const std::string & id = entry.field->id;
    if (id.size() == key_len && std::memcmp(id.data(), key, key_len) == 0) {
      out = entry.layout;
      return true;
    }
  }
  return false;
}

std::int32_t Info::get_int(const InfoLayout & layout) const {
  std::uint32_t idx = layout.offset;
  return parse_int(buf.get(), buf_len, idx, layout.type_size);
}

// Return raw IEEE-754 bits so Python can distinguish BCF_FLOAT_MISSING (`None`) from literal `NaN`.
std::uint32_t Info::get_float_bits(const InfoLayout & layout) const {
  std::uint32_t idx = layout.offset;
  std::uint32_t bits = parse_raw_float(buf.get(), buf_len, idx);
  return (bits == BCF_FLOAT_VECTOR_END) ? BCF_FLOAT_MISSING : bits;
}

std::string Info::get_string(const InfoLayout & layout) const {
  std::uint32_t idx = layout.offset;
  return parse_string(buf.get(), buf_len, idx, layout.n_vals);
}

// Read an integer INFO vector, stopping early if an END_OF_VECTOR (`INT32_MIN + 1`) padding sentinel is hit.
std::vector<std::int32_t> Info::get_ints(const InfoLayout & layout) const {
  std::vector<std::int32_t> vals;
  vals.reserve(layout.n_vals);
  const char * raw = buf.get();
  std::uint32_t idx = layout.offset;
  for (std::uint32_t j = 0; j < layout.n_vals; j++) {
    std::int32_t v = parse_raw_int(raw, buf_len, idx, layout.type_size);
    if (v == INT32_MIN + 1) {
      break;
    }
    vals.push_back(v);
  }
  return vals;
}

// Read a float INFO vector as raw 32-bit words, stopping early at `BCF_FLOAT_VECTOR_END`.
std::vector<std::uint32_t> Info::get_floats_bits(const InfoLayout & layout) const {
  std::vector<std::uint32_t> vals;
  vals.reserve(layout.n_vals);
  const char * raw = buf.get();
  std::uint32_t idx = layout.offset;
  for (std::uint32_t j = 0; j < layout.n_vals; j++) {
    std::uint32_t bits = parse_raw_float(raw, buf_len, idx);
    if (bits == BCF_FLOAT_VECTOR_END) {
      break;
    }
    vals.push_back(bits);
  }
  return vals;
}

} // namespace bcf
