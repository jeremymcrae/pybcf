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
  fields.reserve(n_info);
  Cursor cursor(buf.get(), buf_len, info_offset);
  for (std::uint32_t i = 0; i < n_info; i++) {
    // Decode the INFO tag's header dictionary index and validate it against the header.
    std::int32_t id_idx = cursor.read_typed_int();
    if (id_idx < 0 || static_cast<std::size_t>(id_idx) >= header->info.size() ||
        header->info[id_idx].id.empty()) {
      throw std::invalid_argument("unknown INFO ID index: " + std::to_string(id_idx));
    }
    std::uint32_t uid = static_cast<std::uint32_t>(id_idx);
    const InfoField * field = &header->info[uid];

    // Read the value's type/length descriptor, then step over its payload, which must
    // fit within `buf`.
    TypeDescriptor desc(cursor);
    std::uint32_t offset = cursor.offset();
    cursor.skip(static_cast<std::uint64_t>(desc.n_vals) * desc.type_size,
                "INFO field exceeds BCF shared buffer");

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

    fields.push_back({uid, field, {kind, desc.type_size, offset, desc.n_vals}});
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
  return Cursor(buf.get(), buf_len, layout.offset).read_int(layout.type_size);
}

// Return raw IEEE-754 bits so Python can distinguish BCF_FLOAT_MISSING (`None`) from literal `NaN`.
std::uint32_t Info::get_float_bits(const InfoLayout & layout) const {
  std::uint32_t bits = Cursor(buf.get(), buf_len, layout.offset).read_raw_float();
  return (bits == BCF_FLOAT_VECTOR_END) ? BCF_FLOAT_MISSING : bits;
}

std::string Info::get_string(const InfoLayout & layout) const {
  return Cursor(buf.get(), buf_len, layout.offset).read_string(layout.n_vals);
}

// Read an integer INFO vector, stopping early if an END_OF_VECTOR (`INT32_MIN + 1`) padding sentinel is hit.
std::vector<std::int32_t> Info::get_ints(const InfoLayout & layout) const {
  std::vector<std::int32_t> vals;
  vals.reserve(layout.n_vals);
  Cursor cursor(buf.get(), buf_len, layout.offset);
  for (std::uint32_t j = 0; j < layout.n_vals; j++) {
    std::int32_t v = cursor.read_raw_int(layout.type_size);
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
  Cursor cursor(buf.get(), buf_len, layout.offset);
  for (std::uint32_t j = 0; j < layout.n_vals; j++) {
    std::uint32_t bits = cursor.read_raw_float();
    if (bits == BCF_FLOAT_VECTOR_END) {
      break;
    }
    vals.push_back(bits);
  }
  return vals;
}

} // namespace bcf
