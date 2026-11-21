#ifndef BCF_TYPES_H
#define BCF_TYPES_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

namespace bcf {

// Low 4-bit type codes used in BCF2 type-descriptor bytes (codes 4, 6, 8..15 are reserved).
enum Types : std::uint8_t {
  flag = 0,
  int8 = 1,
  int16 = 2,
  int32 = 3,
  float_ = 5,
  char_ = 7,
};

// Byte size of one element for each BCF atomic type code (indexed by `Types`).
static const std::uint8_t type_sizes[8] = {0, 1, 2, 4, 0, 4, 0, 1};

// BCF2 encodes float MISSING, END_OF_VECTOR, and 5 future reserved values as IEEE-754
// signaling NaNs (`0x7f800001..0x7f800007`), keeping quiet NaN (`0x7fc00000`) for literal NaN.
static const std::uint32_t BCF_FLOAT_MISSING = 0x7f800001u;
static const std::uint32_t BCF_FLOAT_VECTOR_END = 0x7f800002u;
static const std::uint32_t BCF_FLOAT_RESERVED_MAX = 0x7f800007u;

inline bool is_int_type(Types t) {
  return t == int8 || t == int16 || t == int32;
}

/// @brief check whether raw IEEE-754 float bits fall in BCF's reserved sentinel range
inline bool is_bcf_float_missing(std::uint32_t bits) {
  return bits >= BCF_FLOAT_MISSING && bits <= BCF_FLOAT_RESERVED_MAX;
}

// A position in the bytes of a BCF record. Each `read_` method checks that its value lies
// inside the buffer, decodes the value, and moves on to the byte after it.
class Cursor {
  const char * buf;
  std::size_t buf_len;
  std::uint32_t idx;
public:
  Cursor(const char * bytes, std::size_t n_bytes, std::uint32_t start = 0)
      : buf(bytes), buf_len(n_bytes), idx(start) {}

  // Offset in the buffer of the next byte to read.
  std::uint32_t offset() const { return idx; }

  // Move past `count` bytes, throwing `err_msg` if the buffer ends before them.
  void skip(std::uint64_t count, const char * err_msg) {
    if (static_cast<std::uint64_t>(idx) + count > buf_len) {
      throw std::invalid_argument(err_msg);
    }
    idx += static_cast<std::uint32_t>(count);
  }

  /// @brief read one BCF integer, keeping END_OF_VECTOR (`INT32_MIN + 1`) distinct from MISSING
  /// BCF2 reserves the lowest 8 signed values of each integer width (`INT*_MIN .. INT*_MIN + 7`):
  /// `+0` is MISSING (`INT32_MIN`), `+1` is END_OF_VECTOR, and `+2..+7` are reserved (treated as MISSING).
  std::int32_t read_raw_int(std::uint8_t type_size) {
    if (static_cast<std::size_t>(idx) + type_size > buf_len) {
      throw std::invalid_argument("buffer overflow while reading BCF integer");
    }
    std::int32_t val = 0;
    if (type_size == 1) {
      std::int8_t v = 0;
      std::memcpy(&v, buf + idx, 1);
      val = (v == INT8_MIN + 1) ? (INT32_MIN + 1)
          : (v <= INT8_MIN + 7) ? INT32_MIN
          : static_cast<std::int32_t>(v);
    } else if (type_size == 2) {
      std::int16_t v = 0;
      std::memcpy(&v, buf + idx, 2);
      val = (v == INT16_MIN + 1) ? (INT32_MIN + 1)
          : (v <= INT16_MIN + 7) ? INT32_MIN
          : static_cast<std::int32_t>(v);
    } else if (type_size == 4) {
      std::memcpy(&val, buf + idx, 4);
      if (val > INT32_MIN + 1 && val <= INT32_MIN + 7) {
        val = INT32_MIN;
      }
    } else {
      throw std::invalid_argument("invalid BCF integer width: " + std::to_string(type_size));
    }
    idx += type_size;
    return val;
  }

  /// @brief read one BCF integer when both MISSING and END_OF_VECTOR should map to `INT32_MIN`
  std::int32_t read_int(std::uint8_t type_size) {
    std::int32_t val = read_raw_int(type_size);
    return (val == INT32_MIN + 1) ? INT32_MIN : val;
  }

  /// @brief read raw 32-bit float bits, collapsing reserved non-EOV sentinels to `BCF_FLOAT_MISSING`
  std::uint32_t read_raw_float() {
    if (static_cast<std::size_t>(idx) + 4 > buf_len) {
      throw std::invalid_argument("buffer overflow while reading BCF float");
    }
    std::uint32_t bits = 0;
    std::memcpy(&bits, buf + idx, 4);
    idx += 4;
    if (bits != BCF_FLOAT_VECTOR_END && is_bcf_float_missing(bits)) {
      return BCF_FLOAT_MISSING;
    }
    return bits;
  }

  /// @brief read a fixed `size`-byte BCF string slice, stopping early at any NUL padding byte
  std::string read_string(std::uint32_t size) {
    if (static_cast<std::size_t>(idx) + size > buf_len) {
      throw std::invalid_argument("buffer overflow while reading BCF string");
    }
    const char * start = buf + idx;
    const char * nul = static_cast<const char *>(std::memchr(start, '\0', size));
    std::size_t len = (nul != nullptr) ? static_cast<std::size_t>(nul - start) : size;
    idx += size;
    return std::string(start, len);
  }

  /// @brief read one byte of a BCF type descriptor
  std::uint8_t read_byte() {
    if (idx >= buf_len) {
      throw std::invalid_argument("buffer overflow while reading BCF type descriptor");
    }
    return static_cast<std::uint8_t>(buf[idx++]);
  }

  // These two use `TypeDescriptor`, so they are defined below it.

  /// @brief read a type-tagged scalar BCF integer (used for INFO/FORMAT header dictionary IDs)
  std::int32_t read_typed_int();

  /// @brief read a type-tagged BCF string (missing `.` may be encoded as 0-len `char_` `0x07` or `0x00`)
  std::string read_typed_string();
};

// The type descriptor in front of every value in a BCF record, which says what type the
// value's elements are, and how many of them follow:
// - Bits 0..3 of its first byte hold the atomic `Types` code.
// - Bits 4..7 hold the element count (`0..14`), unless equal to `15`, which signals
//   an overflow length stored immediately afterwards as a typed BCF integer.
//
// A descriptor is filled in by its constructor, not returned from a `Cursor` method,
// because returning one by value made reading INFO fields 8% slower.
struct TypeDescriptor {
  Types type = flag;
  std::uint32_t n_vals = 0;
  std::uint8_t type_size = 0;

  /// @brief read the descriptor at `cursor`, leaving the cursor at the value's first element
  explicit TypeDescriptor(Cursor & cursor) {
    std::uint8_t byte = cursor.read_byte();
    std::uint8_t raw_type = byte & 0x0F;
    if (raw_type != flag && raw_type != int8 && raw_type != int16 &&
        raw_type != int32 && raw_type != float_ && raw_type != char_) {
      throw std::invalid_argument("invalid BCF value type: " + std::to_string(raw_type));
    }
    type = static_cast<Types>(raw_type);
    type_size = type_sizes[raw_type];
    n_vals = byte >> 4;

    if (n_vals == 15) {
      // Read the follow-up integer type byte and decode the true vector length.
      std::uint8_t len_byte = cursor.read_byte();
      std::uint8_t len_type = len_byte & 0x0F;
      if ((len_byte >> 4) != 1 || (len_type != int8 && len_type != int16 && len_type != int32)) {
        throw std::invalid_argument("cannot identify number of bytes to read");
      }
      std::int32_t count = cursor.read_raw_int(type_sizes[len_type]);
      if (count < 0) {
        throw std::invalid_argument("negative BCF vector length");
      }
      n_vals = static_cast<std::uint32_t>(count);
    }
  }
};

inline std::int32_t Cursor::read_typed_int() {
  TypeDescriptor desc(*this);
  if (!is_int_type(desc.type) || desc.n_vals != 1) {
    throw std::invalid_argument("expected scalar BCF integer");
  }
  return read_int(desc.type_size);
}

inline std::string Cursor::read_typed_string() {
  TypeDescriptor desc(*this);
  if (desc.type != char_ && !(desc.type == flag && desc.n_vals == 0)) {
    throw std::invalid_argument("expected BCF string value");
  }
  return read_string(desc.n_vals);
}

} // namespace bcf

#endif
