#ifndef BCF_TYPES_H
#define BCF_TYPES_H

#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace bcf {

enum Types : std::uint8_t {
  flag = 0,
  int8 = 1,
  int16 = 2,
  int32 = 3,
  float_ = 5,
  char_ = 7,
};

static const std::uint8_t type_sizes[8] = {0, 1, 2, 4, 0, 4, 0, 1};

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

/// @brief parse a single BCF integer, mapping reserved values to INT32_MIN (missing) and preserving INT32_MIN + 1 (end-of-vector)
inline std::int32_t parse_raw_int(const char * buf, std::size_t buf_len, std::uint32_t & idx, std::uint8_t type_size) {
  if (static_cast<std::size_t>(idx) + type_size > buf_len) {
    throw std::invalid_argument("buffer overflow while reading BCF integer");
  }
  std::int32_t val = 0;
  if (type_size == 1) {
    std::int8_t v = 0;
    std::memcpy(&v, buf + idx, 1);
    val = (v == INT8_MIN + 1) ? (INT32_MIN + 1) : (v <= INT8_MIN + 7) ? INT32_MIN : static_cast<std::int32_t>(v);
  } else if (type_size == 2) {
    std::int16_t v = 0;
    std::memcpy(&v, buf + idx, 2);
    val = (v == INT16_MIN + 1) ? (INT32_MIN + 1) : (v <= INT16_MIN + 7) ? INT32_MIN : static_cast<std::int32_t>(v);
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

/// @brief parse a single integer value from a char array (missing, reserved, and end-of-vector mapped to INT32_MIN)
inline std::int32_t parse_int(const char * buf, std::size_t buf_len, std::uint32_t & idx, std::uint8_t type_size) {
  std::int32_t val = parse_raw_int(buf, buf_len, idx, type_size);
  return (val == INT32_MIN + 1) ? INT32_MIN : val;
}

/// @brief read raw 32-bit float bits, mapping reserved non-EOV sentinels to BCF_FLOAT_MISSING
inline std::uint32_t parse_raw_float(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
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

/// @brief parse a NUL-padded or fixed-length BCF string from a char array
inline std::string parse_string(const char * buf, std::size_t buf_len, std::uint32_t & idx, std::uint32_t size) {
  if (static_cast<std::size_t>(idx) + size > buf_len) {
    throw std::invalid_argument("buffer overflow while reading BCF string");
  }
  const char * start = buf + idx;
  const char * nul = static_cast<const char *>(std::memchr(start, '\0', size));
  std::size_t len = (nul != nullptr) ? static_cast<std::size_t>(nul - start) : size;
  idx += size;
  return std::string(start, len);
}

// class to handle the typed values in BCF files
//
// This takes in a single 8-bit value, which has two sub-components. The
// first four bits correspond to a number for the atomic type, as defined in the
// Types enum above. The other four bits correspond to the number of elements,
// except if the value=15, in which case the number is determined from the
// following typed integer in the buffer.
class Typed {
public:
  Types type = flag;
  std::uint32_t n_vals = 0;
  std::uint8_t type_size = 0;
  Typed() = default;

  Typed(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
    set_type(get_byte(buf, buf_len, idx));
    if (n_vals == 15) {
      std::uint8_t len_byte = get_byte(buf, buf_len, idx);
      std::uint8_t len_type = len_byte & 0x0F;
      if ((len_byte >> 4) > 1 || (len_type != int8 && len_type != int16 && len_type != int32)) {
        throw std::invalid_argument("cannot identify number of bytes to read");
      }
      std::int32_t count = parse_raw_int(buf, buf_len, idx, type_sizes[len_type]);
      if (count < 0) {
        throw std::invalid_argument("negative BCF vector length");
      }
      n_vals = static_cast<std::uint32_t>(count);
    }
  }

  std::uint8_t get_byte(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
    if (idx >= buf_len) {
      throw std::invalid_argument("buffer overflow while reading BCF type descriptor");
    }
    std::uint8_t byte = 0;
    std::memcpy(&byte, buf + idx, 1);
    idx += 1;
    return byte;
  }

  void set_type(std::uint8_t byte) {
    std::uint8_t raw_type = byte & 0x0F;
    if (raw_type != flag && raw_type != int8 && raw_type != int16 &&
        raw_type != int32 && raw_type != float_ && raw_type != char_) {
      throw std::invalid_argument("invalid BCF value type: " + std::to_string(raw_type));
    }
    type = static_cast<Types>(raw_type);
    type_size = type_sizes[raw_type];
    n_vals = byte >> 4;
  }
};

/// @brief read a type-tagged scalar BCF integer (used for INFO/FORMAT header dictionary IDs)
inline std::int32_t parse_typed_int(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
  Typed type_val(buf, buf_len, idx);
  if (!is_int_type(type_val.type) || type_val.n_vals != 1) {
    throw std::invalid_argument("expected scalar BCF integer");
  }
  return parse_int(buf, buf_len, idx, type_val.type_size);
}

/// @brief read a type-tagged BCF string (used for ID, REF, and ALT alleles)
inline std::string parse_typed_string(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
  Typed type_val(buf, buf_len, idx);
  if (type_val.type != char_ && !(type_val.type == flag && type_val.n_vals == 0)) {
    throw std::invalid_argument("expected BCF string value");
  }
  return parse_string(buf, buf_len, idx, type_val.n_vals);
}

} // namespace bcf

#endif
