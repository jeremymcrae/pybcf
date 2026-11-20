#ifndef BCF_TYPES_H
#define BCF_TYPES_H

#include <algorithm>
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

// class to handle the typed values in BCF files
//
// This takes in a single 8-bit value, which has two sub-components. The
// first four bits correspond to a number for the atomic type, as defined in the
// Types enum above. The other four bits correspond to the number of elements,
// except if the value=15, in which case the number is determined from the
// following byte/s in the file.
class Typed {
public:
  Types type = flag;
  std::uint32_t n_vals = 0;
  std::uint8_t type_size = 0;
  Typed() = default;
  // reading values from a char buffer
  Typed(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
    std::uint8_t byte = get_byte(buf, buf_len, idx);
    set_type(byte);
    if (n_vals == 15) {
      // determine the count from the following bytes
      byte = get_byte(buf, buf_len, idx);
      std::uint8_t next = byte & 0x0F;
      if (next == int8) {
        if (static_cast<std::size_t>(idx) + 1 > buf_len) {
          throw std::out_of_range("buffer overflow while reading BCF vector length");
        }
        std::uint8_t count = 0;
        std::memcpy(&count, buf + idx, 1);
        n_vals = count;
        idx += 1;
      } else if (next == int16) {
        if (static_cast<std::size_t>(idx) + 2 > buf_len) {
          throw std::out_of_range("buffer overflow while reading BCF vector length");
        }
        std::uint16_t count = 0;
        std::memcpy(&count, buf + idx, 2);
        n_vals = count;
        idx += 2;
      } else if (next == int32) {
        if (static_cast<std::size_t>(idx) + 4 > buf_len) {
          throw std::out_of_range("buffer overflow while reading BCF vector length");
        }
        std::uint32_t count = 0;
        std::memcpy(&count, buf + idx, 4);
        n_vals = count;
        idx += 4;
      } else {
        throw std::invalid_argument("cannot identify number of bytes to read");
      }
      if (n_vals == 0) {
        type = flag;
      }
    }
  }

  // get a single byte from a char array, and increment the index
  std::uint8_t get_byte(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
    if (idx >= buf_len) {
      throw std::out_of_range("buffer overflow while reading BCF type descriptor");
    }
    std::uint8_t byte = 0;
    std::memcpy(&byte, buf + idx, 1);
    idx += 1;
    return byte;
  }

  // determine the type (and size) from the bit values
  void set_type(std::uint8_t byte) {
    std::uint8_t raw_type = byte & 0x0F;
    if (raw_type != flag && raw_type != int8 && raw_type != int16 &&
        raw_type != int32 && raw_type != float_ && raw_type != char_) {
      throw std::invalid_argument("invalid BCF value type: " + std::to_string(raw_type));
    }
    type = static_cast<Types>(raw_type);
    type_size = type_sizes[raw_type];
    n_vals = byte >> 4;
    if (n_vals == 0) {
      type = flag;
    }
  }
};

/// @brief check whether raw IEEE-754 float bits encode BCF missing or end-of-vector
inline bool is_bcf_float_missing(std::uint32_t bits) {
  return bits == 0x7f800001u || bits == 0x7f800002u;
}

/// @brief parse a single integer value from a char array
/// @param buf char array
/// @param buf_len total length of buf in bytes
/// @param idx offset where integer data starts at
/// @param type_size number of bytes to use for the integer (1, 2 or 4)
/// @return value as 32-bit int (missing and end-of-vector mapped to INT32_MIN)
inline std::int32_t parse_int(const char * buf, std::size_t buf_len, std::uint32_t & idx, std::uint8_t type_size) {
  if (static_cast<std::size_t>(idx) + type_size > buf_len) {
    throw std::out_of_range("buffer overflow while reading BCF integer");
  }
  // There are two types of missing values: 1) the normal missing value (0x80,
  // 0x8000, 0x80000000), and an end-of-vector value (introduced in BCFv2.2:
  // 0x81, 0x8001, 0x80000001).
  std::int32_t val = 0;
  if (type_size == 1) {
    std::int8_t v = 0;
    std::memcpy(&v, buf + idx, 1);
    val = (v == INT8_MIN || v == INT8_MIN + 1) ? INT32_MIN : static_cast<std::int32_t>(v);
  } else if (type_size == 2) {
    std::int16_t v = 0;
    std::memcpy(&v, buf + idx, 2);
    val = (v == INT16_MIN || v == INT16_MIN + 1) ? INT32_MIN : static_cast<std::int32_t>(v);
  } else if (type_size == 4) {
    std::int32_t v = 0;
    std::memcpy(&v, buf + idx, 4);
    val = (v == INT32_MIN || v == INT32_MIN + 1) ? INT32_MIN : v;
  } else {
    throw std::invalid_argument("invalid BCF integer width: " + std::to_string(type_size));
  }

  idx += type_size;
  return val;
}

/// @brief parse a single float value from a char array
/// @param buf char array
/// @param buf_len total length of buf in bytes
/// @param idx offset where float data starts at
/// @return value as float (BCF signaling-NaN missing/EOV mapped to quiet NaN)
inline float parse_float(const char * buf, std::size_t buf_len, std::uint32_t & idx) {
  if (static_cast<std::size_t>(idx) + 4 > buf_len) {
    throw std::out_of_range("buffer overflow while reading BCF float");
  }
  std::uint32_t bits = 0;
  std::memcpy(&bits, buf + idx, 4);
  idx += 4;
  if (is_bcf_float_missing(bits)) {
    return std::nanf("");
  }
  float val = 0.0f;
  std::memcpy(&val, &bits, 4);
  return val;
}

/// @brief parse a string from a char array
/// @param buf char array
/// @param buf_len total length of buf in bytes
/// @param idx offset where string data starts at
/// @param size number of bytes to include in the string
/// @return value as string
inline std::string parse_string(const char * buf, std::size_t buf_len, std::uint32_t & idx, std::uint32_t size) {
  if (static_cast<std::size_t>(idx) + size > buf_len) {
    throw std::out_of_range("buffer overflow while reading BCF string");
  }
  std::string val(buf + idx, size);
  idx += size;
  val.erase(std::find(val.begin(), val.end(), '\0'), val.end());
  return val;
}

} // namespace

#endif
