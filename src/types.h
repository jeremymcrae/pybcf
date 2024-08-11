#ifndef BCF_TYPES_H
#define BCF_TYPES_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#include "gzstream.h"

namespace bcf {

enum Types {
  flag=0,
  int8=1,
  int16=2,
  int32=3,
  float_=5,
  char_=7,
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
  Types type;
  std::uint32_t n_vals=0;
  std::uint8_t type_size=0;
  Typed() {}
  // reading values from a char buffer
  Typed(char *buf, std::uint32_t &idx){
    std::uint8_t byte = get_byte(buf, idx);
    set_type(byte);
    if (n_vals == 15) {
      // determine the count from the following bytes
      byte = get_byte(buf, idx);
      Types next = Types(byte & 0x0F);
      if (next == int8) {
        n_vals = *reinterpret_cast<std::uint8_t *>(&buf[idx]);
      } else if (next == int16) {
        n_vals = *reinterpret_cast<std::uint16_t *>(&buf[idx]);
      } else if (next == int32) {
        n_vals = *reinterpret_cast<std::uint32_t *>(&buf[idx]);
      } else {
        throw std::invalid_argument("cannot identify number of bytes to read");
      }
      idx += type_sizes[next];
    }
  }
  // get a single byte from a char array, and increment the index
  std::uint8_t get_byte(char *buf, std::uint32_t &idx) {
    return *reinterpret_cast<std::uint8_t *>(&buf[idx++]);
  }
  
  // determine the type (and size) from the bit values
  void set_type(std::uint8_t byte) {
    type = Types(byte & 0x0F);
    type_size = type_sizes[type];
    n_vals = byte >> 4;
    if (n_vals == 0) {
      type = Types(0);
    }
  }
};

/// @brief parse a single integer value from a char array
/// @param buf char array
/// @param idx offset where integer data starts at
/// @param type_size number of bytes to use for the integer (1, 2 or 4)
/// @return value as 32-bit int
inline std::int32_t parse_int(char * buf, std::uint32_t & idx, std::uint8_t type_size) {
  std::int32_t val=0;
  // There are two types of missing values: 1) the normal missing value (0x80,
  // 0x8000, 0x80000000), and an end-of-vector value (introduced in BCFv2.2).
  // As an example of the end-of-vector value, for allele depth, samples with 
  // only reference alleles cannot have alt allele depths, so can be assigned 
  // the end-of-vector value at the alt allele positions.
  std::int32_t missing = 1 << ((type_size << 3) - 1); // 8bit: 0x80, 16bit: 0x8000 etc
  std::int32_t not_recorded = missing | 0x1; // 8bit: 0x81, 16bit: 0x8001, 32bit: 0x80000001
  if (type_size == 1) {
    val = *reinterpret_cast<std::int8_t *>(&buf[idx]) & 0x000000FF;
  } else if (type_size == 2) {
    val = *reinterpret_cast<std::int16_t *>(&buf[idx]) & 0x0000FFFF;
  } else {
    val = *reinterpret_cast<std::int32_t *>(&buf[idx]);
  }
  // handle missing data values
  if ((val == missing) || (val == not_recorded)) {
    val = 0x80000000;
  }
  
  idx += type_size;
  return val;
}

/// @brief parse a single float value from a char array
/// @param buf char array
/// @param idx offset where float data starts at
/// @return value as float
inline float parse_float(char * buf, std::uint32_t & idx) {
  float val = *reinterpret_cast<float *>(&buf[idx]);
  idx += 4;
  return val;
}

/// @brief parse a string from a char array
/// @param buf char array
/// @param idx offset where string data starts at
/// @param size number of bytes to include in the string
/// @return value as string
inline std::string parse_string(const char * buf, std::uint32_t & idx, std::uint32_t size) {
  std::string val;
  val.resize(size);
  std::memcpy(&val[0], &buf[idx], size);
  idx += size;
  val.erase(std::find(val.begin(), val.end(), '\0'), val.end());
  return val;
}

} // namespace

#endif
