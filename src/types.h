#ifndef BCF_TYPES_H
#define BCF_TYPES_H

#include <cstdint>

#include "gzstream/gzstream.h"

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
  Typed(std::uint8_t byte, igzstream & infile) {
    set_type(byte);
    if (n_vals == 15) {
      // determine the count from the following bytes
      infile.read(reinterpret_cast<char *>(&byte), 1);
      Types next = Types(byte & 0x0F);
      infile.read(reinterpret_cast<char *>(&n_vals), type_sizes[next]);
    }
  }
  
  // alternative version for reading values from a char buffer
  Typed(char *buf, std::uint32_t &idx){
    std::uint8_t byte = get_byte(buf, idx);
    set_type(byte);
    if (n_vals == 15) {
      // determine the count from the following bytes
      byte = get_byte(buf, idx);
      Types next = Types(byte & 0x0F);
      n_vals = *reinterpret_cast<std::uint8_t *>(&buf[idx]);
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

} // namespace

#endif
