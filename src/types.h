#ifndef BCF_TYPES_H
#define BCF_TYPES_H

#include <cstdint>

#include "gzstream/gzstream.h"

namespace bcf {

enum Types {
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
    type = Types(byte << 4 >> 4);
    type_size = type_sizes[type];
    n_vals = byte >> 4;
    if (n_vals == 15) {
      infile.read(reinterpret_cast<char *>(&n_vals), type_size);
    }
  }
};

} // namespace

#endif
