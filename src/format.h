#ifndef BCF_FORMAT_H
#define BCF_FORMAT_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "gzstream/gzstream.h"
#include "header.h"
#include "types.h"

namespace bcf {

struct FormatType {
  Types type;            // which type the data is
  std::uint8_t size;     // size of individual data entries
  std::uint32_t offset;  // where the data starts in the buffer
  std::uint32_t len;     // number of entries
};

class Format {
  // for each info record, track which type it is, and some index value
  std::unordered_map<std::string, FormatType> keys;
  char * buf={};
  Header * header;
public:
  Format(igzstream & infile, Header & _header, std::uint32_t len, std::uint32_t n_fmt);
  Format() {};
  // ~Format() {
  //   delete[] buf;
  // };
  void parse(std::string & key);

  // bool get(std::string key);

  // keep a number of data stores, for the different value types, which can be
  // one of "Integer, Float, Flag, Character, and String"
  std::vector<float> scalar_floats;                      // type 0
  std::vector<std::int64_t> scalar_ints;                 // type 1
  std::vector<std::string> strings;                     // type 2
  
  std::vector<std::vector<float>> vector_floats;         // type 3
  std::vector<std::vector<std::int64_t>> vector_ints;    // type 4
};

} // namespace

#endif
