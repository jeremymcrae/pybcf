#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "gzstream/gzstream.h"
#include "header.h"

namespace bcf {

struct InfoType {
  std::int8_t type;
  std::uint32_t offset;
};

class Info {
  // for each info record, track which type it is, and some index value
  std::unordered_map<std::string, InfoType> keys;
public:
  Info(igzstream & infile, Header & header, std::uint32_t n_info);
  Info() {};
  bool get(std::string key);

  // keep a number of data stores, for the different value types, which can be
  // one of "Integer, Float, Flag, Character, and String"
  std::vector<float> scalar_floats;                      // type 0
  std::vector<std::int64_t> scalar_ints;                 // type 1
  std::vector<std::string> strings;                      // type 2
  
  std::vector<std::vector<float>> vector_floats;         // type 3
  std::vector<std::vector<std::int64_t>> vector_ints;    // type 4
};

} // namespace

#endif
