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
  
  Header * header;
  std::vector<char> buf;
  bool is_parsed=false;
  std::uint32_t n_info=0;
  
  void parse();
  
  // keep a number of data stores, for the different value types, which can be
  // one of "Integer, Float, Flag, Character, and String"
  std::vector<float> scalar_floats;                      // type 0
  std::vector<std::int32_t> scalar_ints;                 // type 1
  std::vector<std::string> strings;                      // type 2
  
  std::vector<std::vector<float>> vector_floats;         // type 3
  std::vector<std::vector<std::int32_t>> vector_ints;    // type 4
public:
  Info(igzstream & infile, Header & _header, std::uint32_t info_len, std::uint32_t _n_info);
  Info() {};
  InfoType get_type(std::string &key);
  
  std::int32_t get_int(uint32_t offset) {return scalar_ints[offset];};
  float get_float(uint32_t offset) {return scalar_floats[offset];};
  std::string get_string(uint32_t offset) {return strings[offset];};
  
  std::vector<std::int32_t> get_ints(uint32_t offset) {return vector_ints[offset];};
  std::vector<float> get_floats(uint32_t offset) {return vector_floats[offset];};
};

} // namespace

#endif
