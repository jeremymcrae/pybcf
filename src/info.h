#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "gzstream.h"
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
  bool is_parsed=false;
  std::uint32_t n_info=0;
  char * buf;
  
  void parse();
  
  // keep a number of data stores, for the different value types, which can be
  // one of "Integer, Float, Flag, Character, and String"
  std::vector<float> scalar_floats;                      // type 0
  std::vector<std::int32_t> scalar_ints;                 // type 1
  std::vector<std::string> strings;                      // type 2
  
  std::vector<std::vector<float>> vector_floats;         // type 3
  std::vector<std::vector<std::int32_t>> vector_ints;    // type 4
public:
  /// @brief start Info object for a variant. Doesn't parse data until required
  /// @param _buf data to parse for info fields and values. The pointer points to
  ///             the data of the info data
  /// @param _header pointer to bcf header data, so we can interpret info fields
  /// @param _n_info number of info fields
  Info(char * buf, Header * header, std::uint32_t n_info) : header(header), n_info(n_info),  buf(buf) {};
  Info() {};
  std::vector<std::string> get_keys();
  InfoType get_type(std::string &key);
  
  std::int32_t get_int(uint32_t offset) {return scalar_ints[offset];};
  float get_float(uint32_t offset) {return scalar_floats[offset];};
  std::string get_string(uint32_t offset) {return strings[offset];};
  
  std::vector<std::int32_t> get_ints(uint32_t offset) {return vector_ints[offset];};
  std::vector<float> get_floats(uint32_t offset) {return vector_floats[offset];};
};

} // namespace

#endif
