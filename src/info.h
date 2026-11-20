#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "header.h"

namespace bcf {

struct InfoType {
  std::int8_t type;
  std::uint32_t offset;
};

class Info {
  // for each info record, track which type it is, and some index value
  std::unordered_map<std::string, InfoType> keys;

  Header * header = nullptr;
  bool is_parsed = false;
  std::uint32_t n_info = 0;
  std::vector<char> buf;

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
  /// @param data pointer to start of info payload
  /// @param len length of info payload in bytes
  /// @param header pointer to bcf header data, so we can interpret info fields
  /// @param n_info number of info fields
  Info(const char * data, std::size_t len, Header * header, std::uint32_t n_info)
      : header(header), n_info(n_info) {
    if (len > 0) {
      buf.assign(data, data + len);
    }
  }
  Info() = default;
  std::vector<std::string> get_keys();
  InfoType get_type(std::string &key);

  std::int32_t get_int(std::uint32_t offset) { return scalar_ints.at(offset); }
  float get_float(std::uint32_t offset) { return scalar_floats.at(offset); }
  std::string get_string(std::uint32_t offset) { return strings.at(offset); }

  std::vector<std::int32_t> get_ints(std::uint32_t offset) { return vector_ints.at(offset); }
  std::vector<float> get_floats(std::uint32_t offset) { return vector_floats.at(offset); }
};

} // namespace

#endif
