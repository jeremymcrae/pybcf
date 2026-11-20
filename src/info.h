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
  std::unordered_map<std::string, InfoType> keys;

  const Header * header = nullptr;
  bool is_parsed = false;
  std::uint32_t n_info = 0;
  std::vector<char> buf;

  void parse();

  std::vector<float> scalar_floats;                      // type 0
  std::vector<std::int32_t> scalar_ints;                 // type 1
  std::vector<std::string> strings;                      // type 2

  std::vector<std::vector<float>> vector_floats;         // type 3
  std::vector<std::vector<std::int32_t>> vector_ints;    // type 4
public:
  Info(const char * data, std::size_t len, const Header * header, std::uint32_t n_info)
      : header(header), n_info(n_info) {
    if (n_info > 0 && len > 0) {
      buf.assign(data, data + len);
    }
  }
  Info() = default;
  std::vector<std::string> get_keys();
  InfoType get_type(std::string &key);

  std::int32_t get_int(std::uint32_t offset) const { return scalar_ints.at(offset); }
  float get_float(std::uint32_t offset) const { return scalar_floats.at(offset); }
  const std::string & get_string(std::uint32_t offset) const { return strings.at(offset); }

  const std::vector<std::int32_t> & get_ints(std::uint32_t offset) const { return vector_ints.at(offset); }
  const std::vector<float> & get_floats(std::uint32_t offset) const { return vector_floats.at(offset); }
};

} // namespace

#endif
