#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "header.h"

namespace bcf {

enum InfoKind : std::int8_t {
  info_flag = -1,
  info_scalar_float = 0,
  info_scalar_int = 1,
  info_string = 2,
  info_vector_float = 3,
  info_vector_int = 4,
};

struct InfoType {
  InfoKind type;
  std::uint32_t offset;
};

class Info {
  std::unordered_map<std::string, InfoType> keys;

  const Header * header = nullptr;
  bool is_parsed = false;
  std::uint32_t n_info = 0;
  std::vector<char> buf;

  void parse();

  std::vector<float> scalar_floats;
  std::vector<std::int32_t> scalar_ints;
  std::vector<std::string> strings;

  std::vector<std::vector<float>> vector_floats;
  std::vector<std::vector<std::int32_t>> vector_ints;
public:
  Info(const char * data, std::size_t len, const Header * header, std::uint32_t n_info)
      : header(header), n_info(n_info) {
    if (n_info > 0 && len > 0) {
      buf.assign(data, data + len);
    }
  }
  Info() = default;
  std::vector<std::string> get_keys();
  bool get_type(const std::string & key, InfoType & out);

  std::int32_t get_int(std::uint32_t offset) const { return scalar_ints.at(offset); }
  float get_float(std::uint32_t offset) const { return scalar_floats.at(offset); }
  const std::string & get_string(std::uint32_t offset) const { return strings.at(offset); }

  const std::vector<std::int32_t> & get_ints(std::uint32_t offset) const { return vector_ints.at(offset); }
  const std::vector<float> & get_floats(std::uint32_t offset) const { return vector_floats.at(offset); }
};

} // namespace

#endif
