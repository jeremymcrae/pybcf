#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "header.h"

namespace bcf {

enum InfoKind : std::int8_t {
  info_missing = -2,
  info_flag = -1,
  info_scalar_float = 0,
  info_scalar_int = 1,
  info_string = 2,
  info_vector_float = 3,
  info_vector_int = 4,
  info_vector_string = 5,
};

struct InfoType {
  InfoKind type;
  std::uint8_t type_size;
  std::uint32_t offset;
  std::uint32_t n_vals;
};

struct NamedInfo {
  std::string key;
  InfoType info;
};

class Info {
  const Header * header = nullptr;
  bool is_parsed = false;
  std::uint32_t n_info = 0;
  std::vector<NamedInfo> fields;
  std::vector<char> buf;

  void parse();
public:
  Info(const char * data, std::size_t len, const Header * header, std::uint32_t n_info)
      : header(header), n_info(n_info) {
    if (n_info > 0 && len > 0) {
      buf.assign(data, data + len);
    }
  }
  Info() = default;
  std::vector<std::string> get_keys();
  std::size_t size();
  bool get_type(const std::string & key, InfoType & out);

  std::int32_t get_int(const InfoType & type) const;
  float get_float(const InfoType & type) const;
  std::string get_string(const InfoType & type) const;

  std::vector<std::int32_t> get_ints(const InfoType & type) const;
  std::vector<float> get_floats(const InfoType & type) const;
};

} // namespace bcf

#endif
