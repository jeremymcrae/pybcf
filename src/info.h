#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
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
  std::uint32_t id_idx;
  const InfoField * field;
  InfoType info;
};

class Info {
  const Header * header = nullptr;
  bool is_parsed = false;
  std::uint32_t n_info = 0;
  std::uint32_t info_offset = 0;
  std::uint32_t buf_len = 0;
  std::unique_ptr<char[]> buf;
  std::vector<NamedInfo> fields;

  void parse();
public:
  Info(std::unique_ptr<char[]> shared_buf, std::uint32_t shared_len, std::uint32_t offset,
       const Header * hdr, std::uint32_t num_info)
      : header(hdr), n_info(num_info), info_offset(offset), buf_len(shared_len),
        buf(num_info > 0 ? std::move(shared_buf) : nullptr) {}
  Info() = default;
  Info(Info &&) = default;
  Info & operator=(Info &&) = default;

  std::size_t size();
  std::uint32_t key_id_at(std::size_t idx);
  const InfoType & type_at(std::size_t idx);
  bool get_type(const char * key, std::size_t key_len, InfoType & out);

  std::int32_t get_int(const InfoType & type) const;
  std::int32_t get_raw_int(const InfoType & type) const;
  float get_float(const InfoType & type) const;
  std::uint32_t get_raw_float_bits(const InfoType & type) const;
  std::string get_string(const InfoType & type) const;

  std::vector<std::int32_t> get_ints(const InfoType & type) const;
  std::vector<float> get_floats(const InfoType & type) const;
};

} // namespace bcf

#endif
