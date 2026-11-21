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

// High-level Python return shape of an INFO entry, determined from both the BCF
// record's type byte and the tag's `Number`/`Type` metadata in the VCF header.
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

// Byte location and decoded shape of one INFO field's payload inside `Info::buf`.
struct InfoLayout {
  InfoKind kind;
  std::uint8_t type_size;  // element byte width (1, 2, or 4; 0 for zero-length flag)
  std::uint32_t offset;    // byte offset of the value payload inside `buf`
  std::uint32_t n_vals;    // number of elements encoded in the payload
};

// Associates an INFO tag's header dictionary index and metadata with its buffer slice.
struct NamedInfo {
  std::uint32_t id_idx;
  const InfoField * field;
  InfoLayout layout;
};

// Owns a variant's shared record buffer (starting at `info_offset`) and lazily
// indexes/decodes INFO fields only when queried.
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
  // If `num_info == 0`, `shared_buf` is dropped immediately so sites without INFO release memory early.
  Info(std::unique_ptr<char[]> shared_buf, std::uint32_t shared_len, std::uint32_t offset,
       const Header * hdr, std::uint32_t num_info)
      : header(hdr), n_info(num_info), info_offset(offset), buf_len(shared_len),
        buf(num_info > 0 ? std::move(shared_buf) : nullptr) {}
  Info() = default;
  Info(Info &&) = default;
  Info & operator=(Info &&) = default;

  std::size_t size();
  std::vector<NamedInfo> & get_fields();
  bool find(const char * key, std::size_t key_len, InfoLayout & out);

  std::int32_t get_int(const InfoLayout & layout) const;
  std::uint32_t get_float_bits(const InfoLayout & layout) const;
  std::string get_string(const InfoLayout & layout) const;

  std::vector<std::int32_t> get_ints(const InfoLayout & layout) const;
  std::vector<std::uint32_t> get_floats_bits(const InfoLayout & layout) const;
};

} // namespace bcf

#endif
