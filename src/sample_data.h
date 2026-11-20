#ifndef BCF_FORMAT_H
#define BCF_FORMAT_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "gzreader.h"
#include "header.h"
#include "types.h"

namespace bcf {

/// @brief contains details of a FORMAT field in the sample buffer
struct FormatType {
  Types data_type;
  std::uint8_t type_size;
  std::uint32_t offset;
  std::uint32_t n_vals;
  bool is_geno;
};

struct NamedFormat {
  std::string key;
  FormatType format;
};

// class to parse data types for each sample
class SampleData {
  const Header * header = nullptr;
  std::uint32_t n_fmt = 0;
  bool is_parsed = false;
  std::vector<NamedFormat> fields;
  std::vector<char> buf;

  void parse();
  std::vector<std::int32_t> get_geno(const FormatType & type);
  std::vector<std::int32_t> parse_8bit_ints(const FormatType & type) const;
  std::vector<std::int32_t> parse_16bit_ints(const FormatType & type) const;
  std::vector<std::int32_t> parse_32bit_ints(const FormatType & type) const;
public:
  SampleData(GzReader & infile, const Header & header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t n_samples);
  SampleData() = default;
  std::vector<std::string> get_keys();
  std::size_t size();
  bool get_type(const std::string & key, FormatType & out);
  std::vector<std::int32_t> get_ints(const FormatType & type);
  std::vector<float> get_floats(const FormatType & type) const;
  std::vector<std::string> get_strings(const FormatType & type) const;

  std::uint32_t n_samples = 0;
  bool phase_checked = false;
  std::vector<std::uint8_t> phase;
};

} // namespace

#endif
