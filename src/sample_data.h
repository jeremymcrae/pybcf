#ifndef BCF_FORMAT_H
#define BCF_FORMAT_H

#include <cstdint>
#include <string>
#include <unordered_map>
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

// class to parse data types for each sample
class SampleData {
  // for each format record, track which type it is, and some index value
  std::unordered_map<std::string, FormatType> keys;
  std::vector<char> buf;
  std::vector<std::int32_t> get_geno(FormatType &type);
  std::vector<std::int32_t> parse_8bit_ints(FormatType & type);
  std::vector<std::int32_t> parse_16bit_ints(FormatType & type);
  std::vector<std::int32_t> parse_32bit_ints(FormatType & type);
public:
  SampleData(GzReader &infile, const Header &header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples);
  SampleData() = default;
  std::vector<std::string> get_keys();
  bool get_type(const std::string & key, FormatType & out) const;
  std::vector<std::int32_t> get_ints(FormatType & type);
  std::vector<float> get_floats(FormatType & type);
  std::vector<std::string> get_strings(FormatType & type);

  std::uint32_t n_samples = 0;
  bool phase_checked = false;
  std::vector<std::uint8_t> phase;
};

} // namespace

#endif
