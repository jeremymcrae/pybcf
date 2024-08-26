#ifndef BCF_FORMAT_H
#define BCF_FORMAT_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "gzstream.h"
#include "header.h"
#include "types.h"

namespace bcf {

/// @brief contains details of a data type (number of values per person, data type etc)
struct FormatType {
  // which type the data is (1=int8, 2=int16, 3=int32, 5=float, 7=char)
  std::uint8_t data_type;
  // size of individual data entries
  std::uint8_t type_size;
  // where the data starts in the buffer
  std::uint32_t offset;
  // number of entries per person
  std::uint32_t n_vals; 
  // true/false for whether the type is for genotypes
  bool is_geno;
};

// class to parse data types for each sample
class SampleData {
  // for each format record, track which type it is, and some index value
  std::unordered_map<std::string, FormatType> keys;
  std::vector<char> buf;
  Header * header;
  std::vector<std::int32_t> get_geno(FormatType &type);
  std::vector<std::int32_t> parse_8bit_ints(FormatType & type);
  std::vector<std::int32_t> parse_16bit_ints(FormatType & type);
  std::vector<std::int32_t> parse_32bit_ints(FormatType & type);
public:
  SampleData(igzstream &infile, Header &_header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples);
  SampleData(){};
  std::vector<std::string> get_keys();
  FormatType get_type(std::string &key);
  std::vector<std::int32_t> get_ints(FormatType & type);
  std::vector<float> get_floats(FormatType & type);
  std::vector<std::string> get_strings(FormatType & type);

  std::uint32_t n_samples=0;
  bool phase_checked=false;
  std::vector<std::uint8_t> phase;
  std::vector<std::uint8_t> missing;
};

} // namespace

#endif
