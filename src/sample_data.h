#ifndef BCF_FORMAT_H
#define BCF_FORMAT_H

#include <cstddef>
#include <cstdint>
#include <memory>
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
  std::uint32_t id_idx;
  const FormatField * field;
  FormatType format;
};

// class to parse data types for each sample
class SampleData {
  const Header * header = nullptr;
  std::uint32_t n_fmt = 0;
  std::uint32_t buf_len = 0;
  bool is_parsed = false;
  std::vector<NamedFormat> fields;
  std::unique_ptr<char[]> buf;

  void parse();
  void read_geno(const FormatType & type, double * dst) const;
  void read_8bit_ints(const FormatType & type, double * dst) const;
  void read_16bit_ints(const FormatType & type, double * dst) const;
  void read_32bit_ints(const FormatType & type, double * dst) const;
public:
  SampleData(GzReader & infile, const Header & header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t n_samples);
  SampleData() = default;
  SampleData(SampleData &&) = default;
  SampleData & operator=(SampleData &&) = default;

  std::size_t size();
  std::vector<NamedFormat> & get_fields();
  bool get_type(const char * key, std::size_t key_len, FormatType & out);

  void read_ints(const FormatType & type, double * dst) const;
  void read_floats(const FormatType & type, float * dst) const;
  void read_phase(const FormatType & type, std::uint8_t * dst) const;
  std::vector<std::string> get_strings(const FormatType & type) const;

  std::uint32_t n_samples = 0;
};

} // namespace bcf

#endif
