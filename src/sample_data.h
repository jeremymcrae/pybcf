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

// Byte location and layout of one FORMAT field inside a variant's sample buffer.
struct FormatLayout {
  Types data_type;         // BCF atomic type (int8/int16/int32/float_/char_)
  std::uint8_t type_size;  // byte width of a single element (1, 2, or 4)
  std::uint32_t offset;    // start byte offset of sample 0's values in `buf`
  std::uint32_t n_vals;    // number of values stored per sample (max ploidy/vector length)
  bool is_geno;            // true when this field is GT (uses allele/phase bit-packing)
};

// Associates a FORMAT tag's header dictionary index and metadata with its buffer slice.
struct NamedFormat {
  std::uint32_t id_idx;
  const FormatField * field;
  FormatLayout layout;
};

// Holds the raw per-sample ("individual") byte block of a BCF record and decodes
// requested FORMAT fields on demand directly into caller-supplied output buffers.
class SampleData {
  const Header * header = nullptr;
  std::uint32_t n_fmt = 0;
  std::uint32_t buf_len = 0;
  bool is_parsed = false;
  std::vector<NamedFormat> fields;
  std::unique_ptr<char[]> buf;

  void parse();
public:
  SampleData(GzReader & infile, const Header & header, std::uint32_t len, std::uint32_t n_fmt,
             std::uint32_t n_samples);
  SampleData() = default;
  SampleData(SampleData &&) = default;
  SampleData & operator=(SampleData &&) = default;

  std::size_t size();
  std::vector<NamedFormat> & get_fields();
  bool find(const char * key, std::size_t key_len, FormatLayout & out);

  void read_ints(const FormatLayout & layout, double * dst) const;
  void read_floats(const FormatLayout & layout, float * dst) const;
  void read_phase(const FormatLayout & layout, std::uint8_t * dst) const;
  std::vector<std::string> get_strings(const FormatLayout & layout) const;

  std::uint32_t n_samples = 0;
};

} // namespace bcf

#endif
