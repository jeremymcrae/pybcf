#ifndef BCF_FORMAT_H
#define BCF_FORMAT_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "gzstream/gzstream.h"
#include "header.h"
#include "types.h"

namespace bcf {

struct FormatType {
  Types type;            // which type the data is
  std::uint8_t size;     // size of individual data entries
  std::uint32_t offset;  // where the data starts in the buffer
  std::uint32_t len;     // number of entries
};

class Format {
  // for each info record, track which type it is, and some index value
  std::unordered_map<std::string, FormatType> keys;
  char * buf={};
  Header * header;
  std::uint32_t n_samples=0;
public:
  Format(igzstream & infile, Header & _header, std::uint32_t len, std::uint32_t n_fmt, std::uint32_t _n_samples);
  Format() {};
  // ~Format() {
  //   delete[] buf;
  // };
  FormatType get_type(std::string &key);
  std::vector<std::int32_t> get_ints(FormatType & type);
  std::vector<float> get_floats(FormatType & type);
  std::vector<std::string> get_strings(FormatType & type);
};

} // namespace

#endif
