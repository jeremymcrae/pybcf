#ifndef BCF_VARIANT_H
#define BCF_VARIANT_H

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "gzreader.h"
#include "header.h"
#include "info.h"
#include "sample_data.h"
#include "types.h"

namespace bcf {

class Variant {
public:
  Variant(GzReader & infile, const Header & header);
  Variant() = default;

  std::string chrom;
  std::int32_t pos = 0;
  std::int32_t rlen = 0;
  std::string ref;
  std::vector<std::string> alts;
  float qual = std::nanf("");
  std::string varid;
  std::vector<std::string> filters;
  Info info;
  SampleData sample_data;
};

}

#endif
