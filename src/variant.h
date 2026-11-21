#ifndef BCF_VARIANT_H
#define BCF_VARIANT_H

#include <climits>
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
  Variant() = default;
  Variant(Variant &&) = default;
  Variant & operator=(Variant &&) = default;

  bool read(GzReader & infile, const Header & header);
  bool read(GzReader & infile, const Header & header, std::int32_t query_contig_id,
            std::int32_t query_start, std::int32_t query_end);

  std::uint32_t chrom_id = 0;
  std::int32_t pos = 0;
  std::int32_t rlen = 0;
  std::string ref;
  std::vector<std::string> alts;
  std::uint32_t qual_bits = BCF_FLOAT_MISSING;
  std::string varid;
  std::vector<std::uint32_t> filter_ids;
  Info info;
  SampleData sample_data;
};

}

#endif
