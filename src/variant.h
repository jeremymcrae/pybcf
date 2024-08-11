#ifndef BCF_VARIANT_H
#define BCF_VARIANT_H

#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>

#include "gzstream.h"

#include "header.h"
#include "types.h"
#include "info.h"
#include "sample_data.h"

namespace bcf {


class Variant {
  std::vector<char> buf;
  std::int32_t contig_idx=0;
  std::int32_t rlen=0;
  std::uint32_t n_alleles=0;
  std::uint32_t n_info=0;
  std::uint32_t n_fmt=0;
  std::uint32_t n_sample=0;
public:
  Variant(igzstream & infile, Header & header);
  Variant() {};
  
  std::string chrom="";
  std::int32_t pos=0;
  std::string ref="";
  std::vector<std::string> alts;
  float qual=std::nanf("1");
  std::string varid="";
  std::vector<std::string> filters;
  Info info;
  SampleData sample_data;
};

}

#endif
