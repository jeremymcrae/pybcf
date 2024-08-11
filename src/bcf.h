#ifndef BCF_BCF_H
#define BCF_BCF_H

#include <cstdint>
#include <string>
#include <vector>

#include "gzstream.h"

#include "header.h"
#include "variant.h"
#include "index.h"

namespace bcf {

class BCF {
  igzstream infile;
  IndexFile idxfile;
  std::string query_chrom="";
  std::uint32_t query_start=0;
  std::uint32_t query_end=1 << 31;
public:
  BCF(std::string path, std::string index_path="");
  Variant nextvar();
  void set_region(std::string chrom, std::uint32_t start=0, std::uint32_t end=1 << 31);
  Header header;
};

}

#endif
