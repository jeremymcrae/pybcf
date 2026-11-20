#ifndef BCF_BCF_H
#define BCF_BCF_H

#include <climits>
#include <cstdint>
#include <string>
#include <vector>

#include "gzreader.h"

#include "header.h"
#include "variant.h"
#include "index.h"

namespace bcf {

class BCF {
  GzReader infile;
  IndexFile idxfile;
  std::string query_chrom="";
  std::int32_t query_start=0;
  std::int32_t query_end=INT32_MAX;
  bool empty_region=false;
public:
  BCF(std::string path, std::string index_path="");
  bool nextvar(Variant & var);
  void set_region(std::string chrom, std::int32_t start=1, std::int32_t end=INT32_MAX);
  void close();
  Header header;
};

}

#endif
