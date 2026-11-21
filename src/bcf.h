#ifndef BCF_BCF_H
#define BCF_BCF_H

#include <climits>
#include <cstdint>
#include <string>

#include "gzreader.h"
#include "header.h"
#include "index.h"
#include "variant.h"

namespace bcf {

class BCF {
  GzReader infile;
  IndexFile idxfile;
  std::int32_t query_contig_id = -1;
  std::int32_t query_start = 0;
  std::int32_t query_end = INT32_MAX;
  bool empty_region = false;
public:
  BCF(const std::string & path, const std::string & index_path = "");
  bool nextvar(Variant & var);
  void set_region(const std::string & chrom, std::int32_t start = 0, std::int32_t end = INT32_MAX);
  void close();
  Header header;
};

} // namespace bcf

#endif
