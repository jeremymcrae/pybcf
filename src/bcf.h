#ifndef BCF_BCF_H
#define BCF_BCF_H

#include <cstdint>
#include <string>

#include "gzreader.h"
#include "header.h"
#include "index.h"
#include "variant.h"

namespace bcf {

// Top-level BCF2.2 reader supporting sequential iteration and `.csi`-indexed region queries.
class BCF {
  GzReader infile;
  IndexFile idxfile;
  Region region;           // which records to read: all of them, until `set_region` is called
  bool exhausted = false;  // true once EOF or end of the active region is reached
public:
  BCF(const std::string & path, const std::string & index_path = "");
  bool nextvar(Variant & var);
  void set_region(const std::string & chrom, std::int32_t start = 0, std::int32_t end = INT32_MAX);
  void close();
  Header header;
};

} // namespace bcf

#endif
