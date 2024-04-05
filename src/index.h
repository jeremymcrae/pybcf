#ifndef BCF_INDEX_H
#define BCF_INDEX_H

#include <cstdint>
#include <string>
#include <vector>

namespace bcf {

struct Chunk {
    std::uint64_t begin;
    std::uint64_t end;
};

struct Bin {
  std::uint32_t bin;
  std::uint64_t loffset;
  std::vector<Chunk> chunks;
};

class IndexFile {
  std::int32_t min_shift;
  std::int32_t depth;
  std::int32_t l_aux;
  std::vector<std::int8_t> aux;
  std::int32_t n_ref;
  std::vector<std::vector<Bin>> indices;
public:
  IndexFile(std::string path);
  int reg2bins(size_t contig, std::int64_t beg, std::int64_t end);
  std::uint64_t query();
};

}

#endif
