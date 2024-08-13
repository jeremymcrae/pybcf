#ifndef BCF_INDEX_H
#define BCF_INDEX_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bcf {

struct Offsets {
  std::uint64_t u_offset; // uncompressed offset (within bgzf chunk)
  std::uint64_t c_offset; // compressed offset (within overall file)
};

struct Chunk {
    Offsets begin;
    Offsets end;
};

struct Bin {
  Offsets offset;
  std::vector<Chunk> chunks;
};

class IndexFile {
  std::int32_t min_shift;
  std::int32_t depth;
  std::int32_t l_aux;
  std::vector<std::int8_t> aux;
  std::int32_t n_ref;
  std::vector<std::unordered_map<std::uint32_t, Bin>> indices;
public:
  IndexFile(std::string path);
  IndexFile() {};
  int reg2bin(std::int64_t beg, std::int64_t end);
  std::vector<std::uint32_t> reg2bins(std::int64_t beg, std::int64_t end);
  int bin_limit();
  Offsets query(std::uint32_t contig_id, std::int64_t beg);
  bool has_index = false;
};

}

#endif
