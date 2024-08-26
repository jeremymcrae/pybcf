#ifndef BCF_INDEX_H
#define BCF_INDEX_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bcf {

/// stores the file offset for a bzgf chunk containing the start of a bin, and
/// also the offset where the bin begins within the uncompressed chunk
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
  std::vector<std::uint32_t> region_to_bins(std::int64_t beg, std::int64_t end);
  std::uint32_t get_bin_depth(std::uint32_t bin_idx);
  std::uint32_t get_bin_offset(std::uint32_t bin_idx);
  Offsets query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end);
  bool has_index = false;
};

}

#endif
