#ifndef BCF_INDEX_H
#define BCF_INDEX_H

#include <climits>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bcf {

/// stores the file offset for a bgzf chunk containing the start of a bin, and
/// also the offset where the bin begins within the uncompressed chunk
struct Offsets {
  std::uint64_t u_offset; // uncompressed offset (within bgzf chunk)
  std::uint64_t c_offset; // compressed offset (within overall file)
};

struct Chunk {
  std::uint64_t begin;
  std::uint64_t end;
};

struct Bin {
  std::uint64_t loffset;
  std::vector<Chunk> chunks;
};

class IndexFile {
  std::int32_t min_shift = 0;
  std::int32_t depth = 0;
  std::vector<std::unordered_map<std::uint32_t, Bin>> indices;
public:
  IndexFile(std::string path);
  IndexFile() = default;
  std::vector<std::uint32_t> region_to_bins(std::int64_t beg, std::int64_t end) const;
  bool query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end, Offsets & out) const;
  bool has_index = false;
};

}

#endif
