#ifndef BCF_INDEX_H
#define BCF_INDEX_H

#include <climits>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bcf {

// Decoded BGZF virtual offset: `c_offset` is the byte offset of a compressed BGZF
// block in the BCF file, and `u_offset` is the byte offset inside that uncompressed block.
struct Offsets {
  std::uint64_t u_offset; // uncompressed offset (within bgzf chunk)
  std::uint64_t c_offset; // compressed offset (within overall file)
};

// Half-open BGZF virtual-offset span `[begin, end)` of records belonging to a CSI bin.
struct Chunk {
  std::uint64_t begin;
  std::uint64_t end;
};

// One CSI octree bin: `loffset` is the virtual offset of the earliest record overlapping
// the bin's subtree (used to prune upstream chunks), followed by its BGZF chunk list.
struct Bin {
  std::uint64_t loffset;
  std::vector<Chunk> chunks;
};

// Reader and region-query engine for CSIv1 coordinate-sorted index files.
class IndexFile {
  std::int32_t min_shift = 0; // log2 width (in bp) of leaf-level bins (typically 14 -> 16 kb)
  std::int32_t depth = 0;     // octree depth (0 = root only; typically 5)
  std::vector<std::unordered_map<std::uint32_t, Bin>> indices;

  std::int64_t max_coord() const { return 1LL << (min_shift + depth * 3); }
public:
  explicit IndexFile(const std::string & path) { open(path, true); }
  IndexFile() = default;
  bool open(const std::string & path, bool required = true);
  std::vector<std::uint32_t> region_to_bins(std::int64_t beg, std::int64_t end) const;
  bool query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end, Offsets & out) const;
  bool has_index = false;
};

}

#endif
