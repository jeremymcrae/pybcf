#ifndef BCF_INDEX_H
#define BCF_INDEX_H

#include <cstdint>
#include <string>
#include <vector>

#include "gzreader.h"

namespace bcf {

// Half-open BGZF virtual-offset span `[begin, end)` of records belonging to a CSI bin.
struct Chunk {
  std::uint64_t begin;
  std::uint64_t end;
};
static_assert(sizeof(Chunk) == 16, "CSI chunk records are two packed uint64 offsets");

// One CSI octree bin: `id` is its number in the tree, `loffset` is the virtual offset of
// the first record overlapping the start of the bin (used to skip upstream chunks), and
// `chunks` are the file ranges holding the bin's records.
struct Bin {
  std::uint32_t id;
  std::uint64_t loffset;
  std::vector<Chunk> chunks;
};

// Reader and region-query engine for CSIv1 coordinate-sorted index files.
class IndexFile {
  std::int32_t min_shift = 0; // log2 width (in bp) of leaf-level bins (typically 14 -> 16 kb)
  std::int32_t depth = 0;     // octree depth (0 = root only; typically 5)
  std::vector<std::vector<Bin>> indices; // bins of each contig, sorted by `id`

  std::int64_t max_coord() const { return 1LL << (min_shift + depth * 3); }
public:
  bool open(const std::string & path, bool required = true);
  bool query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end, Offsets & out) const;
  bool has_index = false;
};

}

#endif
