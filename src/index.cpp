#include <algorithm>
#include <climits>
#include <stdexcept>
#include <unordered_map>

#include "gzreader.h"
#include "index.h"

namespace bcf {

namespace {

/// @brief split a 64-bit BGZF virtual offset into intra-block (low 16 bits) and file (high 48 bits) offsets
Offsets parse_virtual_offset(std::uint64_t v_offset) {
  return {v_offset & 0xffffULL, v_offset >> 16};
}

/// @brief index of the first CSI 8-ary tree bin at `level` (`sum_{k=0..level-1} 8^k = (8^level - 1) / 7`)
std::int64_t first_bin_at_level(std::int32_t level) {
  return ((1LL << (level * 3)) - 1) / 7;
}

/// @brief parent bin index in the 8-ary CSI binning tree
std::int64_t parent_bin(std::int64_t bin) {
  return (bin - 1) >> 3;
}

} // namespace

/// @brief load a BGZF-compressed CSIv1 index file (https://samtools.github.io/hts-specs/CSIv1.pdf)
/// @param path path to CSI index file
/// @param required whether a missing index file should throw an error
bool IndexFile::open(const std::string & path, bool required) {
  has_index = false;
  indices.clear();

  GzReader infile(path);
  if (!infile.is_open()) {
    if (required) {
      throw std::invalid_argument("cannot open index file at " + path);
    }
    return false;
  }

  // 1. Validate the 4-byte "CSI\x01" magic header and global binning parameters.
  char magic[4];
  infile.read_exact(magic, 4, "doesn't look like a CSI file");
  if (magic[0] != 'C' || magic[1] != 'S' || magic[2] != 'I' || magic[3] != 1) {
    throw std::invalid_argument("doesn't look like a CSI file");
  }

  std::int32_t l_aux = 0;
  infile.read_exact(&min_shift, sizeof(min_shift), "invalid CSI header");
  infile.read_exact(&depth, sizeof(depth), "invalid CSI header");
  infile.read_exact(&l_aux, sizeof(l_aux), "invalid CSI header");
  if (min_shift < 0 || min_shift > 60 || depth < 0 || depth > 19 || (min_shift + depth * 3) >= 62 || l_aux < 0) {
    throw std::invalid_argument("invalid CSI header");
  }
  // Skip optional auxiliary metadata bytes (`l_aux`), which are not needed for BCF queries.
  if (l_aux > 0) {
    std::vector<char> aux(l_aux);
    infile.read_exact(aux.data(), l_aux, "truncated CSI aux data");
  }

  // 2. Read the per-contig bin maps (`n_ref` contigs in header contig-ID order).
  std::int32_t n_ref = 0;
  infile.read_exact(&n_ref, sizeof(n_ref), "invalid CSI reference count");
  if (n_ref < 0) {
    throw std::invalid_argument("invalid CSI reference count");
  }

  std::uint64_t bin_limit = static_cast<std::uint64_t>(first_bin_at_level(depth + 1));
  indices.reserve(std::min<std::int32_t>(n_ref, 4096));

  for (std::int32_t i = 0; i < n_ref; i++) {
    std::unordered_map<std::uint32_t, Bin> bins;
    std::int32_t n_bins = 0;
    infile.read_exact(&n_bins, sizeof(n_bins), "truncated CSI file");
    if (n_bins < 0) {
      throw std::invalid_argument("invalid CSI bin count");
    }

    for (std::int32_t j = 0; j < n_bins; j++) {
      std::uint32_t bin_idx = 0;
      std::uint64_t loffset = 0;
      std::int32_t n_chunks = 0;
      infile.read_exact(&bin_idx, sizeof(bin_idx), "invalid CSI bin record");
      infile.read_exact(&loffset, sizeof(loffset), "invalid CSI bin record");
      infile.read_exact(&n_chunks, sizeof(n_chunks), "invalid CSI bin record");
      if (n_chunks < 0) {
        throw std::invalid_argument("invalid CSI bin record");
      }

      std::vector<Chunk> chunks(n_chunks);
      for (std::int32_t k = 0; k < n_chunks; k++) {
        infile.read_exact(&chunks[k].begin, sizeof(std::uint64_t), "truncated CSI chunk record");
        infile.read_exact(&chunks[k].end, sizeof(std::uint64_t), "truncated CSI chunk record");
      }
      // Skip htslib's per-contig summary pseudo-bin (`bin_limit + 1`).
      if (bin_idx < bin_limit) {
        bins[bin_idx] = {loffset, std::move(chunks)};
      }
    }
    indices.push_back(std::move(bins));
  }
  has_index = true;
  return true;
}

/// @brief enumerate every octree bin (from root level 0 down to `depth`) overlapping 0-based `[beg, end)`
// Each tree level `l` divides `[0, max_coord)` into `8^l` bins of width `2^(min_shift + 3*(depth - l))`.
std::vector<std::uint32_t> IndexFile::region_to_bins(std::int64_t beg, std::int64_t end) const {
  std::vector<std::uint32_t> bins;
  std::int64_t limit = max_coord();
  if (beg >= limit) {
    return bins;
  }
  beg = std::max<std::int64_t>(0, beg);
  end = std::min(limit, std::max(beg + 1, end));

  std::int64_t last_pos = end - 1;
  std::int32_t shift = min_shift + depth * 3;
  for (std::int32_t level = 0; level <= depth; ++level, shift -= 3) {
    std::int64_t level_offset = first_bin_at_level(level);
    std::int64_t first_bin = level_offset + (beg >> shift);
    std::int64_t last_bin = level_offset + (last_pos >> shift);
    for (std::int64_t bin = first_bin; bin <= last_bin; ++bin) {
      bins.push_back(static_cast<std::uint32_t>(bin));
    }
  }
  return bins;
}

/// @brief find the earliest BGZF virtual offset of a chunk that can overlap 0-based `[beg, end)`
bool IndexFile::query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end, Offsets & out) const {
  if (end < beg) {
    throw std::invalid_argument("start is after end: " + std::to_string(beg) + " > " + std::to_string(end));
  }
  if (contig_id >= indices.size() || beg >= max_coord()) {
    return false;
  }
  beg = std::max<std::int64_t>(0, beg);

  const auto & contig_bins = indices[contig_id];
  if (contig_bins.empty()) {
    return false;
  }

  // Walk from the leaf bin containing `beg` up toward the root to find `min_off`:
  // any chunk ending at or before `min_off` lies strictly upstream of `beg` and can be ignored,
  // even if it belongs to a wide ancestor bin that spans `[beg, end)`.
  std::uint64_t min_off = 0;
  std::int64_t bin = first_bin_at_level(depth) + (beg >> min_shift);
  for (std::int32_t level = depth; level >= 0; --level) {
    auto it = contig_bins.find(static_cast<std::uint32_t>(bin));
    if (it != contig_bins.end()) {
      min_off = it->second.loffset;
      break;
    }
    bin = parent_bin(bin);
  }

  // Among all candidate bins overlapping `[beg, end)`, pick the smallest `chunk.begin`
  // whose `chunk.end` extends past `min_off`.
  std::uint64_t best_offset = UINT64_MAX;
  for (std::uint32_t bin_idx : region_to_bins(beg, end)) {
    auto it = contig_bins.find(bin_idx);
    if (it == contig_bins.end()) {
      continue;
    }
    for (const auto & chunk : it->second.chunks) {
      if (chunk.end > min_off && chunk.begin < best_offset) {
        best_offset = chunk.begin;
      }
    }
  }

  if (best_offset == UINT64_MAX) {
    return false;
  }
  out = parse_virtual_offset(best_offset);
  return true;
}

}
