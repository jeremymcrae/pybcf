#include <algorithm>
#include <stdexcept>

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
  // CSI stores bin numbers as uint32, which caps the tree at depth 10.
  if (min_shift < 0 || min_shift > 60 || depth < 0 || depth > 10 ||
      (min_shift + depth * 3) >= 62 || l_aux < 0) {
    throw std::invalid_argument("invalid CSI header");
  }
  // Skip optional auxiliary metadata bytes (`l_aux`), which are not needed for BCF queries.
  infile.skip(static_cast<std::size_t>(l_aux), "truncated CSI aux data");

  // 2. Read the per-contig bin maps (`n_ref` contigs in header contig-ID order).
  std::int32_t n_ref = 0;
  infile.read_exact(&n_ref, sizeof(n_ref), "invalid CSI reference count");
  if (n_ref < 0) {
    throw std::invalid_argument("invalid CSI reference count");
  }

  std::uint64_t bin_limit = static_cast<std::uint64_t>(first_bin_at_level(depth + 1));
  indices.reserve(std::min<std::int32_t>(n_ref, 4096));

  for (std::int32_t i = 0; i < n_ref; i++) {
    std::vector<Bin> bins;
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

      // Grow `chunks` as records are read rather than sizing it from `n_chunks`, so a
      // corrupt count hits end-of-file instead of allocating gigabytes up front.
      std::vector<Chunk> chunks;
      for (std::int32_t k = 0; k < n_chunks; k++) {
        Chunk chunk;
        infile.read_exact(&chunk, sizeof(chunk), "truncated CSI chunk record");
        chunks.push_back(chunk);
      }
      // Skip htslib's per-contig summary pseudo-bin (`bin_limit + 1`).
      if (bin_idx < bin_limit) {
        bins.push_back({bin_idx, loffset, std::move(chunks)});
      }
    }
    // Queries binary-search the bins by number, so sort them once here.
    std::sort(bins.begin(), bins.end(), [](const Bin & a, const Bin & b) { return a.id < b.id; });
    indices.push_back(std::move(bins));
  }
  has_index = true;
  return true;
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
  end = std::min(max_coord(), std::max(beg + 1, end));

  // `bins` is sorted by bin number. Tree level `l` splits `[0, max_coord)` into `8^l` bins
  // of width `2^(min_shift + 3 * (depth - l))`, numbered left to right from
  // `first_bin_at_level(l)`, so the stored bins of one level sit next to each other.
  const std::vector<Bin> & bins = indices[contig_id];
  auto first_bin_from = [&bins](std::int64_t id) {
    return std::lower_bound(bins.begin(), bins.end(), id,
                            [](const Bin & bin, std::int64_t value) { return bin.id < value; });
  };

  // 1. Find `min_off`, a file offset that no record overlapping `[beg, end)` can come before.
  // A bin's `loffset` is the offset of the first record overlapping the start of that bin,
  // so any stored bin starting at or before `beg` gives such a bound, and the latest one
  // gives the tightest. On each level that is the last stored bin at or to the left of
  // the bin containing `beg`.
  std::uint64_t min_off = 0;
  std::int32_t shift = min_shift + depth * 3;
  for (std::int32_t level = 0; level <= depth; ++level, shift -= 3) {
    std::int64_t level_first = first_bin_at_level(level);
    auto next = first_bin_from(level_first + (beg >> shift) + 1);
    if (next != bins.begin() && (next - 1)->id >= level_first) {
      min_off = std::max(min_off, (next - 1)->loffset);
    }
  }

  // 2. Among the stored bins overlapping `[beg, end)`, find the earliest chunk that ends
  // after `min_off`. Only stored bins are visited, so the work is bounded by the size of
  // the index, even when tiny bins make the range of overlapped bins huge.
  std::uint64_t best_offset = UINT64_MAX;
  shift = min_shift + depth * 3;
  for (std::int32_t level = 0; level <= depth; ++level, shift -= 3) {
    std::int64_t first = first_bin_at_level(level) + (beg >> shift);
    std::int64_t last = first_bin_at_level(level) + ((end - 1) >> shift);
    for (auto bin = first_bin_from(first); bin != bins.end() && bin->id <= last; ++bin) {
      for (const auto & chunk : bin->chunks) {
        if (chunk.end > min_off && chunk.begin < best_offset) {
          best_offset = chunk.begin;
        }
      }
    }
  }

  if (best_offset == UINT64_MAX) {
    return false;
  }
  // A chunk can begin before `min_off`; the records in that first part cannot overlap.
  out = parse_virtual_offset(std::max(best_offset, min_off));
  return true;
}

}
