#include <algorithm>
#include <climits>
#include <stdexcept>
#include <unordered_map>

#include "gzreader.h"
#include "index.h"

namespace bcf {

/// @brief parse the uncompressed and compressed parts of a virtual offset
static Offsets parse_virtual_offset(std::uint64_t v_offset) {
  return {v_offset & 0xffffULL, v_offset >> 16};
}

/// @brief load file according to https://samtools.github.io/hts-specs/CSIv1.pdf
/// @param path path to CSI index file
IndexFile::IndexFile(std::string path) {
  GzReader infile(path);
  if (!infile.is_open()) {
    throw std::invalid_argument("cannot open index file at " + path);
  }

  char magic[4];
  infile.read_exact(magic, 4, "doesn't look like a CSI file");
  if (magic[0] != 'C' || magic[1] != 'S' || magic[2] != 'I' || magic[3] != 1) {
    throw std::invalid_argument("doesn't look like a CSI file");
  }

  std::int32_t l_aux = 0;
  infile.read_exact(&min_shift, sizeof(min_shift), "invalid CSI header");
  infile.read_exact(&depth, sizeof(depth), "invalid CSI header");
  infile.read_exact(&l_aux, sizeof(l_aux), "invalid CSI header");
  if (min_shift < 0 || depth < 0 || (min_shift + depth * 3) >= 62 || l_aux < 0) {
    throw std::invalid_argument("invalid CSI header");
  }
  if (l_aux > 0) {
    std::vector<char> aux(l_aux);
    infile.read_exact(aux.data(), l_aux, "truncated CSI aux data");
  }

  std::int32_t n_ref = 0;
  infile.read_exact(&n_ref, sizeof(n_ref), "invalid CSI reference count");
  if (n_ref < 0) {
    throw std::invalid_argument("invalid CSI reference count");
  }

  std::uint64_t bin_limit = ((1ULL << ((depth + 1) * 3)) - 1) / 7;
  indices.reserve(n_ref);

  for (std::int32_t i = 0; i < n_ref; i++) {
    std::unordered_map<std::uint32_t, Bin> bins;
    std::uint32_t n_bins = 0;
    infile.read_exact(&n_bins, sizeof(n_bins), "truncated CSI file");

    for (std::uint32_t j = 0; j < n_bins; j++) {
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
      // Skip htslib's metadata pseudo-bin (bin_limit + 1).
      if (bin_idx < bin_limit) {
        bins[bin_idx] = {loffset, std::move(chunks)};
      }
    }
    indices.push_back(std::move(bins));
  }
  has_index = true;
}

/// @brief calculate the list of bins that may overlap with region [beg,end) (zero-based).
std::vector<std::uint32_t> IndexFile::region_to_bins(std::int64_t beg, std::int64_t end) const {
  std::vector<std::uint32_t> bins;
  std::int64_t max_coord = 1LL << (min_shift + depth * 3);
  if (beg >= max_coord) {
    return bins;
  }
  if (beg < 0) {
    beg = 0;
  }
  if (end <= beg) {
    end = beg + 1;
  }
  if (end > max_coord) {
    end = max_coord;
  }
  std::int32_t s = min_shift + depth * 3;
  std::int64_t t = 0;
  --end;
  for (std::int32_t l = 0; l <= depth; s -= 3, t += 1LL << (l * 3), ++l) {
    std::int64_t b = t + (beg >> s), e = t + (end >> s);
    for (std::int64_t i = b; i <= e; ++i) {
      bins.push_back(static_cast<std::uint32_t>(i));
    }
  }
  return bins;
}

/// @brief find the virtual file offset for the first chunk overlapping [beg, end)
bool IndexFile::query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end, Offsets & out) const {
  if (end < beg) {
    throw std::invalid_argument("start is after end: " + std::to_string(beg) + " > " + std::to_string(end));
  }
  std::int64_t max_coord = 1LL << (min_shift + depth * 3);
  if (contig_id >= indices.size() || beg >= max_coord) {
    return false;
  }
  if (beg < 0) {
    beg = 0;
  }

  const auto & contig_bins = indices[contig_id];
  if (contig_bins.empty()) {
    return false;
  }

  // Walk from the leaf bin containing `beg` up to the root to find `min_off`.
  std::uint64_t min_off = 0;
  std::int64_t b = ((1LL << (depth * 3)) - 1) / 7 + (beg >> min_shift);
  for (std::int32_t d = depth; d >= 0; --d) {
    auto it = contig_bins.find(static_cast<std::uint32_t>(b));
    if (it != contig_bins.end()) {
      min_off = it->second.loffset;
      break;
    }
    b = (b - 1) >> 3;
  }

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
