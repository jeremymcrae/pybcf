
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#include "gzreader.h"

#include "index.h"

namespace bcf {

/// @brief parse the uncompressed and compressed parts of a virtual offset
/// @param v_offset virtual offset (contains compressed offset in highest 48 bits,
///                 and uncompressed offset in lowest 16 bits). The uncompressed
///                 offset can be 65536 at max (2^16), since that is the max size
///                 of an uncompressed bgzip block
/// @return struct containing u_offset and c_offset fields
Offsets parse_virtual_offset(std::uint64_t v_offset) {
  std::uint64_t u_offset = v_offset & 0x000000000000ffff;
  std::uint64_t c_offset = v_offset >> 16;
  return { u_offset, c_offset};
}

/// @brief load file according to https://samtools.github.io/hts-specs/CSIv1.pdf
/// @param path path to CSI index file
IndexFile::IndexFile(std::string path) {
  GzReader infile(path);
  if (!infile.is_open()) {
    throw std::invalid_argument("cannot open index file at " + path);
  }
  
  // check the file header indicates this is a CSI file
  char magic[4];
  infile.read_exact(magic, 4, "doesn't look like a CSI file");
  if (magic[0] != 'C' || magic[1] != 'S' || magic[2] != 'I' || magic[3] != 1) {
    throw std::invalid_argument("doesn't look like a CSI file");
  }
  
  infile.read_exact(&min_shift, sizeof(min_shift), "invalid CSI header");
  infile.read_exact(&depth, sizeof(depth), "invalid CSI header");
  infile.read_exact(&l_aux, sizeof(l_aux), "invalid CSI header");
  if (min_shift < 0 || depth < 0 || (min_shift + depth * 3) >= 62 || l_aux < 0) {
    throw std::invalid_argument("invalid CSI header");
  }
  if (l_aux > 0) {
    aux.resize(l_aux);
    infile.read_exact(aux.data(), l_aux, "truncated CSI aux data");
  }
  infile.read_exact(&n_ref, sizeof(n_ref), "invalid CSI reference count");
  if (n_ref < 0) {
    throw std::invalid_argument("invalid CSI reference count");
  }
  
  std::uint32_t n_bins, bin_idx;
  std::uint64_t v_offset;
  std::int32_t n_chunks;
  Offsets bin_offsets, chunk_begin, chunk_end;
  for (std::int32_t i=0; i < n_ref; i++) {
    std::unordered_map<std::uint32_t, Bin> bins;
    infile.read_exact(&n_bins, sizeof(n_bins), "truncated CSI file");
    
    for (std::uint32_t j=0; j<n_bins; j++) {
      infile.read_exact(&bin_idx, sizeof(bin_idx), "invalid CSI bin record");
      infile.read_exact(&v_offset, sizeof(v_offset), "invalid CSI bin record");
      infile.read_exact(&n_chunks, sizeof(n_chunks), "invalid CSI bin record");
      if (n_chunks < 0) {
        throw std::invalid_argument("invalid CSI bin record");
      }

      bin_offsets = parse_virtual_offset(v_offset);
      std::vector<Chunk> chunks;
      chunks.reserve(n_chunks);
      for (std::int32_t k=0; k<n_chunks; k++) {
        infile.read_exact(&v_offset, sizeof(v_offset), "truncated CSI chunk record");
        chunk_begin = parse_virtual_offset(v_offset);
        infile.read_exact(&v_offset, sizeof(v_offset), "truncated CSI chunk record");
        chunk_end = parse_virtual_offset(v_offset);
        chunks.push_back({chunk_begin, chunk_end});
      }
      bins[bin_idx] = {bin_offsets, std::move(chunks)};
    }
    indices.push_back(std::move(bins));
  }
  has_index = true;
}

/// @brief calculate the list of bins that may overlap with region [beg,end) (zero-based).
///
/// This code is from https://samtools.github.io/hts-specs/CSIv1.pdf, but adapted
/// for being inside a class.
///
/// @param beg start position of region
/// @param end end position of region
/// @return vector of candidate bin indices
std::vector<std::uint32_t> IndexFile::region_to_bins(std::int64_t beg, std::int64_t end) {
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
  std::int32_t l, s = min_shift + depth * 3;
  std::int64_t t;
  for (--end, l = 0, t = 0; l <= depth; s -= 3, t += 1LL << (l * 3), ++l) {
    std::int64_t b = t + (beg >> s), e = t + (end >> s);
    for (std::int64_t i = b; i <= e; ++i) {
      bins.push_back(static_cast<std::uint32_t>(i));
    }
  }
  return bins;
}

/// @brief find the depth for a bin
///
/// We need to know the depth of a bin in order to figure out how wide the bin
/// is, and its start and end. The width can be found via get_width().
///
/// @param bin_idx bin index
/// @return depth of the bin
std::uint32_t IndexFile::get_bin_depth(std::uint32_t bin_idx) {
  std::uint64_t factor = 8;
  std::uint64_t start = 0;
  std::uint64_t end = 0;

  std::int32_t bin_depth = depth;
  while (bin_depth >= 0) {
    if ((bin_idx >= start) && (bin_idx <= end)) {
      return static_cast<std::uint32_t>(bin_depth);
    }
    bin_depth -= 1;
    start = end + 1;
    end = start * factor;
  }
  
  throw std::invalid_argument("couldn't get depth for bin: " + std::to_string(bin_idx));
}

/// @brief find the offset of a the bin relative to the first at the same depth
/// @param bin_idx bin index
/// @return 
std::uint64_t IndexFile::get_bin_offset(std::uint32_t bin_idx) {
  std::uint64_t factor = 8;
  std::uint64_t start = 0;
  std::uint64_t end = 0;

  std::int32_t bin_depth = depth;
  while (bin_depth >= 0) {
    if ((bin_idx >= start) && (bin_idx <= end)) {
      return bin_idx - start;
    }
    bin_depth -= 1;
    start = end + 1;
    end = start * factor;
  }
  
  throw std::invalid_argument("couldn't get offset for bin: " + std::to_string(bin_idx));
}

/// @brief find the actual width of a bin
/// @param depth index depth for a bin
/// @param min_shift minimum bit shift (typically 14)
std::int64_t get_width(std::int32_t depth, std::int32_t min_shift) {
  return 1LL << (min_shift + (depth * 3));
}

/// @brief find the file offsets for the first bin to overlap a given chrom region
/// @param contig_id integer for indexing into the indices vector
/// @param beg start position of region
/// @param end end position of region
/// @return file offsets for the bin as struct with u_offset and c_offset members.
Offsets IndexFile::query(std::uint32_t contig_id, std::int64_t beg, std::int64_t end) {
  if (end < beg) {
    throw std::invalid_argument("start is after end: " + std::to_string(beg) + " > " + std::to_string(end));
  }
  if (contig_id >= indices.size()) {
    throw std::out_of_range("contig ID out of range in CSI index: " + std::to_string(contig_id));
  }
  
  // find the bins which could overlap a position
  auto bins = region_to_bins(beg, end);

  // cull bins which do not exist in the indexfile, and find the bin with the
  // closest start to the position, preferrably immediately upstream of the begin
  // position, but allow using the first downstream, if upstream bins do not exist
  std::uint32_t bin_depth;
  std::int64_t bin_width, bin_start, delta;
  std::int64_t left_delta = 1LL << (min_shift + 3 * depth);
  std::int64_t right_delta = 1LL << (min_shift + 3 * depth);
  std::int64_t left_bin = -1;
  std::int64_t right_bin = -1;
  for (auto & bin_idx: bins) {
    if (indices[contig_id].count(bin_idx) == 0) {
      continue;
    }
    
    bin_depth = get_bin_depth(bin_idx);
    bin_width = get_width(bin_depth, min_shift);
    bin_start = static_cast<std::int64_t>(get_bin_offset(bin_idx)) * bin_width;
    delta = std::abs(bin_start - beg);
    
    // find the bin that starts closest to the begin position
    if (bin_start < beg) {
      if (delta < left_delta) {
        left_bin = bin_idx;
        left_delta = delta;
      };
    } else {
      if (delta < right_delta) {
        right_bin = bin_idx;
        right_delta = delta;
      };
    }
  }
  
  if ((left_bin < 0) && (right_bin < 0)) {
    throw std::out_of_range("cannot find bin for: " + std::to_string(beg) + "-" + std::to_string(end));
  }
  
  if (left_bin >= 0) {
    return indices[contig_id][static_cast<std::uint32_t>(left_bin)].offset;
  } else {
    return indices[contig_id][static_cast<std::uint32_t>(right_bin)].offset;
  }
}

}
