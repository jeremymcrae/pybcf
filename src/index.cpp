
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

#include "gzstream.h"

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
  igzstream infile(path.c_str());
  if (infile.fail()) {
    throw std::invalid_argument("cannot open index file at " + path);
  }
  
  // check the file header indicates this is a bcf file
  char magic[4];
  infile.read(&magic[0], 4);
  if (magic[0] != 'C' || magic[1] != 'S' || magic[2] != 'I' || magic[3] != 1) {
    throw std::invalid_argument("doesn't look like a CSI file");
  }
  
  infile.read(reinterpret_cast<char *>(&min_shift), sizeof(min_shift));
  infile.read(reinterpret_cast<char *>(&depth), sizeof(depth));
  infile.read(reinterpret_cast<char *>(&l_aux), sizeof(l_aux));
  infile.read(reinterpret_cast<char *>(&aux[0]), l_aux);
  infile.read(reinterpret_cast<char *>(&n_ref), sizeof(n_ref));
  
  std::uint32_t n_bins, bin_idx;
  std::uint64_t v_offset;
  std::int32_t n_chunks;
  Offsets bin_offsets, chunk_begin, chunk_end;
  for (std::int32_t i=0; i < n_ref; i++) {
    std::unordered_map<std::uint32_t, Bin> bins;
    infile.read(reinterpret_cast<char *>(&n_bins), sizeof(n_bins));
    
    for (std::uint32_t j=0; j<n_bins; j++) {
      infile.read(reinterpret_cast<char *>(&bin_idx), sizeof(bin_idx));
      infile.read(reinterpret_cast<char *>(&v_offset), sizeof(v_offset));
      infile.read(reinterpret_cast<char *>(&n_chunks), sizeof(n_chunks));

      bin_offsets = parse_virtual_offset(v_offset);
      std::vector<Chunk> chunks;
      for (std::int32_t k=0; k<n_chunks; k++) {
        infile.read(reinterpret_cast<char *>(&v_offset), sizeof(v_offset));
        chunk_begin = parse_virtual_offset(v_offset);
        infile.read(reinterpret_cast<char *>(&v_offset), sizeof(v_offset));
        chunk_end = parse_virtual_offset(v_offset);
        chunks.push_back({chunk_begin, chunk_end});
      }
      bins[bin_idx] = {bin_offsets, chunks};
    }
    indices.push_back(bins);
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
/// @return currently integer, but this should be an iterator instead
std::vector<std::uint32_t> IndexFile::region_to_bins(std::int64_t beg, std::int64_t end) {
  std::vector<std::uint32_t> bins;
  int l, t, n, s = min_shift + depth * 3;
  for (--end, l = n = t = 0; l <= depth; s -= 3, t += 1 << l * 3, ++l) {
    int b = t + (beg >> s), e = t + (end >> s), i;
    for (i = b; i <= e; ++i) {
      // I should use the vector of bins for a contig (chromosome).
      // Maybe an iterator of bins/chunks for the relevant region?
      bins.push_back(i);
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
  std::uint32_t factor = 8;
  std::uint32_t start = 0;
  std::uint32_t end = 0;

  std::uint32_t bin_depth = depth;
  while (bin_depth >= 0) {
    if ((bin_idx >= start) && (bin_idx <= end)) {
      return bin_depth;
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
std::uint32_t IndexFile::get_bin_offset(std::uint32_t bin_idx) {
  std::uint32_t factor = 8;
  std::uint32_t start = 0;
  std::uint32_t end = 0;

  std::uint32_t bin_depth = depth;
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
std::uint32_t get_width(std::int32_t depth, std::int32_t min_shift) {
  return 1 << (min_shift + (depth * 3));
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
  
  // find the bins which could overlap a position
  auto bins = region_to_bins(beg, end);

  // cull bins which do not exist in the indexfile, and find the bin with the
  // closest start to the position, preferrably immediately upstream of the begin
  // position, but allow using the first downstream, if upstream bins do not exist
  std::uint32_t bin_depth, bin_width, bin_start, delta;
  std::uint32_t left_delta = 1 << (min_shift + 3 * depth);
  std::uint32_t right_delta = 1 << (min_shift + 3 * depth);
  std::int32_t left_bin = -1;
  std::int32_t right_bin = -1;
  for (auto & bin_idx: bins) {
    if (indices[contig_id].count(bin_idx) == 0) {
      continue;
    }
    
    bin_depth = get_bin_depth(bin_idx);
    bin_width = get_width(bin_depth, min_shift);
    bin_start = get_bin_offset(bin_idx) * bin_width;
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
    return indices[contig_id][left_bin].offset;
  } else {
    return indices[contig_id][right_bin].offset;
  }
}

}
