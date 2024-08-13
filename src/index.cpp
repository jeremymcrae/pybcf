
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

#include "gzstream.h"

#include "index.h"

namespace bcf {


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
  for (std::uint32_t i=0; i< n_ref; i++) {
    std::unordered_map<std::uint32_t, Bin> bins;
    infile.read(reinterpret_cast<char *>(&n_bins), sizeof(n_bins));
    
    for (std::uint32_t j=0; j<n_bins; j++) {
      infile.read(reinterpret_cast<char *>(&bin_idx), sizeof(bin_idx));
      infile.read(reinterpret_cast<char *>(&v_offset), sizeof(v_offset));
      infile.read(reinterpret_cast<char *>(&n_chunks), sizeof(n_chunks));

      bin_offsets = parse_virtual_offset(v_offset);
      std::vector<Chunk> chunks;
      for (std::uint32_t k=0; k<n_chunks; k++) {
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

/// calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
int IndexFile::reg2bin(std::int64_t beg, std::int64_t end) {
  int l, s = min_shift, t = ((1 << depth * 3) - 1) / 7;
  for (--end, l = depth; l > 0; --l, s += 3, t -= 1 << l * 3) {
    if (beg >> s == end >> s) {
      return t + (beg >> s);
    }
  }
  return 0;
}

/// @brief calculate the list of bins that may overlap with region [beg,end) (zero-based).
///
/// This code is from https://samtools.github.io/hts-specs/CSIv1.pdf, but adapted
/// for being inside a class.
///
/// @param beg start position of region
/// @param end end position of region
/// @return currently integer, but this should be an iterator instead
std::vector<std::uint32_t> IndexFile::reg2bins(std::int64_t beg, std::int64_t end) {
  int max_bins = bin_limit();
  
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

/* calculate maximum bin number -- valid bin numbers range within [0,bin_limit) */
int IndexFile::bin_limit() {
  return ((1 << (depth + 1) * 3) - 1) / 7;
}

Offsets IndexFile::query(std::uint32_t contig_id, std::int64_t beg) {
  // find the bins which could overlap a position
  auto bins = reg2bins(beg, beg);
  
  // cull bins which do not exist in the indexfile, and find the bin with the
  // closest start to the position
  std::int32_t bin = -1;
  for (auto & bin_idx: bins) {
    if (indices[contig_id].count(bin_idx) == 0) {
      continue;
    }
    // just use the highest bin for now, which should be the most precise
    bin = std::max(bin, (std::int32_t) bin_idx);
  }
  
  if (bin < 0) {
    throw std::out_of_range("cannot find bin including position: " + std::to_string(beg));
  }
  
  return indices[contig_id][bin].offset;
}

}

// int main() {
//   bcf::IndexFile indexfile = bcf::IndexFile("/users/jmcrae/apps/pybcf/chr5.110000001-115000000.bcf.csi");
//   std::cout << "bin: " << indexfile.reg2bin(111000000, 112000000) << std::endl;
//   indexfile.query(4, 111000000);
//   return 0;
// }

// g++ -stdlib=libc++ -std=c++11 -lz index.cpp gzstream.cpp; ./a.out
