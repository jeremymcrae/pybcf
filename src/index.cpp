
#include <fstream>
#include <stdexcept>

#include <iostream>

#include "gzstream.h"

#include "index.h"

namespace bcf {


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
  
  // std::cout << "min_shift: " << min_shift << ", depth: " << depth << std::endl;
  
  std::uint32_t n_bins, bin_idx;
  std::uint64_t loffset, start, end;
  std::int32_t n_chunks;
  for (std::uint32_t i=0; i< n_ref; i++) {
    std::vector<Bin> bins;
    infile.read(reinterpret_cast<char *>(&n_bins), sizeof(n_bins));
    // if (n_bins > 0) {
    //   std::cout << "ref #: " << i << std::endl;
    //   std::cout << " - n_bins: " << n_bins << std::endl;
    // }
    for (std::uint32_t j=0; j<n_bins; j++) {
      infile.read(reinterpret_cast<char *>(&bin_idx), sizeof(bin_idx));
      infile.read(reinterpret_cast<char *>(&loffset), sizeof(loffset));
      infile.read(reinterpret_cast<char *>(&n_chunks), sizeof(n_chunks));
      // std::cout << "bin_idx: " << bin_idx << ", virtual offset: " << loffset << ", n_chunks: " << n_chunks << std::endl;
      std::vector<Chunk> chunks;
      for (std::uint32_t k=0; k<n_chunks; k++) {
        infile.read(reinterpret_cast<char *>(&start), sizeof(start));
        infile.read(reinterpret_cast<char *>(&end), sizeof(end));
        chunks.push_back({ start, end });
      }
      bins.push_back({ bin_idx, loffset, chunks });
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
int IndexFile::reg2bins(std::int64_t beg, std::int64_t end) {
  // throw std::invalid_argument("not implemented yet");
  int max_bins = bin_limit();
  
  int l, t, n, s = min_shift + depth * 3;
  for (--end, l = n = t = 0; l <= depth; s -= 3, t += 1 << l * 3, ++l) {
    int b = t + (beg >> s), e = t + (end >> s), i;
    for (i = b; i <= e; ++i) {
      std::cout << "n: " << n+1 << ", i: " << i << std::endl;
      // bins[n++] = i;
      
      // I should use the vector of bins for a contig (chromosome).
      // Maybe an iterator of bins/chunks for the relevant region?
      // indices[contig][n++] = i;
    }
  }
  return n;
}

/* calculate maximum bin number -- valid bin numbers range within [0,bin_limit) */
int IndexFile::bin_limit() {
  return ((1 << (depth + 1) * 3) - 1) / 7;
}

std::uint64_t IndexFile::query(std::uint32_t contig_id, std::int64_t beg) {
  return 0;
}

}

// int main() {
//   bcf::IndexFile indexfile = bcf::IndexFile("/users/jmcrae/apps/pybcf/chr5.110000001-115000000.bcf.csi");
//   indexfile.reg2bins(111000000, 112000000);
//   std::cout << "bin: " << indexfile.reg2bin(111000000, 112000000) << std::endl;
//   return 0;
// }

// g++ -stdlib=libc++ -std=c++11 -lz index.cpp gzstream/gzstream.C
