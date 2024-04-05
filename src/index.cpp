
#include <fstream>
#include <stdexcept>

#include "gzstream/gzstream.h"

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
  
  std::uint32_t n_bins, bin_idx;
  std::uint64_t loffset, start, end;
  std::int32_t n_chunks;
  for (std::uint32_t i=0; i< n_ref; i++) {
    std::vector<Bin> bins;
    infile.read(reinterpret_cast<char *>(&n_bins), sizeof(n_bins));
    for (std::uint32_t j=0; j<n_bins; j++) {
      infile.read(reinterpret_cast<char *>(&bin_idx), sizeof(bin_idx));
      infile.read(reinterpret_cast<char *>(&loffset), sizeof(loffset));
      infile.read(reinterpret_cast<char *>(&n_chunks), sizeof(n_chunks));
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
}

/// @brief calculate the list of bins that may overlap with region [beg,end) (zero-based).
///
/// This code is from https://samtools.github.io/hts-specs/CSIv1.pdf, but adapted
/// for being inside a class.
///
/// @param contig index of contig to check
/// @param beg start position of region
/// @param end end position of region
/// @return currently integer, but this should be an iterator instead
int IndexFile::reg2bins(size_t contig, std::int64_t beg, std::int64_t end) {
  throw std::invalid_argument("not implemented yet");
  
  int l, t, n, s = min_shift + depth * 3;
  for (--end, l = n = t = 0; l <= depth; s -= 3, t += 1 << l * 3, ++l) {
    int b = t + (beg >> s), e = t + (end >> s), i;
    for (i = b; i <= e; ++i) {
      // bins[n++] = i;
      
      // I should use the vector of bins for a contig (chromosome).
      // Maybe an iterator of bins/chunks for the relevant region?
      // indices[contig][n++] = i;
    }
  }
  return n;
}

std::uint64_t IndexFile::query() {
  return 0;
}

}

// int main() {
//   bcf::IndexFile indexfile = bcf::IndexFile("/users/jmcrae/apps/pybcf/test/data/hapmap_3.3.hg38.bcf.csi");
//   return 0;
// }

// g++ -stdlib=libc++ -std=c++11 -lz index.cpp gzstream/gzstream.C
