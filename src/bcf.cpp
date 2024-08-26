
#include <algorithm>
#include <iterator>
#include <stdexcept>

#include "bcf.h"

namespace bcf {


/// @brief open BCF file
/// @param path path to bcf file
/// @param index_path path to index file (if exists)
BCF::BCF(std::string path, std::string index_path) {
  infile.open(path.c_str());
  if (infile.fail()) {
    throw std::invalid_argument("cannot open file at " + path);
  }
  
  // check the file header indicates this is a bcf file
  char magic[5];
  infile.read(&magic[0], 5);
  if (magic[0] != 'B' || magic[1] != 'C' || magic[2] != 'F' || magic[3] != 2 || magic[4] != 2) {
    throw std::invalid_argument("doesn't look like a BCF2.2 file");
  }
  
  std::uint32_t len;
  infile.read(reinterpret_cast<char *>(&len), sizeof(len));
  std::string text(len, ' ');
  infile.read(reinterpret_cast<char *>(&text[0]), len);
  header = Header(text);
  
  // try opening the index file
  if (index_path.size() == 0) {
    index_path = path + ".csi";
  }
  try {
    idxfile = IndexFile(index_path);
  }  catch (const std::invalid_argument& e) {}
  
}

// extract the next variant from the BCF
Variant BCF::nextvar() {
  Variant var = Variant(infile, header);
  if ((query_chrom.size() != 0) && (var.chrom != query_chrom)) {
    throw std::out_of_range("variant not on required chrom");
  }
  
  while (var.pos < query_start) {
    var = Variant(infile, header);
    if ((query_chrom.size() != 0) && (var.chrom != query_chrom)) {
      throw std::out_of_range("variant not on required chrom");
    }
  }
  if (var.pos > query_end) {
    throw std::out_of_range("variant out of bounds");
  }
  return var;
}

/// @brief define genome region to extract variants from
/// @param chrom chromosome/contig name
/// @param start start position (in base pairs)
/// @param end end position (in base-pairs)
void BCF::set_region(std::string chrom, std::int32_t start, std::int32_t end) {
  if (!idxfile.has_index) {
    throw std::invalid_argument("cannot fetch without an index file");
  }
  
  std::uint32_t contig_id = header.get_contig_id(chrom);
  
  Offsets offset = idxfile.query(contig_id, start, end);
  infile.seek(offset);
  
  query_chrom = chrom;
  query_start = start;
  query_end = end;
}

}