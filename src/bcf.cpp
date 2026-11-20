
#include <algorithm>
#include <iterator>
#include <stdexcept>

#include "bcf.h"

namespace bcf {


/// @brief open BCF file
/// @param path path to bcf file
/// @param index_path path to index file (if exists)
BCF::BCF(std::string path, std::string index_path) {
  infile.open(path);
  if (!infile.is_open()) {
    throw std::invalid_argument("cannot open file at " + path);
  }
  
  // check the file header indicates this is a bcf file
  char magic[5];
  infile.read_exact(magic, 5, "doesn't look like a BCF2.2 file");
  if (magic[0] != 'B' || magic[1] != 'C' || magic[2] != 'F' || magic[3] != 2 || magic[4] != 2) {
    throw std::invalid_argument("doesn't look like a BCF2.2 file");
  }
  
  std::uint32_t len = 0;
  infile.read_exact(&len, sizeof(len), "truncated BCF header length");
  std::string text(len, '\0');
  if (len > 0) {
    infile.read_exact(&text[0], len, "truncated BCF header text");
  }
  header = Header(text);
  
  bool explicit_index = !index_path.empty();
  if (!explicit_index) {
    index_path = path + ".csi";
  }
  idxfile.open(index_path, explicit_index);
}

// extract the next variant from the BCF, returning false at EOF or end of region
bool BCF::nextvar(Variant & var) {
  if (empty_region) {
    return false;
  }
  while (var.read(infile, header)) {
    if (!query_chrom.empty() && var.chrom != query_chrom) {
      empty_region = true;
      return false;
    }
    if (static_cast<std::int64_t>(var.pos) + std::max<std::int32_t>(1, var.rlen) - 1 < query_start) {
      continue;
    }
    if (var.pos > query_end) {
      empty_region = true;
      return false;
    }
    return true;
  }
  empty_region = true;
  return false;
}

/// @brief define genome region to extract variants from
/// @param chrom chromosome/contig name
/// @param start start position (in base pairs)
/// @param end end position (in base-pairs)
void BCF::set_region(std::string chrom, std::int32_t start, std::int32_t end) {
  if (!idxfile.has_index) {
    throw std::invalid_argument("cannot fetch without an index file");
  }
  if (end < start) {
    throw std::invalid_argument("start is after end: " + std::to_string(start) + " > " + std::to_string(end));
  }
  
  std::uint32_t contig_id = header.get_contig_id(chrom);
  
  query_chrom = chrom;
  query_start = start;
  query_end = end;

  Offsets offset;
  if (idxfile.query(contig_id, std::max<std::int64_t>(0, start - 1), end, offset)) {
    infile.seek(offset);
    empty_region = false;
  } else {
    empty_region = true;
  }
}

void BCF::close() {
  infile.close();
}

}