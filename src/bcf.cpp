
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
  if (empty_region) {
    throw std::out_of_range("no variants in queried region");
  }
  try {
    Variant var(infile, header);
    if (!query_chrom.empty() && var.chrom != query_chrom) {
      throw std::out_of_range("variant not on required chrom");
    }

    while (static_cast<std::int64_t>(var.pos) + std::max<std::int32_t>(1, var.rlen) - 1 < query_start) {
      var = Variant(infile, header);
      if (!query_chrom.empty() && var.chrom != query_chrom) {
        throw std::out_of_range("variant not on required chrom");
      }
    }
    if (var.pos > query_end) {
      throw std::out_of_range("variant out of bounds");
    }
    return var;
  } catch (const std::out_of_range &) {
    empty_region = true;
    throw;
  }
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
  empty_region = false;
  try {
    Offsets offset = idxfile.query(contig_id, std::max<std::int64_t>(0, start - 1), end);
    infile.seek(offset);
  } catch (const std::out_of_range &) {
    empty_region = true;
  }
}

void BCF::close() {
  infile.close();
}

}