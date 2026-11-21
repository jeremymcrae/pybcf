#include <algorithm>
#include <stdexcept>

#include "bcf.h"

namespace bcf {

/// @brief open a BGZF-compressed BCF2.2 file and load its header and optional CSI index
/// @param path path to bcf file
/// @param index_path explicit CSI index path, or empty to probe `<path>.csi`
BCF::BCF(const std::string & path, const std::string & index_path) {
  infile.open(path);
  if (!infile.is_open()) {
    throw std::invalid_argument("cannot open file at " + path);
  }

  // Verify the 5-byte BCF2.2 magic header ("BCF\x02\x02").
  char magic[5];
  infile.read_exact(magic, 5, "doesn't look like a BCF2.2 file");
  if (magic[0] != 'B' || magic[1] != 'C' || magic[2] != 'F' || magic[3] != 2 || magic[4] != 2) {
    throw std::invalid_argument("doesn't look like a BCF2.2 file");
  }

  // Read the uint32 length of the NUL-terminated ASCII VCF header text, then parse it.
  // The length is untrusted, so grow the text in 1 MB chunks rather than allocating it upfront.
  std::uint32_t len = 0;
  infile.read_exact(&len, sizeof(len), "truncated BCF header length");
  std::string text;
  while (text.size() < len) {
    std::size_t start = text.size();
    text.resize(start + std::min<std::size_t>(len - start, 1 << 20));
    infile.read_exact(&text[start], text.size() - start, "truncated BCF header text");
  }
  header = Header(text);

  // Load `<path>.csi` opportunistically if present, but fail immediately if a
  // user-supplied `index_path` is missing.
  bool explicit_index = !index_path.empty();
  std::string csi_path = explicit_index ? index_path : (path + ".csi");
  idxfile.open(csi_path, explicit_index);
}

// Read the next variant matching the current stream/region state, returning false at EOF or end of region.
bool BCF::nextvar(Variant & var) {
  if (exhausted) {
    return false;
  }
  if (var.read(infile, header, query_contig_id, query_start, query_end)) {
    return true;
  }
  exhausted = true;
  return false;
}

/// @brief seek to a 1-based inclusive genomic region `[start, end]` using the CSI index
/// @param chrom chromosome/contig name
/// @param start 1-based start position in base pairs (0 includes telomeric POS = 0 records)
/// @param end 1-based inclusive end position in base pairs
void BCF::set_region(const std::string & chrom, std::int32_t start, std::int32_t end) {
  if (!idxfile.has_index) {
    throw std::invalid_argument("cannot fetch without an index file");
  }
  if (start < 0) {
    throw std::invalid_argument("start position must be non-negative: " + std::to_string(start));
  }
  if (end < start) {
    throw std::invalid_argument("start is after end: " + std::to_string(start) + " > " + std::to_string(end));
  }

  std::uint32_t contig_id = header.get_contig_id(chrom);

  query_contig_id = static_cast<std::int32_t>(contig_id);
  query_start = start;
  query_end = end;

  // CSI bins use 0-based half-open `[beg, end)` coordinates, so convert 1-based `start`
  // to `max(0, start - 1)`.
  Offsets offset;
  if (idxfile.query(contig_id, std::max<std::int64_t>(0, start - 1), end, offset)) {
    infile.seek(offset);
    exhausted = false;
  } else {
    exhausted = true;
  }
}

void BCF::close() {
  infile.close();
}

} // namespace bcf
