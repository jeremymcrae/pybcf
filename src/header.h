#ifndef BCF_HEADER_H
#define BCF_HEADER_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bcf {

// Metadata for a `##INFO` definition needed to interpret variant payloads.
struct InfoField {
  std::string id;
  bool is_list;  // true when `Number` is not 0 or 1 (e.g. "2", "A", "R", "G", ".")
  bool is_flag;  // true when `Type=Flag` or `Number=0`
};

// Metadata for a `##FORMAT` definition needed during sample decoding.
struct FormatField {
  std::string id;
  bool is_gt;    // true when `ID=GT` (alleles and phasing are bit-packed)
};

// Parsed VCF text header stored at the start of a BCF file.
// `contigs`, `info`, `filters`, and `format` are indexed directly by integer BCF `IDX`
// for O(1) lookup while reading variant records; unused `IDX` slots have an empty `id`.
class Header {
  std::unordered_map<std::string, std::uint32_t> contig_ids;
public:
  Header(const std::string & text);
  Header() = default;
  std::vector<std::string> contigs;
  std::vector<InfoField> info;
  std::vector<std::string> filters;
  std::vector<FormatField> format;
  std::vector<std::string> samples;

  // Compact lists of defined names (skipping empty IDX gaps) for Python `header.*` properties.
  std::vector<std::string> get_contigs() const;
  std::vector<std::string> get_info() const;
  std::vector<std::string> get_filters() const;
  std::vector<std::string> get_formats() const;

  // IDX-aligned tables (preserving empty slots) for O(1) ID -> name lookup in Cython.
  std::vector<std::string> get_info_by_idx() const;
  std::vector<std::string> get_formats_by_idx() const;
  std::uint32_t get_contig_id(const std::string & contig) const;
};

}

#endif
