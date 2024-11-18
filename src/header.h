#ifndef BCF_HEADER_H
#define BCF_HEADER_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace bcf {

struct ContigField {
  std::string id;
};

struct FilterField {
  std::string id;
  std::string description;
};

struct FormatField {
  std::string id;
  std::string number;
  std::string type;
  std::string description;
};

struct InfoField {
  std::string id;
  std::string number;
  std::string type;
  std::string description;
};

class Header {
  std::unordered_set<std::string> valid = {"contig", "INFO", "FILTER", "FORMAT"};
  std::uint32_t idx=0;
  std::uint32_t contig_idx=0;
  bool has_idx_tag=false;
public:
  Header(std::string &text);
  Header() {}
  std::unordered_map<std::uint32_t, ContigField> contigs;
  std::unordered_map<std::uint32_t, InfoField> info;
  std::unordered_map<std::uint32_t, FilterField> filters;
  std::unordered_map<std::uint32_t, FormatField> format;
  std::vector<std::string> samples={};
  std::vector<std::string> get_contigs();
  std::vector<std::string> get_info();
  std::vector<std::string> get_filters();
  std::vector<std::string> get_formats();
  std::uint32_t get_contig_id(std::string contig);
};

}

#endif
