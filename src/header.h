#ifndef BCF_HEADER_H
#define BCF_HEADER_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace bcf {

struct InfoField {
  std::string id;
  bool is_list;
};

class Header {
  std::unordered_map<std::string, std::uint32_t> contig_ids;
public:
  Header(const std::string &text);
  Header() = default;
  std::unordered_map<std::uint32_t, std::string> contigs;
  std::unordered_map<std::uint32_t, InfoField> info;
  std::unordered_map<std::uint32_t, std::string> filters;
  std::unordered_map<std::uint32_t, std::string> format;
  std::vector<std::string> samples;
  std::vector<std::string> get_contigs() const;
  std::vector<std::string> get_info() const;
  std::vector<std::string> get_filters() const;
  std::vector<std::string> get_formats() const;
  std::uint32_t get_contig_id(const std::string & contig) const;
};

}

#endif
