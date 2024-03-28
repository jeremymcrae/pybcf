#ifndef BCF_HEADER_H
#define BCF_HEADER_H

#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace bcf {

struct Contig {
  std::string id;
};

struct Filter {
  std::string id;
  std::string description;
};

struct Format {
  std::string id;
  std::string number;
  std::string type;
  std::string description;
};

struct Info {
  std::string id;
  std::string number;
  std::string type;
  std::string description;
};

class Header {
  std::unordered_set<std::string> valid = {"contig", "INFO", "FILTER", "FORMAT"};
  std::uint32_t idx=0;
public:
  Header(std::string &text);
  Header() {}
  std::unordered_map<std::uint32_t, Contig> contigs;
  std::unordered_map<std::uint32_t, Info> info;
  std::unordered_map<std::uint32_t, Filter> filters;
  std::unordered_map<std::uint32_t, Format> format;
};

}

#endif
