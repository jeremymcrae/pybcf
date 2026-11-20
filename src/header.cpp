#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "header.h"

namespace bcf {

/// @brief strip specified characters from the ends of a string
static std::string trim(const std::string &s, const std::string &vals) {
  size_t start = s.find_first_not_of(vals);
  if (start == std::string::npos) {
    return "";
  }
  size_t end = s.find_last_not_of(vals);
  return s.substr(start, (end + 1) - start);
}

/// @brief parse the comma-separated fields of an info/format/contig line
static std::unordered_map<std::string, std::string> parse_delimited(std::string text) {
  text = trim(text, "<>\r");
  std::unordered_map<std::string, std::string> data;

  bool in_quotes = false;
  size_t start = 0;
  for (size_t i = 0; i <= text.size(); i++) {
    if (i < text.size() && text[i] == '"' && (i == 0 || text[i - 1] != '\\')) {
      in_quotes = !in_quotes;
    } else if (i == text.size() || (text[i] == ',' && !in_quotes)) {
      std::string item = text.substr(start, i - start);
      size_t delim = item.find('=');
      if (delim != std::string::npos) {
        data[item.substr(0, delim)] = trim(item.substr(delim + 1), "\"");
      }
      start = i + 1;
    }
  }
  return data;
}

template <typename Map>
static std::vector<std::uint32_t> sorted_keys(const Map & m) {
  std::vector<std::uint32_t> keys;
  keys.reserve(m.size());
  for (const auto & kv : m) {
    keys.push_back(kv.first);
  }
  std::sort(keys.begin(), keys.end());
  return keys;
}

static std::vector<std::string> sorted_ids(const std::unordered_map<std::uint32_t, std::string> & m) {
  std::vector<std::string> vals;
  vals.reserve(m.size());
  for (std::uint32_t k : sorted_keys(m)) {
    vals.push_back(m.at(k));
  }
  return vals;
}

Header::Header(const std::string & text) {
  filters[0] = "PASS";
  std::uint32_t idx = 1;
  std::uint32_t contig_idx = 0;
  bool has_idx_tag = false;

  std::istringstream lines(text);
  std::string line;
  while (std::getline(lines, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == '\0')) {
      line.pop_back();
    }
    if (line.size() < 2 || line[0] != '#') {
      continue;
    }
    if (line[1] != '#') {
      if (line.rfind("#CHROM", 0) == 0) {
        if (!samples.empty()) {
          throw std::invalid_argument("BCF has two header lines starting with #CHROM!");
        }
        const std::string marker = "FORMAT\t";
        size_t pos = line.find(marker);
        if (pos != std::string::npos) {
          std::istringstream iss(line.substr(pos + marker.size()));
          std::string sample;
          while (std::getline(iss, sample, '\t')) {
            samples.push_back(sample);
          }
        }
      }
      continue;
    }

    std::string body = line.substr(2);
    size_t eq = body.find('=');
    if (eq == std::string::npos) {
      continue;
    }
    std::string kind = body.substr(0, eq);
    if (kind != "contig" && kind != "INFO" && kind != "FILTER" && kind != "FORMAT") {
      continue;
    }

    auto data = parse_delimited(body.substr(eq + 1));
    const std::string & id = data["ID"];
    auto idx_it = data.find("IDX");
    if (kind == "contig") {
      if (idx_it != data.end()) {
        contig_idx = static_cast<std::uint32_t>(std::stoul(idx_it->second));
      }
      contigs[contig_idx] = id;
      contig_ids[id] = contig_idx;
      contig_idx += 1;
    } else {
      if (has_idx_tag && idx_it == data.end()) {
        throw std::invalid_argument("invalid BCF - missing IDX field in " + id);
      }
      if (idx_it != data.end()) {
        idx = static_cast<std::uint32_t>(std::stoul(idx_it->second));
        has_idx_tag = true;
      } else if (kind == "FILTER" && id == "PASS") {
        // PASS is implicitly at index 0 when IDX tags are absent.
        continue;
      }

      if (kind == "INFO") {
        const std::string & num = data["Number"];
        bool is_list = !num.empty() && num != "1" && num != "0";
        info[idx] = {id, is_list};
      } else if (kind == "FORMAT") {
        format[idx] = id;
      } else if (kind == "FILTER" && id != "PASS") {
        filters[idx] = id;
      }
      idx += 1;
    }
  }
}

std::vector<std::string> Header::get_contigs() const { return sorted_ids(contigs); }
std::vector<std::string> Header::get_filters() const { return sorted_ids(filters); }
std::vector<std::string> Header::get_formats() const { return sorted_ids(format); }

std::vector<std::string> Header::get_info() const {
  std::vector<std::string> vals;
  vals.reserve(info.size());
  for (std::uint32_t k : sorted_keys(info)) {
    vals.push_back(info.at(k).id);
  }
  return vals;
}

std::uint32_t Header::get_contig_id(const std::string & contig) const {
  auto it = contig_ids.find(contig);
  if (it == contig_ids.end()) {
    throw std::invalid_argument(contig + " is not in this BCF");
  }
  return it->second;
}

}
