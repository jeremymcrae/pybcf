

#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "header.h"

namespace bcf {

/// @brief strip specified characters form the ends of a string
/// @param s string to trim on either ends
/// @param vals characters to trim from the ends
static std::string trim(const std::string &s, const std::string &vals) {
  size_t start = s.find_first_not_of(vals);
  if (start == std::string::npos) {
    return "";
  }
  size_t end = s.find_last_not_of(vals);
  return s.substr(start, (end + 1) - start);
}

/// @brief split line into key, remainder, unless it is a single field
/// @param line string to split apart
/// @return vector with either one or two entries
static std::vector<std::string> split_line(std::string line) {
  line = trim(line, "#");
  size_t delim_pos = line.find('=');
  std::vector<std::string> data;
  if (delim_pos != std::string::npos) {
    data.push_back(line.substr(0, delim_pos));
    data.push_back(line.substr(delim_pos + 1));
  } else {
    data.push_back(line);
  }
  return data;
}

/// @brief parse the comma-separated fields of an info/format/contig line
/// @param text 
/// @return 
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
        std::string key = item.substr(0, delim);
        std::string value = item.substr(delim + 1);
        data[key] = trim(value, "\"");
      }
      start = i + 1;
    }
  }
  return data;
}

Header::Header(const std::string & text) {
  filters[idx] = {"PASS", "All filters passed"};
  
  std::istringstream lines(text);
  std::string line;
  std::unordered_map<std::string, std::string> data;
  while (std::getline(lines, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == '\0')) {
      line.pop_back();
    }
    if (line.size() < 2 || line[0] != '#') {
      continue;
    }
    if (line[1] != '#') {
      if (line.substr(0, 6).find("CHROM") != std::string::npos) {
        if (samples.size() > 0) {
          throw std::invalid_argument("BCF has two header lines starting with #CHROM!");
        }
        // parse the sample IDs
        std::string item = "FORMAT\t";
        size_t i = line.find(item);
        if (i != std::string::npos) {
          line = line.substr(i + item.size());
          std::istringstream iss(line);
          while (std::getline(iss, item, '\t')) {
            samples.push_back(item);
          }
        }
      }
    } else {
      std::vector<std::string> parsed = split_line(line);
      if (parsed.size() == 2) {
        std::string id = parsed[0];
        std::string remainder = parsed[1];
        
        if (!(valid.count(id) > 0)) {
          continue;
        }
        
        data = parse_delimited(remainder);

        if (id == "contig") {
          if (data.count("IDX") > 0) {
            contig_idx = static_cast<std::uint32_t>(std::stoul(data["IDX"]));
          }
          contigs[contig_idx] = {data["ID"]};
          contig_idx += 1;
        } else {
          // use the IDX field for the map key, if IDX was included
          if (has_idx_tag && data.count("IDX") == 0) {
            throw std::invalid_argument("invalid BCF - missing IDX field in " + data["ID"]);
          }
          if (data.count("IDX") > 0) {
            idx = static_cast<std::uint32_t>(std::stoul(data["IDX"]));
            has_idx_tag = true;
          }
          if (id == "INFO") {
            info[idx] = {data["ID"], data["Number"], data["Type"], 
                         data["Description"]};
          } else if (id == "FORMAT") {
            format[idx] = {data["ID"], data["Number"], data["Type"], 
                           data["Description"]};
          } else if ((id == "FILTER") && (data["ID"] != "PASS")) {
            filters[idx] = {data["ID"], data["Description"]};
          }
          idx += 1;
        }
      }
    }
  }
  
}

/// @brief get contigs/chromosomes defined in BCF header
/// @return vector of chromosome/contig names
std::vector<std::string> Header::get_contigs() {
  std::vector<std::uint32_t> keys;
  keys.reserve(contigs.size());
  for (auto& it : contigs) {
      keys.push_back(it.first);
  }
  std::sort(keys.begin(), keys.end());
  
  std::vector<std::string> vals;
  vals.reserve(contigs.size());
  for (auto & k : keys) {
    vals.push_back(contigs[k].id);
  }
  return vals;
}

/// @brief get info fields used in the BCF header
/// @return vector of info keys (as strings) used in the BCF
std::vector<std::string> Header::get_info() {
  std::vector<std::uint32_t> keys;
  keys.reserve(info.size());
  for (auto& it : info) {
      keys.push_back(it.first);
  }
  std::sort(keys.begin(), keys.end());
  
  std::vector<std::string> vals;
  vals.reserve(info.size());
  for (auto & k : keys) {
    vals.push_back(info[k].id);
  }
  return vals;
}

/// @brief get filters used in the BCF
/// @return vector of unqiue filter values used in the BCF
std::vector<std::string> Header::get_filters() {
  std::vector<std::uint32_t> keys;
  keys.reserve(filters.size());
  for (auto& it : filters) {
      keys.push_back(it.first);
  }
  std::sort(keys.begin(), keys.end());
  
  std::vector<std::string> vals;
  vals.reserve(filters.size());
  for (auto & k : keys) {
    vals.push_back(filters[k].id);
  }
  return vals;
}

/// @brief get format fields used in the BCF header
/// @return vector of format keys (as strings) used in the BCF
std::vector<std::string> Header::get_formats() {
  std::vector<std::uint32_t> keys;
  keys.reserve(format.size());
  for (auto& it : format) {
      keys.push_back(it.first);
  }
  std::sort(keys.begin(), keys.end());
  
  std::vector<std::string> vals;
  vals.reserve(format.size());
  for (auto & k : keys) {
    vals.push_back(format[k].id);
  }
  return vals;
}

/// @brief find the numeric contig ID for a chromosome/contig name
/// @param contig chromosome/contig name to find
/// @return 
std::uint32_t Header::get_contig_id(std::string contig) {
  for (auto & x : contigs) {
    if (x.second.id == contig) {
      return x.first;
    }
  }

  throw std::invalid_argument(contig + " is not in this BCF");
}

}