

#include <sstream>

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
    data.push_back(line.substr(delim_pos + 1, line.size() - delim_pos));
  } else {
    data.push_back(line);
  }
  return data;
}

/// @brief parse the comma-separated fields of an info/format/contig line
/// @param text 
/// @return 
static std::unordered_map<std::string, std::string> parse_delimited(std::string text) {
  text = trim(text, "<>");
  std::unordered_map<std::string, std::string> data;

  size_t delim;
  std::string item, key, value;
  std::istringstream iss(text);

  while (std::getline(iss, item, ',')) {
    delim = item.find('=');
    key = item.substr(0, delim);
    value = item.substr(delim + 1, item.size() - delim);
    data[key] = trim(value, "\"");
  }
  return data;
}

Header::Header(std::string & text) {
  filters[idx] = {"PASS", "All filters passed"};
  
  std::istringstream lines(text);
  std::string line;
  std::unordered_map<std::string, std::string> data;
  bool is_valid = false;
  while (std::getline(lines, line)) {
    if (line[1] != '#') {
      continue;
    } else {
      std::vector<std::string> parsed = split_line(line);
      if (parsed.size() == 2) {
        std::string id = parsed[0];
        std::string remainder = parsed[1];
        is_valid = valid.count(id) > 0;
        if (is_valid) {
          data = parse_delimited(remainder);
        }
        
        if (id == "contig") {
          contigs[idx] = {data["ID"]};
        } else if (id == "INFO") {
          info[idx] = {data["ID"], data["Number"], data["Type"], 
                       data["Description"]};
        } else if (id == "FORMAT") {
          format[idx] = {data["ID"], data["Number"], data["Type"], 
                         data["Description"]};
        } else if (id == "FILTER") {
          filters[idx] = {data["ID"], data["Description"]};
        }
      idx += is_valid;
      }
    }
  }
  
}

}