#include <algorithm>
#include <stdexcept>

#include "header.h"

namespace bcf {

namespace {

// Upper bounds on `IDX=` values, so a malformed header cannot make the lookup tables huge.
// Draft assemblies can have millions of contigs; FILTER/INFO/FORMAT names are far fewer.
const std::uint32_t MAX_CONTIG_IDX = 1u << 22;
const std::uint32_t MAX_DICT_IDX = 1u << 16;

/// @brief strip specified characters from the ends of a string
std::string trim(const std::string & s, const std::string & vals) {
  size_t start = s.find_first_not_of(vals);
  if (start == std::string::npos) {
    return "";
  }
  size_t end = s.find_last_not_of(vals);
  return s.substr(start, (end + 1) - start);
}

/// @brief parse the `IDX` attribute of a header line, which must be a number up to `max_idx`
std::uint32_t parse_idx(const std::string & s, std::uint32_t max_idx) {
  if (s.empty()) {
    throw std::invalid_argument("empty IDX field in BCF header");
  }
  for (unsigned char ch : s) {
    if (ch < '0' || ch > '9') {
      throw std::invalid_argument("invalid IDX field in BCF header: " + s);
    }
  }
  // Numbers of up to 9 digits always fit in 32 bits; longer ones are over any limit.
  if (s.size() > 9 || std::stoul(s) > max_idx) {
    throw std::invalid_argument("IDX field out of range in BCF header: " + s);
  }
  return static_cast<std::uint32_t>(std::stoul(s));
}

// The `Key=Value` attributes of one header line, such as `ID`, `Number`, `Type` and `IDX`.
typedef std::unordered_map<std::string, std::string> Attributes;

/// @brief parse `<Key=Value,...>` attributes from a `##INFO`/`##FORMAT`/`##FILTER`/`##contig` line
/// Tracks `"` quotes and `\` escapes so commas inside `Description="..."` do not split fields.
Attributes parse_attributes(std::string text) {
  text = trim(trim(text, " \t\r"), "<>");
  Attributes attrs;

  bool in_quotes = false;
  size_t start = 0;
  for (size_t i = 0; i <= text.size(); i++) {
    if (in_quotes && i + 1 < text.size() && text[i] == '\\') {
      i++;
      continue;
    }
    if (i < text.size() && text[i] == '"') {
      in_quotes = !in_quotes;
    } else if (i == text.size() || (text[i] == ',' && !in_quotes)) {
      std::string item = text.substr(start, i - start);
      size_t delim = item.find('=');
      if (delim != std::string::npos) {
        attrs[item.substr(0, delim)] = trim(item.substr(delim + 1), "\"");
      }
      start = i + 1;
    }
  }
  return attrs;
}

/// @brief get the sample IDs from the `#CHROM` line: the tab-separated columns after `FORMAT`
std::vector<std::string> parse_samples(const std::string & line) {
  std::vector<std::string> samples;
  const std::string marker = "FORMAT\t";
  size_t pos = line.find(marker);
  if (pos == std::string::npos) {
    return samples;
  }
  pos += marker.size();
  samples.reserve(std::count(line.begin() + pos, line.end(), '\t') + 1);
  while (pos < line.size()) {
    size_t tab = std::min(line.find('\t', pos), line.size());
    samples.emplace_back(line, pos, tab - pos);
    pos = tab + 1;
  }
  return samples;
}

// Numbers the names of FILTER, INFO and FORMAT lines, which share one dictionary. Either
// every line gives its number as `IDX`, or none does. Without `IDX`, `PASS` is 0 and the
// other names count up from 1 in order of first appearance, so a name used by both an
// INFO and a FORMAT line (e.g. `DP`) has a single number.
class StringDictionary {
  std::unordered_map<std::string, std::uint32_t> ids = {{"PASS", 0}};
  bool seen_idx = false;
  bool seen_no_idx = false;
public:
  std::uint32_t index_of(const std::string & id, const Attributes & attrs) {
    auto idx = attrs.find("IDX");
    bool has_idx = idx != attrs.end();
    // Numbers counted here could collide with numbers given as `IDX`, so reject a mix.
    if (has_idx ? seen_no_idx : seen_idx) {
      throw std::invalid_argument("invalid BCF - inconsistent IDX fields in " + id);
    }
    if (has_idx) {
      seen_idx = true;
      return parse_idx(idx->second, MAX_DICT_IDX);
    }
    seen_no_idx = true;
    // A new name takes the next number; `emplace` leaves a name already present alone.
    return ids.emplace(id, static_cast<std::uint32_t>(ids.size())).first->second;
  }
};

template <typename T>
void ensure_slot(std::vector<T> & vec, std::uint32_t idx) {
  if (idx >= vec.size()) {
    vec.resize(static_cast<std::size_t>(idx) + 1);
  }
}

} // namespace

Header::Header(const std::string & text) {
  // Records refer to header lines by number, using two separate dictionaries:
  //   1. Contigs: numbered from 0 in `##contig` order, unless a line gives an `IDX`.
  //   2. Strings: the names of `FILTER`, `INFO` and `FORMAT` lines (see `StringDictionary`).
  std::uint32_t next_contig = 0;
  StringDictionary dictionary;
  bool seen_chrom = false;
  bool seen_pass = false;

  size_t pos = 0;
  while (pos < text.size()) {
    // Take the next line, minus any trailing `\r`, or the `\0` that ends the header text.
    size_t line_end = std::min(text.find('\n', pos), text.size());
    std::string line(text, pos, line_end - pos);
    pos = line_end + 1;
    while (!line.empty() && (line.back() == '\r' || line.back() == '\0')) {
      line.pop_back();
    }
    if (line.size() < 2 || line[0] != '#') {
      continue;
    }

    if (line[1] != '#') {
      // The `#CHROM` line names the columns, which end with the sample IDs.
      if (line == "#CHROM" || line.compare(0, 7, "#CHROM\t") == 0) {
        if (seen_chrom) {
          throw std::invalid_argument("BCF has two header lines starting with #CHROM!");
        }
        seen_chrom = true;
        samples = parse_samples(line);
      }
      continue;
    }

    // A `##` line, e.g. `##INFO=<ID=DP,Number=1,Type=Integer,IDX=3>`. Only the four kinds
    // of line that records refer to by number are needed.
    size_t eq = line.find('=', 2);
    if (eq == std::string::npos) {
      continue;
    }
    std::string kind = line.substr(2, eq - 2);
    if (kind != "contig" && kind != "INFO" && kind != "FILTER" && kind != "FORMAT") {
      continue;
    }
    Attributes attrs = parse_attributes(line.substr(eq + 1));
    std::string id = attrs["ID"];
    if (id.empty()) {
      throw std::invalid_argument("missing ID field in BCF header");
    }

    if (kind == "contig") {
      auto idx = attrs.find("IDX");
      if (idx == attrs.end() && contig_ids.count(id) > 0) {
        continue; // a repeat of an earlier contig line
      }
      std::uint32_t contig_idx =
          (idx == attrs.end()) ? next_contig : parse_idx(idx->second, MAX_CONTIG_IDX);
      next_contig = contig_idx + 1;
      ensure_slot(contigs, contig_idx);
      contigs[contig_idx] = id;
      contig_ids[id] = contig_idx;
      continue;
    }

    if (kind == "FILTER" && id == "PASS" && attrs.count("IDX") == 0) {
      continue; // `PASS` without an `IDX` is number 0, which is filled in below
    }
    std::uint32_t idx = dictionary.index_of(id, attrs);

    // Keep the first definition seen for each number.
    if (kind == "INFO") {
      const std::string & num = attrs["Number"];
      bool is_list = !num.empty() && num != "1" && num != "0";
      bool is_flag = attrs["Type"] == "Flag" || num == "0";
      ensure_slot(info, idx);
      if (info[idx].id.empty()) {
        info[idx] = {id, is_list, is_flag};
      }
    } else if (kind == "FORMAT") {
      ensure_slot(format, idx);
      if (format[idx].id.empty()) {
        format[idx] = {id, id == "GT"};
      }
    } else {
      ensure_slot(filters, idx);
      if (filters[idx].empty()) {
        filters[idx] = id;
      }
      seen_pass = seen_pass || id == "PASS";
    }
  }

  // Ensure `PASS` is always available at filter index 0 unless explicitly placed elsewhere.
  if (!seen_pass && (filters.empty() || filters[0].empty())) {
    ensure_slot(filters, 0);
    filters[0] = "PASS";
  }
}

std::uint32_t Header::get_contig_id(const std::string & contig) const {
  auto it = contig_ids.find(contig);
  if (it == contig_ids.end()) {
    throw std::invalid_argument(contig + " is not in this BCF");
  }
  return it->second;
}

}
