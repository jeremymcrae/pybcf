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

/// @brief parse `<Key=Value,...>` attributes from a `##INFO`/`##FORMAT`/`##FILTER`/`##contig` line
/// Tracks `"` quotes and `\` escapes so commas inside `Description="..."` do not split fields.
std::unordered_map<std::string, std::string> parse_delimited(std::string text) {
  text = trim(trim(text, " \t\r"), "<>");
  std::unordered_map<std::string, std::string> data;

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
        data[item.substr(0, delim)] = trim(item.substr(delim + 1), "\"");
      }
      start = i + 1;
    }
  }
  return data;
}

template <typename T>
void ensure_slot(std::vector<T> & vec, std::uint32_t idx) {
  if (idx >= vec.size()) {
    vec.resize(static_cast<std::size_t>(idx) + 1);
  }
}

} // namespace

Header::Header(const std::string & text) {
  // BCF headers have two separate integer dictionaries:
  //   1. Contigs: numbered from 0 in `##contig` order unless overridden by `IDX=...`.
  //   2. Strings: shared across `FILTER`, `INFO`, and `FORMAT`, with `PASS` implicitly
  //      assigned 0 and subsequent unique names numbered from 1 when `IDX=...` is absent.
  std::uint32_t idx = 1;
  std::uint32_t contig_idx = 0;
  std::unordered_map<std::string, std::uint32_t> dict_ids = {{"PASS", 0}};
  bool seen_idx = false;
  bool seen_no_idx = false;
  bool seen_chrom = false;
  bool seen_pass = false;

  // Scan the header text line by line in place (stripping trailing `\r` and BCF's trailing `\0`).
  size_t pos = 0;
  while (pos < text.size()) {
    size_t nl = text.find('\n', pos);
    size_t line_end = (nl == std::string::npos) ? text.size() : nl;
    size_t clean_end = line_end;
    while (clean_end > pos && (text[clean_end - 1] == '\r' || text[clean_end - 1] == '\0')) {
      clean_end--;
    }
    size_t line_len = clean_end - pos;
    size_t line_start = pos;
    pos = (nl == std::string::npos) ? text.size() : (nl + 1);

    if (line_len < 2 || text[line_start] != '#') {
      continue;
    }
    if (text[line_start + 1] != '#') {
      // Single-'#' `#CHROM` column header: sample IDs are the tab-separated columns after `FORMAT`.
      if ((line_len == 6 && text.compare(line_start, 6, "#CHROM") == 0) ||
          (line_len > 6 && text.compare(line_start, 7, "#CHROM\t") == 0)) {
        if (seen_chrom) {
          throw std::invalid_argument("BCF has two header lines starting with #CHROM!");
        }
        seen_chrom = true;
        const std::string marker = "FORMAT\t";
        size_t fmt_pos = text.find(marker, line_start);
        if (fmt_pos != std::string::npos && fmt_pos + marker.size() <= clean_end) {
          size_t s_pos = fmt_pos + marker.size();
          size_t n_tabs = static_cast<size_t>(
              std::count(text.begin() + s_pos, text.begin() + clean_end, '\t'));
          samples.reserve(n_tabs + 1);
          while (s_pos < clean_end) {
            size_t tab = text.find('\t', s_pos);
            if (tab == std::string::npos || tab > clean_end) {
              tab = clean_end;
            }
            samples.emplace_back(text, s_pos, tab - s_pos);
            s_pos = tab + 1;
          }
        }
      }
      continue;
    }

    // Double-'#' metadata line (e.g. `##INFO=<ID=DP,Number=1,Type=Integer,IDX=3>`).
    size_t body_start = line_start + 2;
    size_t eq = text.find('=', body_start);
    if (eq == std::string::npos || eq >= clean_end) {
      continue;
    }
    std::string kind(text, body_start, eq - body_start);
    if (kind != "contig" && kind != "INFO" && kind != "FILTER" && kind != "FORMAT") {
      continue;
    }

    auto data = parse_delimited(text.substr(eq + 1, clean_end - (eq + 1)));
    auto id_it = data.find("ID");
    if (id_it == data.end() || id_it->second.empty()) {
      throw std::invalid_argument("missing ID field in BCF header");
    }
    const std::string & id = id_it->second;
    auto idx_it = data.find("IDX");
    if (kind == "contig") {
      std::uint32_t entry_idx = 0;
      if (idx_it != data.end()) {
        entry_idx = parse_idx(idx_it->second, MAX_CONTIG_IDX);
        contig_idx = entry_idx + 1;
      } else {
        if (contig_ids.find(id) != contig_ids.end()) {
          continue;
        }
        entry_idx = contig_idx++;
      }
      ensure_slot(contigs, entry_idx);
      contigs[entry_idx] = id;
      contig_ids[id] = entry_idx;
    } else {
      if (idx_it == data.end() && kind == "FILTER" && id == "PASS") {
        // PASS is implicitly at index 0 when IDX tags are absent.
        continue;
      }
      // Reject headers that mix explicit `IDX=` lines with unindexed lines, since
      // implicit numbering could collide with explicit indices.
      bool has_idx = idx_it != data.end();
      if ((seen_idx && !has_idx) || (seen_no_idx && has_idx)) {
        throw std::invalid_argument("invalid BCF - inconsistent IDX fields in " + id);
      }
      std::uint32_t entry_idx = 0;
      if (has_idx) {
        entry_idx = parse_idx(idx_it->second, MAX_DICT_IDX);
        seen_idx = true;
      } else {
        seen_no_idx = true;
        // Reuse the same dictionary ID if a tag name (e.g. `DP`) appears in both INFO and FORMAT.
        auto existing = dict_ids.find(id);
        if (existing != dict_ids.end()) {
          entry_idx = existing->second;
        } else {
          entry_idx = idx++;
          dict_ids[id] = entry_idx;
        }
      }

      // Store the first definition seen for each `(kind, entry_idx)`.
      if (kind == "INFO") {
        const std::string & num = data["Number"];
        bool is_list = !num.empty() && num != "1" && num != "0";
        bool is_flag = data["Type"] == "Flag" || num == "0";
        ensure_slot(info, entry_idx);
        if (info[entry_idx].id.empty()) {
          info[entry_idx] = {id, is_list, is_flag};
        }
      } else if (kind == "FORMAT") {
        ensure_slot(format, entry_idx);
        if (format[entry_idx].id.empty()) {
          format[entry_idx] = {id, id == "GT"};
        }
      } else if (kind == "FILTER") {
        ensure_slot(filters, entry_idx);
        if (filters[entry_idx].empty()) {
          filters[entry_idx] = id;
        }
        if (id == "PASS") {
          seen_pass = true;
        }
      }
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
