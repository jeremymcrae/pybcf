#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <utility>

#include "variant.h"

namespace bcf {

namespace {

// Fixed 32-byte binary header at the start of every BCF2 variant record.
// `l_shared` counts the bytes from `chrom_id` to the end of INFO (24 fixed bytes + variable data).
struct RecordHeader {
  std::uint32_t l_shared;
  std::uint32_t l_indiv;
  std::int32_t chrom_id;
  std::int32_t pos0;
  std::int32_t rlen;
  std::uint32_t qual_bits;
  std::uint32_t n_allele_info;
  std::uint32_t n_fmt_sample;
};
static_assert(sizeof(RecordHeader) == 32, "BCF2 fixed record header must be 32 bytes");

const std::uint32_t FIXED_SHARED_BYTES = 24;

} // namespace

bool Variant::read(GzReader & infile, const Header & header) {
  return read(infile, header, -1, 0, INT32_MAX);
}

/// @brief read the next variant record matching an optional region filter
bool Variant::read(GzReader & infile, const Header & header, std::int32_t query_contig_id,
                   std::int32_t query_start, std::int32_t query_end) {
  while (true) {
    RecordHeader hdr{};
    std::size_t n_read = infile.read(&hdr, sizeof(hdr));
    if (n_read == 0) {
      return false;
    }
    if (n_read != sizeof(hdr)) {
      throw std::invalid_argument("truncated BCF variant header");
    }
    if (hdr.l_shared < FIXED_SHARED_BYTES) {
      throw std::invalid_argument("invalid BCF shared record length");
    }
    if (hdr.chrom_id < 0 || static_cast<std::size_t>(hdr.chrom_id) >= header.contigs.size() ||
        header.contigs[hdr.chrom_id].empty()) {
      throw std::invalid_argument("unknown contig index: " + std::to_string(hdr.chrom_id));
    }
    if (hdr.pos0 < -1 || hdr.pos0 == INT32_MAX || hdr.rlen < 0) {
      throw std::invalid_argument("invalid BCF variant coordinates");
    }

    // Unpack allele/INFO counts (16 bits each) and FORMAT/sample counts (8 + 24 bits).
    std::uint32_t n_alleles = hdr.n_allele_info >> 16;
    std::uint32_t n_info = hdr.n_allele_info & 0xffff;
    if (n_alleles == 0) {
      throw std::invalid_argument(header.contigs[hdr.chrom_id] + ":" +
                                  std::to_string(hdr.pos0 + 1) + " lacks a ref allele");
    }

    std::uint32_t n_sample = hdr.n_fmt_sample & 0xffffff;
    std::uint32_t n_fmt = hdr.n_fmt_sample >> 24;
    if ((n_fmt > 0 || hdr.l_indiv > 0) && n_sample != header.samples.size()) {
      throw std::invalid_argument("BCF record sample count does not match header");
    }

    std::uint32_t shared_len = hdr.l_shared - FIXED_SHARED_BYTES;
    std::int32_t rec_pos = hdr.pos0 + 1;

    // During region queries, stop once past the target contig or end position, and
    // skip upstream records in the same BGZF block that end before `query_start`.
    if (query_contig_id >= 0 && hdr.chrom_id != query_contig_id) {
      return false;
    }
    if (rec_pos > query_end) {
      return false;
    }
    if (static_cast<std::int64_t>(rec_pos) + std::max<std::int32_t>(1, hdr.rlen) - 1 < query_start) {
      infile.skip(static_cast<std::size_t>(shared_len) + static_cast<std::size_t>(hdr.l_indiv),
                  "truncated BCF variant record");
      continue;
    }

    // Read the variable-length shared block (ID, REF, ALT, FILTER, and raw INFO bytes).
    std::unique_ptr<char[]> shared_buf;
    if (shared_len > 0) {
      shared_buf.reset(new char[shared_len]);
      infile.read_exact(shared_buf.get(), shared_len, "truncated BCF variant record");
    }
    const char * buf = shared_buf.get();

    chrom_id = static_cast<std::uint32_t>(hdr.chrom_id);
    pos = rec_pos;
    rlen = hdr.rlen;

    qual_bits = is_bcf_float_missing(hdr.qual_bits) ? BCF_FLOAT_MISSING : hdr.qual_bits;

    std::uint32_t idx = 0;
    varid = parse_typed_string(buf, shared_len, idx);
    ref = parse_typed_string(buf, shared_len, idx);

    alts.resize(n_alleles - 1);
    for (std::uint32_t i = 0; i < (n_alleles - 1); i++) {
      alts[i] = parse_typed_string(buf, shared_len, idx);
    }

    Typed type_val(buf, shared_len, idx);
    if (type_val.n_vals > 0 && !is_int_type(type_val.type)) {
      throw std::invalid_argument("expected integer BCF filter indices");
    }
    filter_ids.clear();
    filter_ids.reserve(type_val.n_vals);
    for (std::uint32_t i = 0; i < type_val.n_vals; i++) {
      std::int32_t filter_idx = parse_int(buf, shared_len, idx, type_val.type_size);
      if (filter_idx == INT32_MIN) {
        continue;
      }
      if (filter_idx < 0 || static_cast<std::size_t>(filter_idx) >= header.filters.size() ||
          header.filters[filter_idx].empty()) {
        throw std::invalid_argument("unknown filter index: " + std::to_string(filter_idx));
      }
      filter_ids.push_back(static_cast<std::uint32_t>(filter_idx));
    }

    // Hand ownership of `shared_buf` to Info so INFO values can be decoded lazily from `idx`.
    info = Info(std::move(shared_buf), shared_len, idx, &header, n_info);
    sample_data = SampleData(infile, header, hdr.l_indiv, n_fmt, n_sample);
    return true;
  }
}

}
