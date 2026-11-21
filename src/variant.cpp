#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <utility>

#include "variant.h"

namespace bcf {

namespace {

// Fixed 32-byte binary header at the start of every BCF2 variant record.
// Every record is split into two blocks:
//   1. "shared" site block (`l_shared` bytes): 24 fixed bytes (`chrom_id..n_fmt_sample`)
//      followed by variable-length ID, REF, ALT, FILTER, and INFO fields.
//   2. "indiv" genotype block (`l_indiv` bytes): per-sample FORMAT values.
struct RecordHeader {
  std::uint32_t l_shared;
  std::uint32_t l_indiv;
  std::int32_t chrom_id;
  std::int32_t pos0;          // 0-based coordinate (-1 corresponds to 1-based VCF POS = 0)
  std::int32_t rlen;
  std::uint32_t qual_bits;
  std::uint32_t n_allele_info; // upper 16 bits = n_allele, lower 16 bits = n_info
  std::uint32_t n_fmt_sample;  // upper 8 bits = n_fmt, lower 24 bits = n_sample
};
static_assert(sizeof(RecordHeader) == 32, "BCF2 fixed record header must be 32 bytes");

const std::uint32_t FIXED_SHARED_BYTES = 24;

} // namespace

bool Variant::read(GzReader & infile, const Header & header) {
  return read(infile, header, -1, 0, INT32_MAX);
}

/// @brief read the next variant record matching an optional 1-based inclusive region filter
bool Variant::read(GzReader & infile, const Header & header, std::int32_t query_contig_id,
                   std::int32_t query_start, std::int32_t query_end) {
  while (true) {
    // 1. Read and validate the fixed 32-byte record header.
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

    // 2. Unpack bit-packed allele/INFO counts (16 + 16 bits) and FORMAT/sample counts (8 + 24 bits).
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

    // 3. During indexed `fetch()` queries (`query_contig_id >= 0`), check coordinates before
    // allocating or decoding the record payload:
    // - Because records are sorted by (contig, pos), moving onto another contig or past
    //   `query_end` means no later record can overlap the region.
    // - A CSI seek lands at the start of a BGZF block that may begin upstream of `query_start`;
    //   skip records whose reference span `[rec_pos, rec_pos + rlen - 1]` ends before `query_start`.
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

    // 4. Read the variable-length shared block (ID, REF, ALT, FILTER, and raw INFO bytes).
    std::unique_ptr<char[]> shared_buf;
    if (shared_len > 0) {
      shared_buf.reset(new char[shared_len]);
      infile.read_exact(shared_buf.get(), shared_len, "truncated BCF variant record");
    }
    const char * buf = shared_buf.get();

    chrom_id = static_cast<std::uint32_t>(hdr.chrom_id);
    pos = rec_pos;
    rlen = hdr.rlen;

    // Normalize any reserved BCF float sentinel in QUAL to `BCF_FLOAT_MISSING` while keeping NaN intact.
    qual_bits = is_bcf_float_missing(hdr.qual_bits) ? BCF_FLOAT_MISSING : hdr.qual_bits;

    // 5. Decode ID, REF, and the `n_alleles - 1` ALT strings sequentially from `shared_buf`.
    std::uint32_t idx = 0;
    varid = parse_typed_string(buf, shared_len, idx);
    ref = parse_typed_string(buf, shared_len, idx);

    alts.resize(n_alleles - 1);
    for (std::uint32_t i = 0; i < (n_alleles - 1); i++) {
      alts[i] = parse_typed_string(buf, shared_len, idx);
    }

    // 6. Decode the FILTER vector of header dictionary indices (skipping any MISSING entries).
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

    // 7. Transfer `shared_buf` into `Info` (which starts at `idx`) and read the raw sample block.
    info = Info(std::move(shared_buf), shared_len, idx, &header, n_info);
    sample_data = SampleData(infile, header, hdr.l_indiv, n_fmt, n_sample);
    return true;
  }
}

}
