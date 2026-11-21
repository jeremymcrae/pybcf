#include <algorithm>
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

  std::uint32_t n_allele() const { return n_allele_info >> 16; }
  std::uint32_t n_info() const { return n_allele_info & 0xffff; }
  std::uint32_t n_fmt() const { return n_fmt_sample >> 24; }
  std::uint32_t n_sample() const { return n_fmt_sample & 0xffffff; }
};
static_assert(sizeof(RecordHeader) == 32, "BCF2 fixed record header must be 32 bytes");

const std::uint32_t FIXED_SHARED_BYTES = 24;

// Read the fixed-size header of the next record, and check it against the file's header.
// Returns false at the end of the file.
bool read_record_header(GzReader & infile, const Header & header, RecordHeader & hdr) {
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
  if (hdr.n_allele() == 0) {
    throw std::invalid_argument(header.contigs[hdr.chrom_id] + ":" +
                                std::to_string(hdr.pos0 + 1) + " lacks a ref allele");
  }
  if ((hdr.n_fmt() > 0 || hdr.l_indiv > 0) && hdr.n_sample() != header.samples.size()) {
    throw std::invalid_argument("BCF record sample count does not match header");
  }
  return true;
}

// Records are sorted by contig, then position. So once a record is on another contig, or
// starts after the region, no later record can overlap the region.
bool is_past(const RecordHeader & hdr, const Region & region) {
  return (region.contig_id >= 0 && hdr.chrom_id != region.contig_id) || hdr.pos0 + 1 > region.end;
}

// A `fetch()` starts reading where the index says the first overlapping record could be,
// which can be upstream of the region. The records to skip there are the ones whose
// reference span, `[pos, pos + rlen - 1]`, ends before the region starts.
bool ends_before(const RecordHeader & hdr, const Region & region) {
  std::int64_t pos = static_cast<std::int64_t>(hdr.pos0) + 1;
  return pos + std::max<std::int32_t>(1, hdr.rlen) - 1 < region.start;
}

// Find the next record that overlaps `region`, leaving `infile` just after the record's
// fixed-size header. Records before the region are skipped without decoding them.
// Returns false if no record is left in the region.
bool find_record(GzReader & infile, const Header & header, const Region & region,
                 RecordHeader & hdr) {
  while (read_record_header(infile, header, hdr) && !is_past(hdr, region)) {
    if (!ends_before(hdr, region)) {
      return true;
    }
    infile.skip(static_cast<std::size_t>(hdr.l_shared - FIXED_SHARED_BYTES) + hdr.l_indiv,
                "truncated BCF variant record");
  }
  return false;
}

// Decode FILTER, a vector of indices into the header's filters. MISSING entries are skipped.
std::vector<std::uint32_t> parse_filters(Cursor & cursor, const Header & header) {
  TypeDescriptor desc(cursor);
  if (desc.n_vals > 0 && !is_int_type(desc.type)) {
    throw std::invalid_argument("expected integer BCF filter indices");
  }
  std::vector<std::uint32_t> filter_ids;
  for (std::uint32_t i = 0; i < desc.n_vals; i++) {
    std::int32_t filter_idx = cursor.read_int(desc.type_size);
    if (filter_idx == INT32_MIN) {
      continue;
    }
    if (filter_idx < 0 || static_cast<std::size_t>(filter_idx) >= header.filters.size() ||
        header.filters[filter_idx].empty()) {
      throw std::invalid_argument("unknown filter index: " + std::to_string(filter_idx));
    }
    filter_ids.push_back(static_cast<std::uint32_t>(filter_idx));
  }
  return filter_ids;
}

} // namespace

/// @brief read the next variant record that overlaps `region`, returning false if there are none
bool Variant::read(GzReader & infile, const Header & header, const Region & region) {
  RecordHeader hdr{};
  if (!find_record(infile, header, region, hdr)) {
    return false;
  }
  chrom_id = static_cast<std::uint32_t>(hdr.chrom_id);
  pos = hdr.pos0 + 1;
  // Normalize any reserved BCF float sentinel in QUAL to `BCF_FLOAT_MISSING` while keeping NaN intact.
  qual_bits = is_bcf_float_missing(hdr.qual_bits) ? BCF_FLOAT_MISSING : hdr.qual_bits;

  // Read the rest of the shared block: ID, REF, ALT, FILTER, and the raw INFO bytes.
  std::uint32_t shared_len = hdr.l_shared - FIXED_SHARED_BYTES;
  std::unique_ptr<char[]> shared_buf;
  if (shared_len > 0) {
    shared_buf.reset(new char[shared_len]);
    infile.read_exact(shared_buf.get(), shared_len, "truncated BCF variant record");
  }

  // Decode the fields in front of INFO, which follow one another in the buffer.
  Cursor cursor(shared_buf.get(), shared_len);
  varid = cursor.read_typed_string();
  ref = cursor.read_typed_string();
  alts.resize(hdr.n_allele() - 1);
  for (std::string & alt : alts) {
    alt = cursor.read_typed_string();
  }
  filter_ids = parse_filters(cursor, header);

  // `Info` takes over the buffer, where INFO starts at the cursor, and `SampleData` reads
  // the per-sample block. Both decode their fields only when asked to.
  info = Info(std::move(shared_buf), shared_len, cursor.offset(), &header, hdr.n_info());
  sample_data = SampleData(infile, header, hdr.l_indiv, hdr.n_fmt(), hdr.n_sample());
  return true;
}

}
