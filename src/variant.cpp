#include <cstring>
#include <stdexcept>

#include "variant.h"

namespace bcf {

namespace {

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

/// @brief read the next variant record from a BCF stream
/// @param infile file to read variant data from
/// @param header for interpreting contig name, info and sample fields
/// @return false at clean EOF, true when a variant was read
bool Variant::read(GzReader & infile, const Header & header) {
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

  std::uint32_t shared_len = hdr.l_shared - FIXED_SHARED_BYTES;
  std::vector<char> buf(shared_len);
  if (shared_len > 0) {
    infile.read_exact(buf.data(), shared_len, "truncated BCF variant record");
  }

  auto contig_it = (hdr.chrom_id >= 0) ? header.contigs.find(static_cast<std::uint32_t>(hdr.chrom_id))
                                       : header.contigs.end();
  if (contig_it == header.contigs.end()) {
    throw std::invalid_argument("unknown contig index: " + std::to_string(hdr.chrom_id));
  }
  if (hdr.pos0 < -1 || hdr.pos0 == INT32_MAX || hdr.rlen < 0) {
    throw std::invalid_argument("invalid BCF variant coordinates");
  }
  chrom = contig_it->second;
  pos = hdr.pos0 + 1; // convert 0-based BCF coordinate to 1-based VCF coordinate
  rlen = hdr.rlen;

  qual = std::nanf("");
  if (!is_bcf_float_missing(hdr.qual_bits)) {
    std::memcpy(&qual, &hdr.qual_bits, sizeof(qual));
  }

  std::uint32_t n_alleles = hdr.n_allele_info >> 16;
  std::uint32_t n_info = hdr.n_allele_info & 0xffff;
  if (n_alleles == 0) {
    throw std::invalid_argument(chrom + ":" + std::to_string(pos) + " lacks a ref allele");
  }

  std::uint32_t n_sample = hdr.n_fmt_sample & 0xffffff;
  std::uint32_t n_fmt = hdr.n_fmt_sample >> 24;
  if ((n_fmt > 0 || hdr.l_indiv > 0) && n_sample != header.samples.size()) {
    throw std::invalid_argument("BCF record sample count does not match header");
  }

  std::uint32_t idx = 0;
  varid = parse_typed_string(buf.data(), buf.size(), idx);
  ref = parse_typed_string(buf.data(), buf.size(), idx);

  alts.resize(n_alleles - 1);
  for (std::uint32_t i = 0; i < (n_alleles - 1); i++) {
    alts[i] = parse_typed_string(buf.data(), buf.size(), idx);
  }

  Typed type_val(buf.data(), buf.size(), idx);
  if (type_val.n_vals > 0 && !is_int_type(type_val.type)) {
    throw std::invalid_argument("expected integer BCF filter indices");
  }
  filters.resize(type_val.n_vals);
  for (std::uint32_t i = 0; i < type_val.n_vals; i++) {
    std::int32_t filter_idx = parse_int(buf.data(), buf.size(), idx, type_val.type_size);
    auto filt_it = (filter_idx >= 0) ? header.filters.find(static_cast<std::uint32_t>(filter_idx))
                                     : header.filters.end();
    if (filt_it == header.filters.end()) {
      throw std::invalid_argument("unknown filter index: " + std::to_string(filter_idx));
    }
    filters[i] = filt_it->second;
  }

  info = Info(buf.data() + idx, buf.size() - idx, &header, n_info);
  sample_data = SampleData(infile, header, hdr.l_indiv, n_fmt, n_sample);
  return true;
}

}
