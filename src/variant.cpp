#include <climits>
#include <cstring>
#include <stdexcept>

#include "variant.h"

namespace bcf {

/// @brief start a Variant object
/// @param infile file to read variant data from
/// @param header for interpreting contig name, info and sample fields
Variant::Variant(GzReader & infile, Header & header) {
  std::uint32_t l_shared = 0;
  std::size_t n_read = infile.read(&l_shared, sizeof(std::uint32_t));
  if (n_read == 0) {
    throw std::out_of_range("reached end of file");
  }
  if (n_read != sizeof(std::uint32_t)) {
    throw std::invalid_argument("truncated BCF variant header");
  }
  if (l_shared < 24 || l_shared > UINT32_MAX - 4) {
    throw std::invalid_argument("invalid BCF shared record length");
  }

  std::uint32_t metadata_len = l_shared + 4;
  buf.resize(metadata_len);
  infile.read_exact(buf.data(), metadata_len, "truncated BCF variant record");

  std::uint32_t sampledata_len = 0;
  std::uint32_t idx = 0;
  std::memcpy(&sampledata_len, buf.data() + idx, 4);
  idx += 4;
  std::memcpy(&contig_idx, buf.data() + idx, 4);
  idx += 4;
  std::int32_t pos0 = 0;
  std::memcpy(&pos0, buf.data() + idx, 4);
  pos = pos0 + 1; // convert to 1-based coordinate
  idx += 4;
  std::memcpy(&rlen, buf.data() + idx, 4);
  idx += 4;

  std::uint32_t qual_bits = 0;
  std::memcpy(&qual_bits, buf.data() + idx, 4);
  if (!is_bcf_float_missing(qual_bits)) {
    std::memcpy(&qual, &qual_bits, 4);
  }
  idx += 4;

  if (contig_idx < 0 || header.contigs.count(static_cast<std::uint32_t>(contig_idx)) == 0) {
    throw std::invalid_argument("unknown contig index: " + std::to_string(contig_idx));
  }
  chrom = header.contigs.at(static_cast<std::uint32_t>(contig_idx)).id;

  std::uint32_t n_allele_info = 0;
  std::memcpy(&n_allele_info, buf.data() + idx, 4);
  idx += 4;
  n_alleles = n_allele_info >> 16;
  n_info = n_allele_info & 0xffff;

  if (n_alleles == 0) {
    throw std::invalid_argument(chrom + ":" + std::to_string(pos) + " lacks a ref allele");
  }

  std::uint32_t n_fmt_sample = 0;
  std::memcpy(&n_fmt_sample, buf.data() + idx, 4);
  idx += 4;
  n_sample = n_fmt_sample & 0xffffff;
  n_fmt = n_fmt_sample >> 24;

  // get variant ID
  Typed type_val(buf.data(), buf.size(), idx);
  varid = parse_string(buf.data(), buf.size(), idx, type_val.n_vals);

  // get ref allele. We previously raised an error if no ref allele exists
  type_val = Typed(buf.data(), buf.size(), idx);
  ref = parse_string(buf.data(), buf.size(), idx, type_val.n_vals);

  // get alt alleles
  alts.resize(n_alleles - 1);
  for (std::uint32_t i = 0; i < (n_alleles - 1); i++) {
    type_val = Typed(buf.data(), buf.size(), idx);
    alts[i] = parse_string(buf.data(), buf.size(), idx, type_val.n_vals);
  }

  // read the filter fields
  type_val = Typed(buf.data(), buf.size(), idx);
  filters.resize(type_val.n_vals);
  for (std::uint32_t i = 0; i < type_val.n_vals; i++) {
    std::int32_t filter_idx = parse_int(buf.data(), buf.size(), idx, type_val.type_size);
    if (filter_idx < 0 || header.filters.count(static_cast<std::uint32_t>(filter_idx)) == 0) {
      throw std::invalid_argument("unknown filter index: " + std::to_string(filter_idx));
    }
    filters[i] = header.filters.at(static_cast<std::uint32_t>(filter_idx)).id;
  }

  // prepare the info fields and format fields
  info = Info(buf.data() + idx, buf.size() - idx, &header, n_info);
  sample_data = SampleData(infile, header, sampledata_len, n_fmt, n_sample);
}

}
