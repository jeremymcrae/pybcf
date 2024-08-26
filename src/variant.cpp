

#include "variant.h"

namespace bcf {

/// @brief start a Varint object
/// @param infile file to read varint data from (as input gzstream)
/// @param header for interpreting contig name, info and sample fields
Variant::Variant(igzstream & infile,  Header & header) {

  std::uint32_t metadata_len=0;
  infile.read(reinterpret_cast<char *>(&metadata_len), sizeof(std::uint32_t));
  metadata_len += 4;
  
  if (infile.eof()) {
    throw std::out_of_range("reached end of file");
  }

  buf.resize(metadata_len);
  infile.read(reinterpret_cast<char *>(&buf[0]), metadata_len);

  std::uint32_t sampledata_len;
  std::uint32_t idx = 0;
  sampledata_len = *reinterpret_cast<std::uint32_t *>(&buf[idx]);
  idx += 4;
  contig_idx = *reinterpret_cast<std::int32_t *>(&buf[idx]);
  idx += 4;
  pos = *reinterpret_cast<std::int32_t *>(&buf[idx]) + 1; // convert to 1-based coordinate
  idx += 4;
  rlen = *reinterpret_cast<std::int32_t *>(&buf[idx]);
  idx += 4;
  
  if (*reinterpret_cast<std::uint32_t *>(&buf[idx]) != 0x7f800001) {
    qual = *reinterpret_cast<float *>(&buf[idx]);
  }
  idx += 4;

  chrom = header.contigs[contig_idx].id;
  
  std::uint32_t n_allele_info = *reinterpret_cast<std::uint32_t *>(&buf[idx]);;
  idx += 4;
  n_alleles = n_allele_info >> 16;
  n_info = n_allele_info & 0xffff;
  
  if (n_alleles == 0) {
    throw std::invalid_argument(chrom + ":" + std::to_string(pos) + " lacks a ref allele");
  }
  
  std::uint32_t n_fmt_sample = *reinterpret_cast<std::uint32_t *>(&buf[idx]);
  idx += 4;
  n_sample = n_fmt_sample & 0xffffff;
  n_fmt = n_fmt_sample >> 24;

  // get variant ID
  Typed type_val;
  type_val = {&buf[0], idx};
  varid = parse_string(&buf[0], idx, type_val.n_vals);
  
  // get ref allele. We previously raised an error if no ref allele exists
  type_val = {&buf[0], idx};
  ref = parse_string(&buf[0], idx, type_val.n_vals);

  // get alt alleles
  alts.resize(n_alleles - 1);
  for (std::uint32_t i = 0; i < (n_alleles - 1); i++) {
    type_val = {&buf[0], idx};
    alts[i] = parse_string(&buf[0], idx, type_val.n_vals);
  }

  // read the filter fields
  type_val = {&buf[0], idx};
  filters.resize(type_val.n_vals);
  for (std::uint32_t i = 0; i < type_val.n_vals; i++) {
    filters[i] = header.filters[parse_int(&buf[0], idx, type_val.type_size)].id;
  }
  
  // prepare the info fields and format fields
  info = Info(&buf[idx], &header, n_info);
  sample_data = SampleData(infile, header, sampledata_len, n_fmt, n_sample);
}



}