

#include "variant.h"

namespace bcf {

Variant::Variant(igzstream & infile,  Header & header) {

  char buf[32];
  infile.read(reinterpret_cast<char *>(&buf[0]), 32);
  
  if (infile.eof()) {
    throw std::out_of_range("reached end of file");
  }

  std::uint32_t metadata_len;
  std::uint32_t sampledata_len;
  
  metadata_len = *reinterpret_cast<std::uint32_t *>(&buf[0]);
  sampledata_len = *reinterpret_cast<std::uint32_t *>(&buf[4]);
  contig_idx = *reinterpret_cast<std::int32_t *>(&buf[8]);
  pos = *reinterpret_cast<std::int32_t *>(&buf[12]) + 1; // convert to 1-based coordinate
  rlen = *reinterpret_cast<std::int32_t *>(&buf[16]);
  
  if (*reinterpret_cast<std::uint32_t *>(&buf[20]) != 0x7f800001) {
    qual = *reinterpret_cast<float *>(&buf[20]);
  }

  chrom = header.contigs[contig_idx].id;
  
  std::uint32_t n_allele_info = *reinterpret_cast<std::uint32_t *>(&buf[24]);;
  n_alleles = n_allele_info >> 16;
  n_info = n_allele_info & 0xffff;
  
  if (n_alleles == 0) {
    throw std::invalid_argument(chrom + ":" + std::to_string(pos) + " lacks a ref allele");
  }
  
  std::uint32_t n_fmt_sample = *reinterpret_cast<std::uint32_t *>(&buf[28]);
  n_sample = n_fmt_sample & 0xffffff;
  n_fmt = n_fmt_sample >> 24;

  // get variant ID
  std::uint8_t typing;
  infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
  Typed type_val = {typing, infile};
  varid.resize(type_val.n_vals);
  infile.read(reinterpret_cast<char *>(&varid[0]), type_val.n_vals);

  // get ref allele. We previously raised an error if no ref allele exists
  infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
  type_val = {typing, infile};
  ref.resize(type_val.n_vals);
  infile.read(reinterpret_cast<char *>(&ref[0]), type_val.n_vals);

  // get alt alleles
  alts.resize(n_alleles - 1);
  std::string allele;
  for (std::uint32_t i = 0; i < (n_alleles - 1); i++) {
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    allele.resize(type_val.n_vals);
    infile.read(reinterpret_cast<char *>(&allele[0]), type_val.n_vals);
    alts[i] = allele;
  }

  // read the filter fields
  infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
  type_val = {typing, infile};
  filters.resize(type_val.n_vals);
  std::uint32_t filter_idx;
  for (std::uint32_t i = 0; i < type_val.n_vals; i++) {
    infile.read(reinterpret_cast<char *>(&filter_idx), type_val.type_size);
    filters[i] = header.filters[filter_idx].id;
  }
  
  // read the info fields. TODO - find out a way to skip this if not required
  info = Info(infile, header, metadata_len, n_info);
  sample_data = SampleData(infile, header, sampledata_len, n_fmt, n_sample);
}



}