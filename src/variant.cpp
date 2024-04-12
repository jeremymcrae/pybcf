

#include "variant.h"

namespace bcf {

Variant::Variant(igzstream & infile,  Header & header) {

  infile.read(reinterpret_cast<char *>(&metadata_len), sizeof(std::uint32_t));
  infile.read(reinterpret_cast<char *>(&sampledata_len), sizeof(std::uint32_t));
  infile.read(reinterpret_cast<char *>(&contig_idx), sizeof(std::int32_t));
  infile.read(reinterpret_cast<char *>(&pos), sizeof(std::int32_t));
  pos += 1; // convert to 1-based coordinate
  infile.read(reinterpret_cast<char *>(&rlen), sizeof(std::int32_t));
  
  if (infile.eof()) {
    throw std::out_of_range("reached end of file");
  }
  
  std::uint32_t bytes;
  infile.read(reinterpret_cast<char *>(&bytes), sizeof(std::uint32_t));
  if (bytes != 0x7f800001) {
    qual = *reinterpret_cast<float *>(&bytes);
  }

  chrom = header.contigs[contig_idx].id;
  
  std::uint32_t n_allele_info;
  infile.read(reinterpret_cast<char *>(&n_allele_info), sizeof(std::uint32_t));
  n_alleles = n_allele_info >> 16;
  n_info = n_allele_info & 0xffff;
  
  alts.resize(n_alleles - 1);

  std::uint32_t n_fmt_sample;
  infile.read(reinterpret_cast<char *>(&n_fmt_sample), sizeof(std::uint32_t));
  n_sample = n_fmt_sample & 0xffffff;
  n_fmt = n_fmt_sample >> 24;

  // get variant ID
  std::uint8_t typing;
  infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
  Typed type_val = {typing, infile};
  varid.resize(type_val.n_vals);
  infile.read(reinterpret_cast<char *>(&varid[0]), type_val.n_vals);

  // get all alleles
  std::string allele;
  for (std::uint32_t i = 0; i < n_alleles; i++) {
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    allele.resize(type_val.n_vals);
    infile.read(reinterpret_cast<char *>(&allele[0]), type_val.n_vals);
    if (i == 0) {
      ref = allele;
    } else {
      alts[i-1] = allele;
    }
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
  info = Info(infile, header, n_info);
  sample_data = SampleData(infile, header, sampledata_len, n_fmt, n_sample);
}



}