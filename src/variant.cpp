

#include "variant.h"

namespace bcf {

static std::string parse_string(const char * buf, std::uint32_t & idx, std::uint32_t size) {
  std::string val;
  val.resize(size);
  std::memcpy(&val[0], &buf[idx], size);
  idx += size;
  return val;
}

static std::int32_t parse_int(char * buf, std::uint32_t & idx, std::uint8_t type_size) {
  std::int32_t val=0;
  if (type_size == 1) {
    val = *reinterpret_cast<std::int8_t *>(&buf[idx]) & 0x000000FF;
    if (val == 0x80) { val = 0x80000000; }  // handle missing data value
  } else if (type_size == 2) {
    val = *reinterpret_cast<std::int16_t *>(&buf[idx]) & 0x0000FFFF;
    if (val == 0x8000) { val = 0x80000000; }  // handle missing data value
  } else {
    val = *reinterpret_cast<std::int32_t *>(&buf[idx]);
  }
  idx += type_size;
  return val;
}

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
  info = Info(&buf[0], &header, idx, n_info);
  sample_data = SampleData(infile, header, sampledata_len, n_fmt, n_sample);
}



}