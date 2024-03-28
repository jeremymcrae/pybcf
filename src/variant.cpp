

#include "variant.h"

namespace bcf {

enum Types {
  int8=1,
  int16=2,
  int32=3,
  float_=5,
  char_=7,
};

class Typed {
public:
  Types type;
  std::uint32_t size;
  Typed(std::uint8_t byte, igzstream & infile) {
    type = Types(byte << 4);
    size = byte >> 4;
    if (size == 15) {
      switch (type) {
        case int8:
          infile.read(reinterpret_cast<char *>(&size), 1);
          break;
        case int16:
          infile.read(reinterpret_cast<char *>(&size), 2);
          break;
        case int32:
          infile.read(reinterpret_cast<char *>(&size), 4);
          break;
        case float_:
          infile.read(reinterpret_cast<char *>(&size), 4);
          break;
        case char_:
          infile.read(reinterpret_cast<char *>(&size), 1);
          break;
      }
    }
  }
  std::int8_t type_size() {
    if (type == int8) {return 1;}
    if (type == int16) {return 2;}
    if (type == int32) {return 4;}
    if (type == float_) {return 4;}
    if (type == char_) {return 1;}
    return 1;
  }
};

Variant::Variant(igzstream & infile,  Header & _header) {
  header = _header;

  infile.read(reinterpret_cast<char *>(&metadata_len), sizeof(std::uint32_t));
  infile.read(reinterpret_cast<char *>(&sampledata_len), sizeof(std::uint32_t));
  infile.read(reinterpret_cast<char *>(&chrom), sizeof(std::int32_t));
  infile.read(reinterpret_cast<char *>(&pos), sizeof(std::int32_t));
  infile.read(reinterpret_cast<char *>(&rlen), sizeof(std::int32_t));
  infile.read(reinterpret_cast<char *>(&qual), sizeof(float));
  
  std::uint32_t n_allele_info;
  infile.read(reinterpret_cast<char *>(&n_allele_info), sizeof(std::uint32_t));
  n_alleles = n_allele_info << 16;
  n_info = n_allele_info >> 16;
  
  alts.resize(n_alleles - 1);

  std::uint32_t n_fmt_sample;
  infile.read(reinterpret_cast<char *>(&n_fmt_sample), sizeof(std::uint32_t));
  n_fmt = n_fmt_sample << 16;
  n_sample = n_fmt_sample >> 16;
  
  // get variant ID
  std::uint8_t typing;
  infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
  Typed type_val = {typing, infile};
  varid.resize(type_val.size);
  infile.read(reinterpret_cast<char *>(&varid), type_val.size);

  // get all alleles
  std::string allele;
  for (std::uint32_t i = 0; i < n_alleles; i++) {
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    allele.resize(type_val.size);
    infile.read(reinterpret_cast<char *>(&allele), type_val.size);
    if (i == 0) {
      ref = allele;
    } else {
      alts[i-1] = allele;
    }
  }

  // read the filter fields
  infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
  type_val = {typing, infile};
  std::uint32_t n_filters = type_val.size;
  filters.resize(n_filters);
  std::uint32_t filter_idx;
  for (std::uint32_t i = 0; i < n_filters; i++) {
    infile.read(reinterpret_cast<char *>(&filter_idx), type_val.type_size());
    filters[i] = header.filters[filter_idx].id;
  }
  
  // read the info fields. TODO - find out a way to skip this if not required
  std::uint32_t info_idx;
  std::string key;
  for (std::uint32_t i = 0; i < n_info; i++) {
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    infile.read(reinterpret_cast<char *>(&info_idx), type_val.type_size());
    key = header.info[info_idx].id;
    
    // now parse value
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
  }
}



}