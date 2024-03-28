

#include <stdexcept>

#include "bcf.h"

namespace bcf {


BCF::BCF(std::string path) {
  std::ifstream tmp(path, std::ios::in | std::ios::binary);
  if (!tmp.is_open()) {
    throw std::invalid_argument("cannot open file at " + path);
  }
  igzstream infile(path.c_str());
  
  // check the file header indicates this is a bcf file
  char magic[5];
  infile.read(&magic[0], 5);
  if (magic[0] != 'B' || magic[1] != 'C' || magic[2] != 'F' || magic[3] != '\2' || magic[4] != '\1') {
    throw std::invalid_argument("");
  }
  
  std::uint32_t header_len;
  handle.read(reinterpret_cast<char *>(&header_len), sizeof(header_len));
  
  std::string header_text(header_len, ' ');
  infile.read(reinterpret_cast<char *>(&header_text), header_len);
  header = Header(header_text);
}

Variant BCF::next() {
  return Variant(infile, header);
}


}