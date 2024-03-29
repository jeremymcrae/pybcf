
#include <algorithm>
#include <iterator>
#include <stdexcept>

#include "bcf.h"

namespace bcf {


BCF::BCF(std::string path) {
  infile.open(path.c_str());
  if (infile.fail()) {
    throw std::invalid_argument("cannot open file at " + path);
  }
  
  // check the file header indicates this is a bcf file
  char magic[5];
  infile.read(&magic[0], 5);
  if (magic[0] != 'B' || magic[1] != 'C' || magic[2] != 'F' || magic[3] != 2 || magic[4] != 2) {
    throw std::invalid_argument("doesn't look like a BCF2.2 file");
  }
  
  std::uint32_t len;
  infile.read(reinterpret_cast<char *>(&len), sizeof(len));
  std::string text(len, ' ');
  infile.read(reinterpret_cast<char *>(&text[0]), len);
  header = Header(text);
}

Variant BCF::nextvar() {
  return Variant(infile, header);
}


}