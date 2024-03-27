

#include "bcf.h"

namespace bcf {


BCF::BCF(std::string path) {
  std::ifstream tmp(path, std::ios::in | std::ios::binary);
  if (!tmp.is_open()) {
    throw std::invalid_argument("cannot open chain file at " + path);
  }
  igzstream infile(path.c_str());
}

Variant BCF::next() {
  
}


}