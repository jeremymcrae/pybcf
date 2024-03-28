#ifndef BCF_BCF_H
#define BCF_BCF_H

#include <cstdint>
#include <string>
#include <vector>

#include "gzstream/gzstream.h"

#include "header.h"
#include "variant.h"

namespace bcf {

class BCF {
  std::ifstream handle;
  igzstream infile;
  Header header;
  std::string curr_line = "";
public:
  BCF(std::string path);
  Variant next();
};

}

#endif
