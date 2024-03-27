#ifndef BCF_BCF_H
#define BCF_BCF_H

#include <vector>
#include <string>

#include "gzstream/gzstream.h"

#include "variant.h"

namespace bcf {

class BCF {
  std::ifstream handle;
  igzstream infile;
public:
  BCF(std::string path);
  Variant next();
  std::string header="";
};

}

#endif
