

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
  std::string header="";
  Variant next();
}

}