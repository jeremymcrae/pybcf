

#include "header.h"

namespace bcf {

Header::Header(std::string & text) {
  lines.push_back(text);
}

}