#ifndef BCF_HEADER_H
#define BCF_HEADER_H

#include <string>
#include <vector>

namespace bcf {

class Header {
    std::vector<std::string> lines;
public:
  Header(std::string & text);
  Header() {}
};

}

#endif
