#ifndef BCF_INFO_H
#define BCF_INFO_H

#include <cstdint>
#include <string>
#include <vector>

#include "gzstream/gzstream.h"
#include "header.h"

namespace bcf {

class Info {
public:
  Info(igzstream & infile, Header & header, std::uint32_t n_info);
};

} // namespace

#endif
