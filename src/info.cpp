
#include <iostream>

#include "info.h"
#include "types.h"

namespace bcf {


Info::Info(igzstream & infile, Header & header, std::uint32_t n_info) {

  std::uint8_t typing;
  Typed type_val;

  // read the info fields. TODO - find out a way to skip this if not required
  std::uint32_t info_idx;
  std::string key;
  for (std::uint32_t i = 0; i < n_info; i++) {
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    infile.read(reinterpret_cast<char *>(&info_idx), type_val.type_size);
    key = header.info[info_idx].id;

    std::cout << "info key: " << key << std::endl;

    // now parse value
    infile.read(reinterpret_cast<char *>(&typing), sizeof(std::uint8_t));
    type_val = {typing, infile};
    if (type_val.n_vals == 1)
    {
    }
  }
}


}