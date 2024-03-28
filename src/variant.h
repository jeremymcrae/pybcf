#ifndef BCF_VARIANT_H
#define BCF_VARIANT_H

#include "gzstream/gzstream.h"

namespace bcf {


class Variant {
  
public:
    Variant(igzstream & infile);
};

}

#endif
