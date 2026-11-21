#ifndef BCF_VARIANT_H
#define BCF_VARIANT_H

#include <cstdint>
#include <string>
#include <vector>

#include "gzreader.h"
#include "header.h"
#include "info.h"
#include "sample_data.h"
#include "types.h"

namespace bcf {

// The region of a `fetch()` query, in 1-based positions that include both ends. The
// default covers every record in the file.
struct Region {
  std::int32_t contig_id = -1;  // index of the contig in the header, or -1 for any contig
  std::int32_t start = 0;
  std::int32_t end = INT32_MAX;
};

// Decoded site-level fields of a single BCF2 variant record, plus lazy `Info`
// and `SampleData` containers holding the raw INFO and per-sample byte blocks.
class Variant {
public:
  bool read(GzReader & infile, const Header & header, const Region & region);

  std::uint32_t chrom_id = 0;                 // contig dictionary index (`header.contigs[chrom_id]`)
  std::int32_t pos = 0;                       // 1-based VCF POS (0 for telomeric records)
  std::string ref;                            // reference allele sequence
  std::vector<std::string> alts;              // alternate alleles (excludes REF)
  std::uint32_t qual_bits = BCF_FLOAT_MISSING; // raw IEEE-754 float bits for QUAL
  std::string varid;                          // variant ID string (empty if missing `.`)
  std::vector<std::uint32_t> filter_ids;      // filter dictionary indices (`header.filters[...]`)
  Info info;
  SampleData sample_data;
};

}

#endif
