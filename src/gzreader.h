#ifndef BCF_GZREADER_H
#define BCF_GZREADER_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "zlib.h"

namespace bcf {

// Decoded BGZF virtual offset: `c_offset` is the byte offset of a compressed BGZF
// block in the BCF file, and `u_offset` is the byte offset inside that uncompressed block.
struct Offsets {
  std::uint64_t u_offset; // uncompressed offset (within bgzf chunk)
  std::uint64_t c_offset; // compressed offset (within overall file)
};

// Buffered binary reader over a zlib gzFile with BGZF virtual-offset seeking.
class GzReader {
  gzFile file = nullptr;
  std::string path;
  std::vector<char> buffer;
  std::size_t buf_len = 0;
  std::size_t buf_pos = 0;

  std::size_t read_raw(char * dst, std::size_t count);
  bool fill();
public:
  GzReader() = default;
  explicit GzReader(const std::string & filepath) { open(filepath); }
  GzReader(const GzReader &) = delete;
  GzReader & operator=(const GzReader &) = delete;
  ~GzReader() { close(); }

  void open(const std::string & filepath);
  void close();
  bool is_open() const { return file != nullptr; }

  // Read up to `count` bytes into `dst`; returns fewer than `count` only at EOF.
  std::size_t read(void * dst, std::size_t count);

  // Read `count` bytes into `dst`, throwing std::invalid_argument if EOF is hit early.
  void read_exact(void * dst, std::size_t count, const char * err_msg = "unexpected end of file");

  // Advance past `count` uncompressed bytes, throwing std::invalid_argument on early EOF.
  void skip(std::size_t count, const char * err_msg = "unexpected end of file");

  // Seek to a BGZF virtual offset (compressed block offset + uncompressed intra-block offset).
  void seek(Offsets offset);
};

} // namespace bcf

#endif
