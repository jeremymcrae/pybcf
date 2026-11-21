#ifndef BCF_GZREADER_H
#define BCF_GZREADER_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "bgzf.h"
#include "zlib.h"

namespace bcf {

// Buffered binary reader over a compressed BCF or CSI file, with BGZF virtual-offset seeking.
// BGZF files, the usual case, are decompressed one block at a time by `BgzfFile`.
// Plain gzip and uncompressed files fall back to zlib's streaming `gzread`, and can't seek.
class GzReader {
  std::string path;
  BgzfFile bgzf;             // reads the file, if it is BGZF
  gzFile gz = nullptr;       // reads any other file, through zlib's `gzread`
  std::vector<char> buffer;  // decompressed bytes not yet returned
  std::size_t buf_len = 0;
  std::size_t buf_pos = 0;

  std::size_t read_gz(char * dst, std::size_t count);
  bool fill();
public:
  GzReader() = default;
  explicit GzReader(const std::string & filepath) { open(filepath); }
  GzReader(const GzReader &) = delete;
  GzReader & operator=(const GzReader &) = delete;
  ~GzReader() { close(); }

  void open(const std::string & filepath);
  void close();
  bool is_open() const { return bgzf.is_open() || gz != nullptr; }

  // Read up to `count` bytes into `dst`; returns fewer than `count` only at EOF.
  std::size_t read(void * dst, std::size_t count);

  // Read `count` bytes into `dst`, throwing std::invalid_argument if EOF is hit early.
  void read_exact(void * dst, std::size_t count, const char * err_msg = "unexpected end of file");

  // Advance past `count` uncompressed bytes, throwing std::invalid_argument on early EOF.
  void skip(std::size_t count, const char * err_msg = "unexpected end of file");

  // Seek to a BGZF virtual offset. Its upper 48 bits are the offset in the file of a
  // compressed block, and its lower 16 bits are an offset in that block's uncompressed bytes.
  // Throws std::invalid_argument for files without BGZF compression.
  void seek(std::uint64_t virtual_offset);
};

} // namespace bcf

#endif
