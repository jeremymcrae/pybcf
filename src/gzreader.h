#ifndef BCF_GZREADER_H
#define BCF_GZREADER_H

#include <cstddef>
#include <string>
#include <vector>
#include "zlib.h"

#include "index.h"

namespace bcf {

// Buffered binary reader over a zlib gzFile with BGZF virtual-offset seeking.
class GzReader {
  gzFile file = nullptr;
  std::string path;
  std::vector<char> buffer;
  std::size_t buf_len = 0;
  std::size_t buf_pos = 0;

  void check_error(int status);
public:
  GzReader() = default;
  explicit GzReader(const std::string & path) { open(path); }
  GzReader(const GzReader &) = delete;
  GzReader & operator=(const GzReader &) = delete;
  ~GzReader() { close(); }

  void open(const std::string & path);
  void close();
  bool is_open() const { return file != nullptr; }

  // Read up to `count` bytes into `dst`; returns fewer than `count` only at EOF.
  std::size_t read(void * dst, std::size_t count);

  // Read `count` bytes into `dst`, throwing std::invalid_argument if EOF is hit early.
  void read_exact(void * dst, std::size_t count, const char * err_msg = "unexpected end of file");

  // Seek to a BGZF virtual offset (compressed block offset + uncompressed intra-block offset).
  void seek(Offsets offset);
};

} // namespace bcf

#endif
