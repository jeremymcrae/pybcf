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

// Buffered binary reader over a compressed BCF or CSI file, with BGZF virtual-offset seeking.
// BGZF files, the usual case, are decompressed one block at a time with zlib's `inflate`.
// Plain gzip and uncompressed files fall back to zlib's streaming `gzread`, and can't seek.
//
// A BGZF file is read at offsets this reader tracks itself, rather than at the file's own
// position. Copies of a reader in forked processes share that position, and so would
// otherwise move it under each other.
class GzReader {
  std::string path;
  int bgzf = -1;                     // BGZF only: the open file
  std::uint64_t file_pos = 0;        // BGZF only: offset in the file of the next bytes to read
  std::vector<unsigned char> raw;    // BGZF only: compressed bytes read from the file
  std::size_t raw_len = 0;
  std::size_t raw_pos = 0;           // BGZF only: start of the next block in `raw`
  std::size_t compressed_size = 0;   // BGZF only: size of the next block in `raw`
  std::size_t block_size = 0;        // BGZF only: uncompressed size of the next block
  z_stream strm{};                   // BGZF only: inflate state, reset for each block
  gzFile file = nullptr;             // other files: zlib's streaming reader
  std::vector<char> buffer;          // decompressed bytes not yet returned
  std::size_t buf_len = 0;
  std::size_t buf_pos = 0;

  bool load(std::size_t count);
  bool next_block();
  void inflate_block(char * dst);
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
  bool is_open() const { return bgzf >= 0 || file != nullptr; }

  // Read up to `count` bytes into `dst`; returns fewer than `count` only at EOF.
  std::size_t read(void * dst, std::size_t count);

  // Read `count` bytes into `dst`, throwing std::invalid_argument if EOF is hit early.
  void read_exact(void * dst, std::size_t count, const char * err_msg = "unexpected end of file");

  // Advance past `count` uncompressed bytes, throwing std::invalid_argument on early EOF.
  void skip(std::size_t count, const char * err_msg = "unexpected end of file");

  // Seek to a BGZF virtual offset (compressed block offset + uncompressed intra-block offset).
  // Throws std::invalid_argument for files without BGZF compression.
  void seek(Offsets offset);
};

} // namespace bcf

#endif
