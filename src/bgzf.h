#ifndef BCF_BGZF_H
#define BCF_BGZF_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "zlib.h"

namespace bcf {

// The most bytes a BGZF block can hold, both compressed and uncompressed.
const std::size_t BGZF_MAX_BLOCK_SIZE = 65536;

// Open a file for reading, or standard input if `path` is "-", returning its file
// descriptor, or -1 on failure.
int open_file(const std::string & path);
void close_file(int fd);

// Reads a BGZF file one block at a time. BGZF is gzip written as a series of small
// blocks, each of which can be decompressed on its own, so reading can start at any block.
//
// The file is read at offsets this class tracks itself, rather than at the file's own
// position. Copies of a reader in forked processes share that position, and so would
// otherwise move it under each other.
class BgzfFile {
  std::string path;                  // for error messages
  int fd = -1;
  std::uint64_t file_pos = 0;        // offset in the file of the next bytes to read
  std::vector<unsigned char> raw;    // compressed bytes read from the file
  std::size_t raw_len = 0;
  std::size_t raw_pos = 0;           // start of the next block in `raw`
  std::size_t compressed_size = 0;   // size of the next block in `raw`
  std::size_t uncompressed_size = 0; // size of the next block once decompressed
  z_stream strm{};                   // inflate state, reset for each block

  bool load(std::size_t count);
public:
  BgzfFile() = default;
  BgzfFile(const BgzfFile &) = delete;
  BgzfFile & operator=(const BgzfFile &) = delete;
  ~BgzfFile() { close(); }

  // Take over an open file, if it starts with a BGZF block. Otherwise return false,
  // leaving the file open and at its start. A pipe is never taken over.
  bool open(int open_fd, const std::string & filepath);
  void close();
  bool is_open() const { return fd >= 0; }

  // Find the next block, returning false at the end of the file.
  bool next_block();

  // Size of the block found by `next_block`, once decompressed.
  std::size_t block_size() const { return uncompressed_size; }

  // Decompress the block found by `next_block` into `dst`, which must hold `block_size()` bytes.
  void inflate_block(char * dst);

  // Carry on from the block that starts `offset` bytes into the file.
  void seek(std::uint64_t offset);
};

} // namespace bcf

#endif
