#include "bgzf.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>

#include <fcntl.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace bcf {

namespace {

// Every BGZF block starts with an 18-byte gzip header and ends with its CRC32 and
// uncompressed size (4 bytes each).
const std::size_t HEADER_SIZE = 18;
const std::size_t FOOTER_SIZE = 8;

// Read up to `count` bytes starting `pos` bytes into the file, without moving the file's
// own position. Returns how many bytes were read, or -1 on failure, as for a pipe.
long long read_at(int fd, unsigned char * dst, std::size_t count, std::uint64_t pos) {
#ifdef _WIN32
  // Windows has no `pread`, so move to `pos`, read, then move back.
  __int64 start = _lseeki64(fd, 0, SEEK_CUR);
  if (start < 0 || _lseeki64(fd, static_cast<__int64>(pos), SEEK_SET) < 0) {
    return -1;
  }
  int got = _read(fd, dst, static_cast<unsigned int>(count));
  _lseeki64(fd, start, SEEK_SET);
  return got;
#else
  return pread(fd, dst, count, static_cast<off_t>(pos));
#endif
}

// A BGZF block is a gzip member whose header has a 'BC' extra field holding the block's size.
bool is_bgzf_header(const unsigned char * h) {
  return h[0] == 0x1f && h[1] == 0x8b && h[2] == 8 && (h[3] & 4) && h[10] == 6 && h[11] == 0 &&
         h[12] == 'B' && h[13] == 'C' && h[14] == 2 && h[15] == 0;
}

} // namespace

int open_file(const std::string & path) {
  // "-" means standard input. It is duplicated, so closing the file leaves stdin open.
#ifdef _WIN32
  if (path == "-") {
    int fd = _dup(0);
    if (fd >= 0) {
      _setmode(fd, _O_BINARY);
    }
    return fd;
  }
  return _open(path.c_str(), _O_RDONLY | _O_BINARY);
#else
  if (path == "-") {
    return ::dup(STDIN_FILENO);
  }
  return ::open(path.c_str(), O_RDONLY);
#endif
}

void close_file(int fd) {
#ifdef _WIN32
  _close(fd);
#else
  ::close(fd);
#endif
}

bool BgzfFile::open(int open_fd, const std::string & filepath) {
  close();
  // Checking the first bytes reads them at an offset, which leaves the file's position at
  // the start. A pipe always fails the check, since a pipe can't be read at an offset.
  unsigned char head[HEADER_SIZE];
  if (read_at(open_fd, head, sizeof(head), 0) != static_cast<long long>(sizeof(head)) ||
      !is_bgzf_header(head)) {
    return false;
  }
  // A window of 15 + 16 makes inflate expect a gzip header, and check each block's CRC32 and size.
  if (inflateInit2(&strm, 15 + 16) != Z_OK) {
    return false;
  }
  fd = open_fd;
  path = filepath;
  raw.resize(BGZF_MAX_BLOCK_SIZE);
  return true;
}

void BgzfFile::close() {
  if (fd >= 0) {
    inflateEnd(&strm);
    close_file(fd);
    fd = -1;
  }
  std::vector<unsigned char>().swap(raw);
  raw_pos = 0;
  raw_len = 0;
  file_pos = 0;
}

// Read more of the file, until `raw` holds `count` compressed bytes from `raw_pos` on.
// Returns false if the file ends first. `raw` has room for the largest BGZF block.
bool BgzfFile::load(std::size_t count) {
  if (raw_len - raw_pos >= count) {
    return true;
  }
  // Move the bytes not yet used to the front, to make room to read more after them.
  raw_len -= raw_pos;
  std::memmove(raw.data(), raw.data() + raw_pos, raw_len);
  raw_pos = 0;
  while (raw_len < count) {
    long long got = read_at(fd, raw.data() + raw_len, raw.size() - raw_len, file_pos);
    if (got < 0) {
      throw std::invalid_argument("cannot read " + path);
    }
    if (got == 0) {
      return false;
    }
    raw_len += static_cast<std::size_t>(got);
    file_pos += static_cast<std::uint64_t>(got);
  }
  return true;
}

bool BgzfFile::next_block() {
  bool has_header = load(HEADER_SIZE);
  if (!has_header && raw_pos == raw_len) {
    return false;
  }
  if (!has_header || !is_bgzf_header(&raw[raw_pos])) {
    throw std::invalid_argument("invalid BGZF block header in " + path);
  }
  // The header's 'BC' field holds the block's total size minus 1, and the footer ends
  // with the uncompressed size.
  compressed_size = (raw[raw_pos + 16] | (raw[raw_pos + 17] << 8)) + 1;
  if (compressed_size < HEADER_SIZE + FOOTER_SIZE) {
    throw std::invalid_argument("invalid BGZF block size in " + path);
  }
  if (!load(compressed_size)) {
    throw std::invalid_argument("truncated BGZF block in " + path);
  }
  std::uint32_t isize = 0;
  std::memcpy(&isize, &raw[raw_pos + compressed_size - 4], sizeof(isize));
  if (isize > BGZF_MAX_BLOCK_SIZE) {
    throw std::invalid_argument("invalid BGZF block size in " + path);
  }
  uncompressed_size = isize;
  return true;
}

// The whole block is decompressed in a single `inflate` call.
void BgzfFile::inflate_block(char * dst) {
  inflateReset(&strm);
  strm.next_in = &raw[raw_pos];
  strm.avail_in = static_cast<unsigned int>(compressed_size);
  strm.next_out = reinterpret_cast<unsigned char *>(dst);
  strm.avail_out = static_cast<unsigned int>(uncompressed_size);
  raw_pos += compressed_size;
  if (inflate(&strm, Z_FINISH) != Z_STREAM_END || strm.avail_in != 0) {
    throw std::invalid_argument(std::string(strm.msg ? strm.msg : "corrupt BGZF block") + " in " + path);
  }
}

// Drop the bytes read ahead, so the next read starts at `offset`.
void BgzfFile::seek(std::uint64_t offset) {
  raw_pos = 0;
  raw_len = 0;
  file_pos = offset;
}

} // namespace bcf
