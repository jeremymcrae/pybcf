#include "gzreader.h"

#include <algorithm>
#include <cstdint>
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

// 64 KB is the most a BGZF block can decompress to. For other files it also amortizes
// `gzread` calls across many small 32-byte variant headers and site payloads.
const std::size_t BUFFER_SIZE = 65536;

// Every BGZF block starts with an 18-byte gzip header and ends with its CRC32 and
// uncompressed size (4 bytes each).
const std::size_t BGZF_HEADER_SIZE = 18;
const std::size_t BGZF_FOOTER_SIZE = 8;

int open_file(const std::string & path) {
#ifdef _WIN32
  return _open(path.c_str(), _O_RDONLY | _O_BINARY);
#else
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

// Open BGZF files for reading a block at a time, and anything else (plain gzip,
// uncompressed, or a pipe) through `gzread`. The reader stays closed if the file can't be opened.
void GzReader::open(const std::string & filepath) {
  close();
  path = filepath;
  int fd = open_file(path);
  if (fd < 0) {
    return;
  }
  // Checking for BGZF leaves the file's position at the start, ready for `gzread` if the
  // check fails. A pipe always fails it, since a pipe can't be read at an offset.
  unsigned char head[BGZF_HEADER_SIZE];
  bool is_bgzf = read_at(fd, head, sizeof(head), 0) == static_cast<long long>(sizeof(head)) &&
                 is_bgzf_header(head);
  // A window of 15 + 16 makes inflate expect a gzip header, and check each block's CRC32 and size.
  if (is_bgzf && inflateInit2(&strm, 15 + 16) == Z_OK) {
    bgzf = fd;
    raw.resize(BUFFER_SIZE);
  } else {
    // `gzdopen` takes over the open file, and `gzclose_r` closes it.
    file = gzdopen(fd, "rb");
    if (file == nullptr) {
      close_file(fd);
      return;
    }
    gzbuffer(file, BUFFER_SIZE);
  }
  buffer.resize(BUFFER_SIZE);
}

void GzReader::close() {
  buf_pos = 0;
  buf_len = 0;
  raw_pos = 0;
  raw_len = 0;
  file_pos = 0;
  std::vector<char>().swap(buffer);
  std::vector<unsigned char>().swap(raw);
  if (bgzf >= 0) {
    inflateEnd(&strm);
    close_file(bgzf);
    bgzf = -1;
  }
  if (file != nullptr) {
    gzFile f = file;
    file = nullptr;
    gzclose_r(f);
  }
}

// Read more of the file, until `raw` holds `count` compressed bytes from `raw_pos` on.
// Returns false if the file ends first. `raw` has room for the largest BGZF block.
bool GzReader::load(std::size_t count) {
  if (raw_len - raw_pos >= count) {
    return true;
  }
  // Move the bytes not yet used to the front, to make room to read more after them.
  raw_len -= raw_pos;
  std::memmove(raw.data(), raw.data() + raw_pos, raw_len);
  raw_pos = 0;
  while (raw_len < count) {
    long long got = read_at(bgzf, raw.data() + raw_len, raw.size() - raw_len, file_pos);
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

// Find the next BGZF block's compressed bytes in `raw`, returning false at the end of the file.
bool GzReader::next_block() {
  bool has_header = load(BGZF_HEADER_SIZE);
  if (!has_header && raw_pos == raw_len) {
    return false;
  }
  if (!has_header || !is_bgzf_header(&raw[raw_pos])) {
    throw std::invalid_argument("invalid BGZF block header in " + path);
  }
  // The header's 'BC' field holds the block's total size minus 1, and the footer ends
  // with the uncompressed size.
  compressed_size = (raw[raw_pos + 16] | (raw[raw_pos + 17] << 8)) + 1;
  if (compressed_size < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE) {
    throw std::invalid_argument("invalid BGZF block size in " + path);
  }
  if (!load(compressed_size)) {
    throw std::invalid_argument("truncated BGZF block in " + path);
  }
  std::uint32_t isize = 0;
  std::memcpy(&isize, &raw[raw_pos + compressed_size - 4], sizeof(isize));
  if (isize > BUFFER_SIZE) {
    throw std::invalid_argument("invalid BGZF block size in " + path);
  }
  block_size = isize;
  return true;
}

// Decompress the next block into `dst`, which must hold `block_size` bytes, in a single
// `inflate` call.
void GzReader::inflate_block(char * dst) {
  inflateReset(&strm);
  strm.next_in = &raw[raw_pos];
  strm.avail_in = static_cast<unsigned int>(compressed_size);
  strm.next_out = reinterpret_cast<unsigned char *>(dst);
  strm.avail_out = static_cast<unsigned int>(block_size);
  raw_pos += compressed_size;
  if (inflate(&strm, Z_FINISH) != Z_STREAM_END || strm.avail_in != 0) {
    throw std::invalid_argument(std::string(strm.msg ? strm.msg : "corrupt BGZF block") + " in " + path);
  }
}

// Decompress up to `count` bytes into `dst`, returning 0 once the end of the file is reached.
std::size_t GzReader::read_raw(char * dst, std::size_t count) {
  if (file == nullptr) {
    throw std::invalid_argument("cannot read from closed file: " + path);
  }
  int got = gzread(file, dst, static_cast<unsigned int>(count));
  if (got > 0) {
    return static_cast<std::size_t>(got);
  }
  // Nothing was read, which is an error unless the file ended cleanly.
  int errnum = 0;
  const char * msg = gzerror(file, &errnum);
  if (got < 0 || (errnum != Z_OK && errnum != Z_STREAM_END)) {
    // zlib's message starts with its own name for the file, so swap that for the path.
    const char * detail = std::strstr(msg, ": ");
    throw std::invalid_argument(path + ": " + (detail ? detail + 2 : msg));
  }
  return 0;
}

// Refill `buffer` with the next bytes of the file, returning false if there are none.
// A BGZF block can be empty (such as the end-of-file marker), leaving `buffer` empty.
bool GzReader::fill() {
  buf_pos = 0;
  buf_len = 0;
  if (bgzf >= 0) {
    if (!next_block()) {
      return false;
    }
    inflate_block(buffer.data());
    buf_len = block_size;
    return true;
  }
  buf_len = read_raw(buffer.data(), BUFFER_SIZE);
  return buf_len > 0;
}

// Serve small reads out of `buffer` with `memcpy`, while reading large payloads (such as
// multi-thousand-sample blocks) straight into `dst`: whole BGZF blocks that fit, or
// `gzread` calls of at least `BUFFER_SIZE` for other files.
std::size_t GzReader::read(void * dst, std::size_t count) {
  char * out = static_cast<char *>(dst);
  std::size_t total = 0;
  while (total < count) {
    std::size_t remaining = count - total;
    if (buf_pos < buf_len) {
      std::size_t take = std::min(remaining, buf_len - buf_pos);
      std::memcpy(out + total, buffer.data() + buf_pos, take);
      buf_pos += take;
      total += take;
    } else if (bgzf >= 0) {
      if (!next_block()) {
        break;
      }
      if (block_size <= remaining) {
        inflate_block(out + total);
        total += block_size;
      } else {
        inflate_block(buffer.data());
        buf_pos = 0;
        buf_len = block_size;
      }
    } else if (remaining >= BUFFER_SIZE) {
      // `gzread` takes a 32-bit count, so read at most 1 GB per call.
      std::size_t got = read_raw(out + total, std::min<std::size_t>(remaining, 1u << 30));
      if (got == 0) {
        break;
      }
      total += got;
    } else if (!fill()) {
      break;
    }
  }
  return total;
}

void GzReader::read_exact(void * dst, std::size_t count, const char * err_msg) {
  if (read(dst, count) != count) {
    throw std::invalid_argument(err_msg);
  }
}

// Advance past `count` uncompressed bytes without allocating or copying into a caller buffer.
void GzReader::skip(std::size_t count, const char * err_msg) {
  while (count > 0) {
    if (buf_pos == buf_len && !fill()) {
      throw std::invalid_argument(err_msg);
    }
    std::size_t take = std::min(count, buf_len - buf_pos);
    buf_pos += take;
    count -= take;
  }
}

// A BGZF virtual offset points at a compressed block, and a position in its uncompressed
// bytes, so read from that block next, then skip into it.
void GzReader::seek(std::uint64_t virtual_offset) {
  if (bgzf < 0) {
    throw std::invalid_argument("can only seek within a bgzipped file: " + path);
  }
  std::size_t within_block = virtual_offset & 0xffff;
  buf_pos = 0;
  buf_len = 0;
  raw_pos = 0;
  raw_len = 0;
  file_pos = virtual_offset >> 16;
  if (within_block > 0 && (!fill() || within_block > buf_len)) {
    throw std::invalid_argument("cannot seek within uncompressed bgzf block");
  }
  buf_pos = within_block;
}

} // namespace bcf
