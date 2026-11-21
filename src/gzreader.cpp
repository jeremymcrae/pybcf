#include "gzreader.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace bcf {

namespace {

// 64 KB holds any BGZF block once decompressed. For other files it also amortizes
// `gzread` calls across many small 32-byte variant headers and site payloads.
const std::size_t BUFFER_SIZE = BGZF_MAX_BLOCK_SIZE;

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
  if (!bgzf.open(fd, path)) {
    // `gzdopen` takes over the open file, and `gzclose_r` closes it.
    gz = gzdopen(fd, "rb");
    if (gz == nullptr) {
      close_file(fd);
      return;
    }
    gzbuffer(gz, BUFFER_SIZE);
  }
  buffer.resize(BUFFER_SIZE);
}

void GzReader::close() {
  buf_pos = 0;
  buf_len = 0;
  std::vector<char>().swap(buffer);
  bgzf.close();
  if (gz != nullptr) {
    gzFile f = gz;
    gz = nullptr;
    gzclose_r(f);
  }
}

// Decompress up to `count` bytes into `dst` with `gzread`, returning 0 once the end of
// the file is reached.
std::size_t GzReader::read_gz(char * dst, std::size_t count) {
  if (gz == nullptr) {
    throw std::invalid_argument("cannot read from closed file: " + path);
  }
  int got = gzread(gz, dst, static_cast<unsigned int>(count));
  if (got > 0) {
    return static_cast<std::size_t>(got);
  }
  // Nothing was read, which is an error unless the file ended cleanly.
  int errnum = 0;
  const char * msg = gzerror(gz, &errnum);
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
  if (bgzf.is_open()) {
    if (!bgzf.next_block()) {
      return false;
    }
    bgzf.inflate_block(buffer.data());
    buf_len = bgzf.block_size();
    return true;
  }
  buf_len = read_gz(buffer.data(), BUFFER_SIZE);
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
    } else if (bgzf.is_open()) {
      if (!bgzf.next_block()) {
        break;
      }
      if (bgzf.block_size() <= remaining) {
        bgzf.inflate_block(out + total);
        total += bgzf.block_size();
      } else {
        bgzf.inflate_block(buffer.data());
        buf_pos = 0;
        buf_len = bgzf.block_size();
      }
    } else if (remaining >= BUFFER_SIZE) {
      // `gzread` takes a 32-bit count, so read at most 1 GB per call.
      std::size_t got = read_gz(out + total, std::min<std::size_t>(remaining, 1u << 30));
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
  if (!bgzf.is_open()) {
    throw std::invalid_argument("can only seek within a bgzipped file: " + path);
  }
  std::size_t within_block = virtual_offset & 0xffff;
  buf_pos = 0;
  buf_len = 0;
  bgzf.seek(virtual_offset >> 16);
  if (within_block > 0 && (!fill() || within_block > buf_len)) {
    throw std::invalid_argument("cannot seek within uncompressed bgzf block");
  }
  buf_pos = within_block;
}

} // namespace bcf
