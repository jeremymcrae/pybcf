#include "gzreader.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace bcf {

namespace {

// 64 KB is the most a BGZF block can decompress to. For other files it also amortizes
// `gzread` calls across many small 32-byte variant headers and site payloads.
const std::size_t BUFFER_SIZE = 65536;

// Every BGZF block starts with an 18-byte gzip header and ends with its CRC32 and
// uncompressed size (4 bytes each).
const std::size_t BGZF_HEADER_SIZE = 18;
const std::size_t BGZF_FOOTER_SIZE = 8;

bool seek_file(std::FILE * f, std::uint64_t pos) {
#ifdef _WIN32
  return _fseeki64(f, static_cast<__int64>(pos), SEEK_SET) == 0;
#else
  return fseeko(f, static_cast<off_t>(pos), SEEK_SET) == 0;
#endif
}

// A BGZF block is a gzip member whose header has a 'BC' extra field holding the block's size.
bool is_bgzf_header(const unsigned char * h) {
  return h[0] == 0x1f && h[1] == 0x8b && h[2] == 8 && (h[3] & 4) && h[10] == 6 && h[11] == 0 &&
         h[12] == 'B' && h[13] == 'C' && h[14] == 2 && h[15] == 0;
}

} // namespace

// Open BGZF files for reading a block at a time, and anything else (plain gzip or
// uncompressed) through `gzread`. The reader stays closed if the file can't be opened.
void GzReader::open(const std::string & filepath) {
  close();
  path = filepath;
  std::FILE * f = std::fopen(path.c_str(), "rb");
  if (f == nullptr) {
    return;
  }
  // Each block takes two small `fread` calls, so buffer the compressed bytes 64 KB at a time.
  std::setvbuf(f, nullptr, _IOFBF, BUFFER_SIZE);
  unsigned char head[BGZF_HEADER_SIZE];
  bool is_bgzf = std::fread(head, 1, sizeof(head), f) == sizeof(head) && is_bgzf_header(head);
  // A window of 15 + 16 makes inflate expect a gzip header, and check each block's CRC32 and size.
  if (is_bgzf && seek_file(f, 0) && inflateInit2(&strm, 15 + 16) == Z_OK) {
    bgzf = f;
  } else {
    std::fclose(f);
    file = gzopen(path.c_str(), "rb");
    if (file == nullptr) {
      return;
    }
    gzbuffer(file, BUFFER_SIZE);
  }
  buffer.resize(BUFFER_SIZE);
}

void GzReader::close() {
  buf_pos = 0;
  buf_len = 0;
  std::vector<char>().swap(buffer);
  std::vector<unsigned char>().swap(block);
  if (bgzf != nullptr) {
    inflateEnd(&strm);
    std::fclose(bgzf);
    bgzf = nullptr;
  }
  if (file != nullptr) {
    gzFile f = file;
    file = nullptr;
    gzclose_r(f);
  }
}

// Read the next BGZF block's compressed bytes into `block`, returning false at the end of the file.
bool GzReader::next_block() {
  block.resize(BGZF_HEADER_SIZE);
  std::size_t got = std::fread(block.data(), 1, BGZF_HEADER_SIZE, bgzf);
  if (got == 0 && !std::ferror(bgzf)) {
    return false;
  }
  if (got != BGZF_HEADER_SIZE || !is_bgzf_header(block.data())) {
    throw std::invalid_argument("invalid BGZF block header in " + path);
  }
  // The header's 'BC' field holds the block's total size minus 1, and the footer ends
  // with the uncompressed size.
  std::size_t total = (block[16] | (block[17] << 8)) + 1;
  if (total < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE) {
    throw std::invalid_argument("invalid BGZF block size in " + path);
  }
  block.resize(total);
  std::size_t rest = total - BGZF_HEADER_SIZE;
  if (std::fread(block.data() + BGZF_HEADER_SIZE, 1, rest, bgzf) != rest) {
    throw std::invalid_argument("truncated BGZF block in " + path);
  }
  std::uint32_t isize = 0;
  std::memcpy(&isize, &block[total - 4], sizeof(isize));
  if (isize > BUFFER_SIZE) {
    throw std::invalid_argument("invalid BGZF block size in " + path);
  }
  block_size = isize;
  return true;
}

// Decompress `block` into `dst`, which must hold `block_size` bytes, in a single `inflate` call.
void GzReader::inflate_block(char * dst) {
  inflateReset(&strm);
  strm.next_in = block.data();
  strm.avail_in = static_cast<unsigned int>(block.size());
  strm.next_out = reinterpret_cast<unsigned char *>(dst);
  strm.avail_out = static_cast<unsigned int>(block_size);
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
    throw std::invalid_argument(msg && *msg ? msg : "error reading compressed file: " + path);
  }
  return 0;
}

// Refill `buffer` with the next bytes of the file, returning false if there are none.
// A BGZF block can be empty (such as the end-of-file marker), leaving `buffer` empty.
bool GzReader::fill() {
  buf_pos = 0;
  buf_len = 0;
  if (bgzf != nullptr) {
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
    } else if (bgzf != nullptr) {
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

// A BGZF virtual offset points at a compressed block (`c_offset`) and a position in its
// uncompressed bytes (`u_offset`), so move the file to that block, then skip into it.
void GzReader::seek(Offsets offset) {
  if (bgzf == nullptr) {
    throw std::invalid_argument("can only seek within a bgzipped file: " + path);
  }
  buf_pos = 0;
  buf_len = 0;
  if (!seek_file(bgzf, offset.c_offset)) {
    throw std::invalid_argument("cannot seek within this gzfile");
  }
  if (offset.u_offset > 0 && (!fill() || offset.u_offset > buf_len)) {
    throw std::invalid_argument("cannot seek within uncompressed bgzf block");
  }
  buf_pos = static_cast<std::size_t>(offset.u_offset);
}

} // namespace bcf
