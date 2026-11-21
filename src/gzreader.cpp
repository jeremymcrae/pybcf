#include "gzreader.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace bcf {

namespace {

// 64 KB matches the uncompressed size of a single BGZF block, amortizing `gzread` calls
// across many small 32-byte variant headers and site payloads.
const std::size_t BUFFER_SIZE = 65536;

int open_readonly(const char * filepath) {
#ifdef _WIN32
  return _open(filepath, _O_RDONLY | _O_BINARY);
#else
  return ::open(filepath, O_RDONLY);
#endif
}

bool seek_fd(int fd, std::uint64_t pos) {
#ifdef _WIN32
  return _lseeki64(fd, static_cast<__int64>(pos), SEEK_SET) >= 0;
#else
  return ::lseek(fd, static_cast<off_t>(pos), SEEK_SET) >= 0;
#endif
}

void close_fd(int fd) {
#ifdef _WIN32
  _close(fd);
#else
  ::close(fd);
#endif
}

} // namespace

void GzReader::open(const std::string & filepath) {
  close();
  path = filepath;
  file = gzopen(path.c_str(), "rb");
  if (file != nullptr) {
    buffer.resize(BUFFER_SIZE);
    gzbuffer(file, BUFFER_SIZE);
  }
}

void GzReader::close() {
  buf_pos = 0;
  buf_len = 0;
  std::vector<char>().swap(buffer);
  if (file != nullptr) {
    gzFile f = file;
    file = nullptr;
    gzclose_r(f);
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
bool GzReader::fill() {
  buf_pos = 0;
  buf_len = 0;
  buf_len = read_raw(buffer.data(), BUFFER_SIZE);
  return buf_len > 0;
}

// Serve small reads out of `buffer` with `memcpy`, while reading large payloads
// (`>= BUFFER_SIZE`, such as multi-thousand-sample blocks) straight into `dst`.
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

// BGZF files are concatenations of independent gzip members:
// 1. Seek the underlying OS file descriptor directly to the compressed block start (`c_offset`).
// 2. Wrap that positioned descriptor in a fresh `gzFile` via `gzdopen` (which takes ownership of `fd`).
// 3. Advance `u_offset` uncompressed bytes inside that BGZF block to reach the target record.
void GzReader::seek(Offsets offset) {
  close();

  int fd = open_readonly(path.c_str());
  if (fd < 0 || !seek_fd(fd, offset.c_offset)) {
    if (fd >= 0) {
      close_fd(fd);
    }
    throw std::invalid_argument("cannot seek within this gzfile");
  }

  file = gzdopen(fd, "rb");
  if (file == nullptr) {
    close_fd(fd);
    throw std::invalid_argument("cannot seek within this gzfile");
  }
  buffer.resize(BUFFER_SIZE);
  gzbuffer(file, BUFFER_SIZE);

  if (offset.u_offset > 0) {
    if (gzseek(file, static_cast<z_off_t>(offset.u_offset), SEEK_SET) < 0) {
      close();
      throw std::invalid_argument("cannot seek within uncompressed bgzf block");
    }
  }
}

} // namespace bcf
