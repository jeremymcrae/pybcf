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

void GzReader::check_error(int status) {
  int errnum = 0;
  const char * msg = gzerror(file, &errnum);
  if (status < 0 || (errnum != Z_OK && errnum != Z_STREAM_END)) {
    throw std::invalid_argument(msg && *msg ? msg : "error reading compressed file: " + path);
  }
}

std::size_t GzReader::read(void * dst, std::size_t count) {
  if (file == nullptr) {
    throw std::invalid_argument("cannot read from closed file: " + path);
  }
  char * out = static_cast<char *>(dst);
  std::size_t total = 0;

  while (total < count) {
    std::size_t avail = buf_len - buf_pos;
    if (avail > 0) {
      std::size_t take = std::min(count - total, avail);
      std::memcpy(out + total, buffer.data() + buf_pos, take);
      buf_pos += take;
      total += take;
      continue;
    }

    std::size_t rem = count - total;
    if (rem >= BUFFER_SIZE) {
      unsigned int chunk = static_cast<unsigned int>(std::min<std::size_t>(rem, 1u << 30));
      int got = gzread(file, out + total, chunk);
      if (got <= 0) {
        check_error(got);
        break;
      }
      total += static_cast<std::size_t>(got);
    } else {
      int got = gzread(file, buffer.data(), BUFFER_SIZE);
      buf_pos = 0;
      if (got <= 0) {
        buf_len = 0;
        check_error(got);
        break;
      }
      buf_len = static_cast<std::size_t>(got);
    }
  }
  return total;
}

void GzReader::read_exact(void * dst, std::size_t count, const char * err_msg) {
  if (read(dst, count) != count) {
    throw std::invalid_argument(err_msg);
  }
}

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
