#include "gzreader.h"

#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace bcf {

void GzReader::open(const std::string & filepath) {
  close();
  path = filepath;
  file = gzopen(path.c_str(), "rb");
  if (file != nullptr) {
    gzbuffer(file, sizeof(buffer));
  }
}

void GzReader::close() {
  buf_pos = 0;
  buf_len = 0;
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
      std::memcpy(out + total, buffer + buf_pos, take);
      buf_pos += take;
      total += take;
      continue;
    }

    std::size_t rem = count - total;
    if (rem >= sizeof(buffer)) {
      unsigned int chunk = static_cast<unsigned int>(std::min<std::size_t>(rem, 1u << 30));
      int got = gzread(file, out + total, chunk);
      if (got <= 0) {
        check_error(got);
        break;
      }
      total += static_cast<std::size_t>(got);
    } else {
      int got = gzread(file, buffer, sizeof(buffer));
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

#ifdef _WIN32
  int fd = _open(path.c_str(), _O_RDONLY | _O_BINARY);
  if (fd < 0 || _lseeki64(fd, static_cast<__int64>(offset.c_offset), SEEK_SET) < 0) {
    if (fd >= 0) {
      _close(fd);
    }
    throw std::invalid_argument("cannot seek within this gzfile");
  }
#else
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0 || ::lseek(fd, static_cast<off_t>(offset.c_offset), SEEK_SET) < 0) {
    if (fd >= 0) {
      ::close(fd);
    }
    throw std::invalid_argument("cannot seek within this gzfile");
  }
#endif

  file = gzdopen(fd, "rb");
  if (file == nullptr) {
#ifdef _WIN32
    _close(fd);
#else
    ::close(fd);
#endif
    throw std::invalid_argument("cannot seek within this gzfile");
  }
  gzbuffer(file, sizeof(buffer));

  if (offset.u_offset > 0) {
    if (gzseek(file, static_cast<z_off_t>(offset.u_offset), SEEK_SET) < 0) {
      close();
      throw std::invalid_argument("cannot seek within uncompressed bgzf block");
    }
  }
}

} // namespace bcf
