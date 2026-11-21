#include <cmath>
#include <cstring>
#include <stdexcept>

#include "sample_data.h"
#include "simd.h"
#include "types.h"

namespace bcf {

// Read the raw per-sample byte block from the BCF stream into an uninitialized buffer.
// Field descriptors and sample values are only parsed later if `samples` is accessed.
SampleData::SampleData(GzReader & infile, const Header & _header, std::uint32_t len,
                       std::uint32_t _n_fmt, std::uint32_t _n_samples)
    : header(&_header), n_fmt(_n_fmt), buf_len(len), n_samples(_n_samples) {
  if (len > 0) {
    buf.reset(new char[len]);
    infile.read_exact(buf.get(), len, "truncated BCF sample data");
  }
}

// Walk the sample buffer once to locate each FORMAT field.
// Unlike VCF (which groups values by sample), BCF stores FORMAT data column-wise:
// each field has a header dictionary key ID, a per-sample type/count descriptor,
// and then `n_samples * n_vals * type_size` contiguous value bytes.
void SampleData::parse() {
  if (is_parsed) {
    return;
  }
  if (n_fmt == 0) {
    is_parsed = true;
    return;
  }
  fields.reserve(n_fmt);
  Cursor cursor(buf.get(), buf_len);
  for (std::uint32_t i = 0; i < n_fmt; i++) {
    // Read the FORMAT tag's dictionary ID and look up its definition in the header.
    std::int32_t format_idx = cursor.read_typed_int();
    if (format_idx < 0 || static_cast<std::size_t>(format_idx) >= header->format.size() ||
        header->format[format_idx].id.empty()) {
      throw std::invalid_argument("unknown FORMAT index: " + std::to_string(format_idx));
    }
    std::uint32_t uid = static_cast<std::uint32_t>(format_idx);
    const FormatField * fmt_field = &header->format[uid];
    bool is_gt = fmt_field->is_gt;

    // Read the per-sample element type and per-sample vector length (`n_vals`).
    TypeDescriptor desc(cursor);
    if (desc.type == flag && desc.n_vals != 0) {
      throw std::invalid_argument("invalid BCF FORMAT flag length");
    }
    if (is_gt && !is_int_type(desc.type)) {
      throw std::invalid_argument("GT field must be encoded as integer");
    }

    // Record where the values start, then step over the contiguous block of them for
    // all samples, which must fit inside `buf`, to reach the next FORMAT field.
    FormatLayout layout = {desc.type, desc.type_size, cursor.offset(), desc.n_vals, is_gt};
    std::uint64_t field_bytes = static_cast<std::uint64_t>(desc.n_vals) *
                                static_cast<std::uint64_t>(desc.type_size) *
                                static_cast<std::uint64_t>(n_samples);
    cursor.skip(field_bytes, "FORMAT field exceeds BCF sample buffer");
    fields.push_back({uid, fmt_field, layout});
  }
  is_parsed = true;
}

std::size_t SampleData::size() {
  parse();
  return fields.size();
}

std::vector<NamedFormat> & SampleData::get_fields() {
  parse();
  return fields;
}

// Look up a FORMAT field by name without allocating a temporary std::string. If a malformed
// record repeats a key, the first copy wins, matching htslib.
bool SampleData::find(const char * key, std::size_t key_len, FormatLayout & out) {
  parse();
  for (const auto & entry : fields) {
    const std::string & id = entry.field->id;
    if (id.size() == key_len && std::memcmp(id.data(), key, key_len) == 0) {
      out = entry.layout;
      return true;
    }
  }
  return false;
}

// Decode an integer FORMAT field (including GT) into `dst` as float64, so missing and
// padding values can be NaN. GT alleles become 0-based allele indices.
void SampleData::read_ints(const FormatLayout & layout, double * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(layout.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(layout.offset) + max_n * layout.type_size > buf_len) {
    throw std::invalid_argument("FORMAT integer data exceeds buffer");
  }
  std::uint64_t n = 0;

#ifdef BCF_HAVE_SIMD
  if (simd_supported()) {
    const char * src = buf.get() + layout.offset;
    if (layout.is_gt) {
      // Only 8-bit GT (fewer than 63 ALT alleles) has a SIMD kernel.
      n = (layout.type_size == 1) ? simd_gt(src, dst, max_n) : 0;
    } else if (layout.type_size == 1) {
      n = simd_int8(src, dst, max_n);
    } else if (layout.type_size == 2) {
      n = simd_int16(src, dst, max_n);
    } else {
      n = simd_int32(src, dst, max_n);
    }
  }
#endif

  Cursor cursor(buf.get(), buf_len,
                layout.offset + static_cast<std::uint32_t>(n * layout.type_size));
  const double qnan = std::nan("");
  for (; n < max_n; n++) {
    std::int32_t v = cursor.read_int(layout.type_size);
    if (layout.is_gt) {
      v = (v < 0) ? -1 : (v >> 1) - 1;
      dst[n] = (v < 0) ? qnan : static_cast<double>(v);
    } else {
      dst[n] = (v == INT32_MIN) ? qnan : static_cast<double>(v);
    }
  }
}

// Compute per-sample phasing flags (0 = unphased, 1 = phased) from GT without decoding
// allele numbers. Uses a SIMD fast path for standard diploid 8-bit GT and a general
// scalar loop for haploid, polyploid, 16/32-bit GT, and tail samples.
void SampleData::read_phase(const FormatLayout & layout, std::uint8_t * dst) const {
  if (n_samples == 0) {
    return;
  }
  if (layout.n_vals == 0) {
    std::memset(dst, 0, n_samples);
    return;
  }
  std::uint64_t max_n = static_cast<std::uint64_t>(layout.n_vals) * n_samples;
  if (static_cast<std::uint64_t>(layout.offset) + max_n * layout.type_size > buf_len) {
    throw std::invalid_argument("FORMAT GT data exceeds buffer");
  }
  const char * raw = buf.get();
  std::uint32_t s = 0;

#ifdef BCF_HAVE_SIMD
  if (simd_supported() && layout.n_vals == 2 && layout.type_size == 1) {
    s = simd_phase(raw + layout.offset, dst, n_samples);
  }
#endif

  // In BCF, allele 0 carries no phase relative to a prior allele; a sample is phased if
  // every non-padding allele after the first (`i > 0`, `v >= 0`) has its low phase bit set.
  Cursor cursor(raw, buf_len, layout.offset + s * layout.n_vals * layout.type_size);
  for (; s < n_samples; s++) {
    std::uint8_t is_phased = 1;
    for (std::uint32_t i = 0; i < layout.n_vals; i++) {
      std::int32_t v = cursor.read_int(layout.type_size);
      if (v >= 0 && i > 0) {
        is_phased &= static_cast<std::uint8_t>(v & 1);
      }
    }
    dst[s] = is_phased;
  }
}

// Decode 32-bit IEEE float FORMAT values into `dst`, replacing BCF's signaling-NaN
// missing/padding sentinels (`0x7f800001..0x7f800007`) with quiet NaN.
void SampleData::read_floats(const FormatLayout & layout, float * dst) const {
  std::uint64_t max_n = static_cast<std::uint64_t>(layout.n_vals) * n_samples;
  if (max_n == 0) {
    return;
  }
  if (static_cast<std::uint64_t>(layout.offset) + max_n * sizeof(float) > buf_len) {
    throw std::invalid_argument("FORMAT float data exceeds buffer");
  }
  const char * src = buf.get() + layout.offset;
  std::uint64_t n = 0;

#ifdef BCF_HAVE_SIMD
  if (simd_supported()) {
    n = simd_floats(src, dst, max_n);
  }
#endif

  const float qnan = std::nanf("");
  for (; n < max_n; n++) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, src + n * sizeof(float), sizeof(bits));
    if (is_bcf_float_missing(bits)) {
      dst[n] = qnan;
    } else {
      std::memcpy(dst + n, &bits, sizeof(float));
    }
  }
}

// Extract per-sample strings: BCF stores FORMAT strings as fixed `n_vals`-byte slices
// per sample, padded with NUL ('\0') bytes for shorter or missing strings.
std::vector<std::string> SampleData::get_strings(const FormatLayout & layout) const {
  std::vector<std::string> vals(n_samples);
  if (n_samples == 0 || layout.n_vals == 0) {
    return vals;
  }
  Cursor cursor(buf.get(), buf_len, layout.offset);
  for (std::uint32_t n = 0; n < n_samples; n++) {
    vals[n] = cursor.read_string(layout.n_vals);
  }
  return vals;
}

} // namespace bcf
